import argparse
import csv
import json
import math
import os
import sqlite3
from contextlib import closing
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlsplit


ROOT = Path(__file__).resolve().parent
CSV_FILE = ROOT / "RCSB_MASTER_DATA_EXPANDED_ANNOTATED.csv"
DATABASE_FILE = ROOT / "data" / "pdb_database.sqlite"
NUMERICAL_FIELDS = {
    "FRET_Diffusion_Length_Angstroms",
    "DET_Diffusion_Length_Angstroms",
    "TaxID",
}


def clean_record(row):
    record = {}
    for key, raw_value in row.items():
        value = (raw_value or "").strip()
        if key not in NUMERICAL_FIELDS:
            record[key] = value
        elif not value:
            record[key] = None
        else:
            try:
                number = float(value)
            except ValueError:
                record[key] = value
            else:
                if not math.isfinite(number):
                    record[key] = value
                elif key == "TaxID":
                    record[key] = int(number)
                else:
                    record[key] = number
    return record


def database_is_current():
    if not DATABASE_FILE.exists():
        return False

    source_stat = CSV_FILE.stat()
    try:
        with closing(sqlite3.connect(DATABASE_FILE)) as connection:
            metadata = dict(
                connection.execute("SELECT key, value FROM metadata")
            )
    except sqlite3.DatabaseError:
        return False

    return (
        metadata.get("source_size") == str(source_stat.st_size)
        and metadata.get("source_mtime_ns") == str(source_stat.st_mtime_ns)
    )


def build_database():
    if not CSV_FILE.is_file():
        raise FileNotFoundError(f"Required source CSV was not found: {CSV_FILE}")

    DATABASE_FILE.parent.mkdir(parents=True, exist_ok=True)
    temporary_file = DATABASE_FILE.with_name(DATABASE_FILE.name + ".tmp")
    temporary_file.unlink(missing_ok=True)
    source_stat = CSV_FILE.stat()

    try:
        with closing(sqlite3.connect(temporary_file)) as connection, connection:
            connection.execute(
                "CREATE TABLE records ("
                "pdb_id TEXT PRIMARY KEY COLLATE NOCASE, "
                "record_json TEXT NOT NULL)"
            )
            connection.execute(
                "CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)"
            )

            with CSV_FILE.open(
                "r", encoding="utf-8-sig", newline=""
            ) as csv_file:
                reader = csv.DictReader(csv_file)
                if not reader.fieldnames or "PDB_ID" not in reader.fieldnames:
                    raise ValueError("Source CSV must include a PDB_ID column.")

                batch = []
                for row in reader:
                    pdb_id = (row.get("PDB_ID") or "").strip().upper()
                    if not pdb_id:
                        continue

                    record = clean_record(row)
                    batch.append(
                        (
                            pdb_id,
                            json.dumps(record, ensure_ascii=False, allow_nan=False),
                        )
                    )
                    if len(batch) >= 1000:
                        connection.executemany(
                            "INSERT INTO records (pdb_id, record_json) VALUES (?, ?) "
                            "ON CONFLICT(pdb_id) DO UPDATE SET "
                            "record_json = excluded.record_json",
                            batch,
                        )
                        batch.clear()

                if batch:
                    connection.executemany(
                        "INSERT INTO records (pdb_id, record_json) VALUES (?, ?) "
                        "ON CONFLICT(pdb_id) DO UPDATE SET "
                        "record_json = excluded.record_json",
                        batch,
                    )

            connection.executemany(
                "INSERT INTO metadata (key, value) VALUES (?, ?)",
                (
                    ("source_size", str(source_stat.st_size)),
                    ("source_mtime_ns", str(source_stat.st_mtime_ns)),
                ),
            )
            record_count = connection.execute(
                "SELECT COUNT(*) FROM records"
            ).fetchone()[0]

        os.replace(temporary_file, DATABASE_FILE)
    finally:
        temporary_file.unlink(missing_ok=True)

    print(f"Indexed {record_count:,} unique PDB records from {CSV_FILE.name}.")


def ensure_database():
    if database_is_current():
        with closing(sqlite3.connect(DATABASE_FILE)) as connection:
            count = connection.execute(
                "SELECT COUNT(*) FROM records"
            ).fetchone()[0]
        print(f"Using current index with {count:,} unique PDB records.")
        return
    print(f"Building the PDB index from {CSV_FILE.name}...")
    build_database()


class PDBRequestHandler(SimpleHTTPRequestHandler):
    def send_json(self, status, payload):
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        request = urlsplit(self.path)
        if request.path == "/api/status":
            with closing(sqlite3.connect(DATABASE_FILE)) as connection:
                count = connection.execute(
                    "SELECT COUNT(*) FROM records"
                ).fetchone()[0]
            self.send_json(200, {"ready": True, "recordCount": count})
            return

        if request.path == "/api/record":
            pdb_id = parse_qs(request.query).get("id", [""])[0].strip().upper()
            if len(pdb_id) != 4 or not pdb_id.isascii() or not pdb_id.isalnum():
                self.send_json(
                    400, {"error": "Enter a valid four-character PDB ID."}
                )
                return

            with closing(sqlite3.connect(DATABASE_FILE)) as connection:
                row = connection.execute(
                    "SELECT record_json FROM records WHERE pdb_id = ?",
                    (pdb_id,),
                ).fetchone()
            if row is None:
                self.send_json(404, {"error": f"No record found for {pdb_id}."})
                return
            self.send_json(200, json.loads(row[0]))
            return

        super().do_GET()


def main():
    parser = argparse.ArgumentParser(
        description="Serve the PDB Diffusion Database and its indexed lookup API."
    )
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    args = parser.parse_args()

    ensure_database()
    handler = partial(PDBRequestHandler, directory=str(ROOT))
    server = ThreadingHTTPServer((args.host, args.port), handler)
    print(f"Serving PDB Diffusion Database at http://{args.host}:{args.port}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down.")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
