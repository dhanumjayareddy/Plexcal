import argparse
import csv
from email import policy
from email.parser import BytesParser
import json
import math
import os
import re
import sqlite3
import tempfile
from contextlib import closing
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import parse_qs, quote, urlsplit
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parent
CSV_FILE = ROOT / "RCSB_MASTER_DATA_EXPANDED_ANNOTATED.csv"
DATABASE_FILE = ROOT / "data" / "pdb_database.sqlite"
MAX_PDB_UPLOAD_BYTES = 25 * 1024 * 1024
COHERENT_PAIR_METRICS = (
    "diffusion_length_angstrom",
    "total_pairs",
    "coherent_pairs",
    "coherent_pair_fraction",
)
NUMERICAL_FIELDS = {
    "FRET_Diffusion_Length_Angstroms",
    "DET_Diffusion_Length_Angstroms",
    "TaxID",
}


def _download_structure(source, identifier):
    identifier = identifier.strip().upper()
    if source == "rcsb":
        if not re.fullmatch(r"[A-Z0-9]{4}", identifier):
            raise ValueError("Enter a valid four-character RCSB PDB ID.")
        structure_url = f"https://files.rcsb.org/download/{identifier}.pdb"
        filename = f"{identifier}.pdb"
    elif source == "alphafold":
        accession_match = re.fullmatch(
            r"(?:AF-)?([A-Z][0-9][A-Z0-9]{3,8}[0-9])(?:-F1)?"
            r"(?:-MODEL_V[0-9]+)?(?:\.PDB)?",
            identifier,
        )
        if not accession_match:
            raise ValueError(
                "Enter a UniProt accession or AlphaFold ID such as AF-P05067-F1."
            )
        accession = accession_match.group(1)
        metadata_request = Request(
            f"https://alphafold.ebi.ac.uk/api/prediction/{quote(accession)}",
            headers={"User-Agent": "Plexcal/1.0"},
        )
        with urlopen(metadata_request, timeout=30) as response:
            metadata = json.loads(response.read(1024 * 1024))
        if (
            not isinstance(metadata, list)
            or not metadata
            or not isinstance(metadata[0], dict)
        ):
            raise ValueError(f"No AlphaFold structure found for {accession}.")
        structure_url = metadata[0].get("pdbUrl")
        parsed_structure_url = (
            urlsplit(structure_url) if isinstance(structure_url, str) else None
        )
        if (
            parsed_structure_url is None
            or parsed_structure_url.scheme != "https"
            or parsed_structure_url.hostname != "alphafold.ebi.ac.uk"
            or not parsed_structure_url.path.startswith("/files/")
            or not parsed_structure_url.path.lower().endswith(".pdb")
        ):
            raise ValueError(f"AlphaFold returned no valid PDB file for {accession}.")
        identifier = accession
        filename = f"AF-{accession}-F1.pdb"
    else:
        raise ValueError("Choose upload, RCSB PDB ID, or AlphaFold structure ID.")

    structure_request = Request(
        structure_url, headers={"User-Agent": "Plexcal/1.0"}
    )
    with urlopen(structure_request, timeout=60) as response:
        pdb_bytes = response.read(MAX_PDB_UPLOAD_BYTES + 1)
    if len(pdb_bytes) > MAX_PDB_UPLOAD_BYTES:
        raise ValueError("The downloaded PDB file exceeds the 25 MiB limit.")
    if not any(
        line.startswith((b"ATOM  ", b"HETATM")) for line in pdb_bytes.splitlines()
    ):
        raise ValueError(f"No PDB atom records were downloaded for {identifier}.")
    return pdb_bytes, identifier, filename


def _finite_metric(result, field):
    value = result.get(field)
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    if field in {"total_pairs", "coherent_pairs"} and number.is_integer():
        return int(number)
    return number


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
            self.send_json(
                200,
                {
                    "ready": True,
                    "recordCount": count,
                    "fretCalculationAvailable": True,
                    "detCalculationAvailable": True,
                },
            )
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

    def do_POST(self):
        request = urlsplit(self.path)
        if request.path not in {"/api/calculate", "/api/calculate-fret"}:
            self.send_json(404, {"error": "API endpoint not found."})
            return
        legacy_fret_endpoint = request.path == "/api/calculate-fret"

        try:
            content_length = int(self.headers.get("Content-Length", ""))
        except ValueError:
            self.send_json(411, {"error": "A valid Content-Length is required."})
            return

        if content_length <= 0:
            self.send_json(400, {"error": "Upload a non-empty PDB file."})
            return
        if content_length > MAX_PDB_UPLOAD_BYTES + 1024 * 1024:
            self.send_json(413, {"error": "The upload exceeds the 25 MiB limit."})
            return
        if not self.headers.get("Content-Type", "").lower().startswith(
            "multipart/form-data"
        ):
            self.send_json(
                415, {"error": "Upload the PDB file as multipart form data."}
            )
            return

        body = self.rfile.read(content_length)
        if len(body) != content_length:
            self.send_json(400, {"error": "The upload was incomplete."})
            return

        try:
            message = BytesParser(policy=policy.default).parsebytes(
                b"Content-Type: "
                + self.headers["Content-Type"].encode("ascii")
                + b"\r\nMIME-Version: 1.0\r\n\r\n"
                + body
            )
        except (UnicodeEncodeError, ValueError) as error:
            self.send_json(400, {"error": f"Could not read the upload: {error}"})
            return
        if not message.is_multipart():
            self.send_json(400, {"error": "The upload form is malformed."})
            return

        upload_part = None
        source = "upload"
        identifier = ""
        for part in message.iter_parts():
            field_name = part.get_param("name", header="content-disposition")
            if field_name == "pdbFile":
                upload_part = part
            elif field_name == "inputSource":
                source = (part.get_payload(decode=True) or b"").decode(
                    "ascii", errors="replace"
                ).strip().lower()
            elif field_name == "structureId":
                identifier = (part.get_payload(decode=True) or b"").decode(
                    "ascii", errors="replace"
                ).strip()

        if legacy_fret_endpoint and source != "upload":
            self.send_json(400, {"error": "The legacy FRET endpoint accepts uploads only."})
            return

        if source == "upload":
            if upload_part is None:
                self.send_json(400, {"error": "Choose a PDB file to calculate."})
                return

            pdb_bytes = upload_part.get_payload(decode=True) or b""
            if not pdb_bytes:
                self.send_json(400, {"error": "The selected file is empty."})
                return
            if len(pdb_bytes) > MAX_PDB_UPLOAD_BYTES:
                self.send_json(
                    413, {"error": "The PDB file must be 25 MiB or smaller."}
                )
                return
            if not any(
                line.startswith((b"ATOM  ", b"HETATM"))
                for line in pdb_bytes.splitlines()
            ):
                self.send_json(
                    400, {"error": "The uploaded file has no PDB atom records."}
                )
                return
            original_name = Path(
                upload_part.get_filename() or "uploaded.pdb"
            ).name
            identifier = original_name
        elif source in {"rcsb", "alphafold"}:
            if upload_part is not None:
                self.send_json(
                    400, {"error": "Provide either a PDB upload or a structure ID."}
                )
                return
            if not identifier:
                self.send_json(400, {"error": "Enter a structure ID to calculate."})
                return
            try:
                pdb_bytes, identifier, original_name = _download_structure(
                    source, identifier
                )
            except ValueError as error:
                self.send_json(400, {"error": str(error)})
                return
            except (HTTPError, URLError, TimeoutError, OSError, json.JSONDecodeError) as error:
                self.send_json(
                    502, {"error": f"Could not download the requested PDB structure: {error}"}
                )
                return
        else:
            self.send_json(400, {"error": "Choose upload, RCSB, or AlphaFold input."})
            return

        temporary_path = None
        try:
            with tempfile.NamedTemporaryFile(
                prefix="pdb-calculation-", suffix=".pdb", delete=False
            ) as temporary_file:
                temporary_file.write(pdb_bytes)
                temporary_path = Path(temporary_file.name)

            from process_single_pdb import process_single_pdb

            fret_result = process_single_pdb(temporary_path)
            if not isinstance(fret_result, dict):
                raise TypeError("The FRET processor must return a metric dictionary.")

            if legacy_fret_endpoint:
                fret_length = _finite_metric(fret_result, "diffusion_length_angstrom")
                if fret_length is None:
                    self.send_json(
                        422,
                        {
                            "error": (
                                "The FRET calculation failed for this PDB. Check "
                                "that aromatic residue ring atoms and coordinates are present."
                            )
                        },
                    )
                else:
                    self.send_json(
                        200,
                        {
                            "fileName": original_name,
                            "diffusionLengthAngstroms": fret_length,
                        },
                    )
                return

            from det_code_clean import process_single_pdb as process_det_pdb

            det_result = process_det_pdb(temporary_path)
            if not isinstance(det_result, dict):
                raise TypeError("The DET processor must return a metric dictionary.")

            def metrics_for(result, method):
                metrics = {
                    field: _finite_metric(result, field)
                    for field in COHERENT_PAIR_METRICS
                }
                error = None
                if metrics["diffusion_length_angstrom"] is None:
                    error = (
                        "Check that aromatic residue ring atoms and coordinates "
                        "are present."
                        if method == "fret"
                        else "Check that the structure contains valid TRP, TYR, "
                        "or PHE residue coordinates."
                    )
                if error:
                    metrics["error"] = error
                return metrics

            fret_metrics = metrics_for(fret_result, "fret")
            det_metrics = metrics_for(det_result, "det")
            response = {
                "input": {
                    "source": source,
                    "identifier": identifier,
                    "filename": original_name,
                },
                "fret": fret_metrics,
                "det": det_metrics,
            }
            failed = (
                fret_metrics["diffusion_length_angstrom"] is None
                or det_metrics["diffusion_length_angstrom"] is None
            )
            self.send_json(422 if failed else 200, response)
        except (OSError, ImportError, TypeError, ValueError) as error:
            self.send_json(
                500, {"error": f"Could not run the requested calculation: {error}"}
            )
        finally:
            if temporary_path is not None:
                temporary_path.unlink(missing_ok=True)


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
