import csv
import json
from pathlib import Path

from server import CSV_FILE, ROOT, clean_record


PAGES_DIR = ROOT / "docs"
SHARDS_DIR = PAGES_DIR / "data" / "shards"
MANIFEST_FILE = PAGES_DIR / "data" / "index.json"
ASSETS = {
    ROOT / "index.html": PAGES_DIR / "index.html",
    ROOT / "style.css": PAGES_DIR / "style.css",
    ROOT / "github_pages.js": PAGES_DIR / "script.js",
}
SHARD_TARGET_BYTES = 900_000


def write_json(path, payload):
    temporary_path = path.with_name(path.name + ".tmp")
    with temporary_path.open("w", encoding="utf-8", newline="\n") as output:
        json.dump(
            payload,
            output,
            ensure_ascii=False,
            separators=(",", ":"),
            allow_nan=False,
        )
    temporary_path.replace(path)


def main():
    if not CSV_FILE.is_file():
        raise FileNotFoundError(f"Required source CSV was not found: {CSV_FILE}")

    PAGES_DIR.mkdir(parents=True, exist_ok=True)
    SHARDS_DIR.mkdir(parents=True, exist_ok=True)
    for source, target in ASSETS.items():
        if not source.is_file():
            raise FileNotFoundError(f"Required static site file was not found: {source}")
        target.write_bytes(source.read_bytes())

    shards = {}
    with CSV_FILE.open("r", encoding="utf-8-sig", newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        if not reader.fieldnames or "PDB_ID" not in reader.fieldnames:
            raise ValueError("Source CSV must include a PDB_ID column.")

        for row in reader:
            pdb_id = (row.get("PDB_ID") or "").strip().upper()
            if len(pdb_id) != 4 or not pdb_id.isascii() or not pdb_id.isalnum():
                continue
            prefix = pdb_id[:2].lower()
            shards.setdefault(prefix, {})[pdb_id] = clean_record(row)

    total_bytes = 0
    for prefix, records in sorted(shards.items()):
        shard_file = SHARDS_DIR / f"{prefix}.json"
        write_json(shard_file, records)
        shard_size = shard_file.stat().st_size
        total_bytes += shard_size
        if shard_size > SHARD_TARGET_BYTES:
            raise ValueError(
                f"Shard {shard_file.name} is {shard_size:,} bytes; "
                f"expected at most {SHARD_TARGET_BYTES:,} bytes."
            )

    expected_shards = {f"{prefix}.json" for prefix in shards}
    for old_shard in SHARDS_DIR.glob("*.json"):
        if old_shard.name not in expected_shards:
            old_shard.unlink()

    write_json(
        MANIFEST_FILE,
        {
            "recordCount": sum(map(len, shards.values())),
            "shards": {prefix: len(records) for prefix, records in sorted(shards.items())},
        },
    )
    (PAGES_DIR / ".nojekyll").touch()
    print(
        f"Built {PAGES_DIR.relative_to(ROOT)} with "
        f"{sum(map(len, shards.values())):,} unique records in {len(shards)} shards."
    )
    print(f"Total shard size: {total_bytes / (1024 * 1024):.2f} MiB.")
    print(f"Largest shard: {max((path.stat().st_size for path in SHARDS_DIR.glob('*.json')), default=0):,} bytes.")


if __name__ == "__main__":
    main()
