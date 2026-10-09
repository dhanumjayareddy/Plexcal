"""Run the existing single-PDB calculation across a folder in parallel."""

from __future__ import annotations

import argparse
import csv
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from process_single_pdb import process_single_pdb


PDB_DIRECTORY = Path(__file__).resolve().parent / "all_proteins_pdb"
DEFAULT_WORKERS = max(1, min(4, os.cpu_count() or 1))


def _process_one_pdb(pdb_path: Path) -> tuple[str, dict]:
    result = process_single_pdb(pdb_path, 0.0)
    return pdb_path.name, result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Calculate diffusion lengths for every PDB in a folder."
    )
    parser.add_argument(
        "--pdb-dir",
        type=Path,
        default=PDB_DIRECTORY,
        help=f"Folder containing PDB files (default: {PDB_DIRECTORY})",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=DEFAULT_WORKERS,
        help=f"Parallel worker processes (default: {DEFAULT_WORKERS})",
    )
    args = parser.parse_args()

    if args.workers < 1:
        parser.error("--workers must be at least 1")
    if not args.pdb_dir.is_dir():
        parser.error(f"PDB folder does not exist: {args.pdb_dir}")

    pdb_files = sorted(args.pdb_dir.glob("*.pdb"))
    if not pdb_files:
        parser.error(f"No .pdb files found in: {args.pdb_dir}")

    results: dict[str, dict] = {}
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(_process_one_pdb, path) for path in pdb_files]
        for completed, future in enumerate(as_completed(futures), start=1):
            pdb_name, result = future.result()
            results[pdb_name] = result
            if completed % 25 == 0 or completed == len(pdb_files):
                print(f"Processed {completed}/{len(pdb_files)} PDB files.", flush=True)

    output_path = args.pdb_dir / "diffusion_lengths.csv"
    with output_path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.writer(output_file)
        writer.writerow(
            (
                "pdb_file",
                "diffusion_length_angstrom",
                "total_pairs",
                "coherent_pairs",
                "coherent_pair_fraction",
            )
        )
        writer.writerows(
            (
                path.name,
                results[path.name]["diffusion_length_angstrom"],
                results[path.name]["total_pairs"],
                results[path.name]["coherent_pairs"],
                results[path.name]["coherent_pair_fraction"],
            )
            for path in pdb_files
        )

    print(f"Saved {len(results)} results to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
