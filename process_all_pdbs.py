"""Run the FRET and DET PDB calculations in parallel and save both results."""

from __future__ import annotations

import argparse
import csv
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# Keep these two modules side-by-side with this script:
#   process_single_pdb.py  -> FRET calculation
#   det_code_clean.py      -> DET calculation
from process_single_pdb import process_single_pdb as process_fret_pdb
from det_code_clean import process_single_pdb as process_det_pdb


PDB_DIRECTORY = Path(r"D:\DJREDDY\EGFR Isha\EGFR_mutants_20260903")
DEFAULT_WORKERS = max(1, min(4, os.cpu_count() or 1))


def _process_one_pdb(pdb_path: Path) -> tuple[str, dict, dict]:
    """Run both models on one PDB file."""
    fret_result = process_fret_pdb(pdb_path)
    det_result = process_det_pdb(pdb_path)
    return pdb_path.name, fret_result, det_result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Calculate FRET and DET diffusion lengths for every PDB in a folder."
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

    results: dict[str, tuple[dict, dict]] = {}
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(_process_one_pdb, path) for path in pdb_files]
        for completed, future in enumerate(as_completed(futures), start=1):
            pdb_name, fret_result, det_result = future.result()
            results[pdb_name] = (fret_result, det_result)
            if completed % 25 == 0 or completed == len(pdb_files):
                print(f"Processed {completed}/{len(pdb_files)} PDB files.", flush=True)

    output_path = args.pdb_dir / "fret_det_diffusion_lengths.csv"
    columns = (
        "pdb_file",
        "fret_diffusion_length_angstrom",
        "fret_total_pairs",
        "fret_coherent_pairs",
        "fret_coherent_pair_fraction",
        "det_diffusion_length_angstrom",
        "det_total_pairs",
        "det_coherent_pairs",
        "det_coherent_pair_fraction",
    )

    with output_path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.writer(output_file)
        writer.writerow(columns)

        for path in pdb_files:
            fret_result, det_result = results[path.name]
            writer.writerow(
                (
                    path.name,
                    fret_result.get("diffusion_length_angstrom"),
                    fret_result.get("total_pairs"),
                    fret_result.get("coherent_pairs"),
                    fret_result.get("coherent_pair_fraction"),
                    det_result.get("diffusion_length_angstrom"),
                    det_result.get("total_pairs"),
                    det_result.get("coherent_pairs"),
                    det_result.get("coherent_pair_fraction"),
                )
            )

    print(f"Saved {len(results)} FRET + DET results to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
