"""Python translation of det_code_clean.m."""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist


RESIDUE_TYPES = ("TRP", "TYR", "PHE")
TRANSFER_PARAMETERS = {
    ("TRP", "TRP"): 0.00058,
    ("TYR", "TYR"): 0.0012,
    ("PHE", "PHE"): 0.00098,
    ("TRP", "TYR"): 0.000487,
    ("TYR", "TRP"): 0.0011,
    ("PHE", "TYR"): 0.0012,
    ("PHE", "TRP"): 0.0039,
    ("TRP", "PHE"): 0.00031,
    ("TYR", "PHE"): 0.00098,
}
LIFETIMES_MICROSECONDS = {"TRP": 7.3, "TYR": 5.6, "PHE": 3.1}
EXCITATION_COUNT = 10_000
MAX_TIME_STEPS = 1_000


def _read_pdb_atoms(
    pdb_filename: str | Path,
) -> list[tuple[str, int, str, float, float, float]]:
    atoms = []
    saw_model = False

    with Path(pdb_filename).open("r", encoding="ascii", errors="replace") as pdb_file:
        for line in pdb_file:
            record = line[0:6].strip().upper()
            if record == "MODEL":
                if saw_model:
                    break
                saw_model = True
                continue
            if record == "ENDMDL" and saw_model:
                break
            if record == "END":
                break
            if record not in {"ATOM", "HETATM"}:
                continue

            try:
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip()
                residue_number = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
            except (ValueError, IndexError) as error:
                raise ValueError(f"Invalid PDB atom record: {line.rstrip()}") from error

            atoms.append((residue_name, residue_number, atom_name, x, y, z))

    return atoms


def _group_residues(
    atoms: list[tuple[str, int, str, float, float, float]],
) -> tuple[list[str], list[tuple[list[str], np.ndarray]]]:
    residue_types: list[str] = []
    groups: list[tuple[list[str], np.ndarray]] = []

    for residue_type in RESIDUE_TYPES:
        matching_atoms = [
            atom for atom in atoms if residue_type in atom[0].upper()
        ]
        residue_numbers = sorted({atom[1] for atom in matching_atoms})

        for residue_number in residue_numbers:
            residue_atoms = [
                atom for atom in matching_atoms if atom[1] == residue_number
            ]
            atom_names = [atom[2] for atom in residue_atoms]
            coordinates = np.asarray(
                [[atom[3], atom[4], atom[5]] for atom in residue_atoms],
                dtype=float,
            ).reshape((-1, 3))
            residue_types.append(residue_type)
            groups.append((atom_names, coordinates))

    return residue_types, groups


def _mean_omit_nan(values: np.ndarray) -> np.ndarray:
    means = np.empty(values.shape[1], dtype=float)
    for column in range(values.shape[1]):
        valid_values = values[:, column]
        valid_values = valid_values[~np.isnan(valid_values)]
        means[column] = np.mean(valid_values) if len(valid_values) else np.nan
    return means


def _build_transfer_parameter_matrix(residue_types: list[str]) -> np.ndarray:
    count = len(residue_types)
    matrix = np.zeros((count, count), dtype=float)
    for row, donor_type in enumerate(residue_types):
        for column, acceptor_type in enumerate(residue_types):
            matrix[row, column] = TRANSFER_PARAMETERS[
                (donor_type, acceptor_type)
            ]
    return matrix


def _calculate_transfer_rates(
    distances_matrix_m: np.ndarray,
    transfer_parameter_matrix: np.ndarray,
) -> np.ndarray:
    det_coupling_extrapolated_ev = (3.22e-2 / (1.4 * 1.4)) * np.exp(
        -(1e10 * 2.5579) * (distances_matrix_m - 4.15e-10)
    )
    det_coupling_extrapolated_cm_minus_1 = (
        det_coupling_extrapolated_ev * (1.602e-19) * (5.03 * 10**22)
    )
    det_coupling_extrapolated_cm_minus_1 = (
        det_coupling_extrapolated_cm_minus_1
        - np.diag(np.diag(det_coupling_extrapolated_cm_minus_1))
    )

    hbar_cm_minus_1 = 5.29e-12
    speed_of_light_cm_per_second = 3e10
    det_rates_per_second = (
        1.0 / (speed_of_light_cm_per_second * hbar_cm_minus_1**2)
    ) * (
        transfer_parameter_matrix * det_coupling_extrapolated_cm_minus_1**2
    )

    cutoff_distance_m = math.inf
    distances_matrix_cutoff_m = distances_matrix_m.copy()
    distances_matrix_cutoff_m[
        distances_matrix_cutoff_m > cutoff_distance_m
    ] = np.nan
    det_rates_per_second_cutoff = det_rates_per_second.copy()
    det_rates_per_second_cutoff[np.isnan(distances_matrix_cutoff_m)] = np.nan

    det_rates_per_microsecond = (10**-6) * det_rates_per_second_cutoff
    det_rates_per_microsecond[
        np.isinf(det_rates_per_microsecond)
        | np.isnan(det_rates_per_microsecond)
    ] = 0.0
    return det_rates_per_microsecond


def _simulate_diffusion_length(
    distances_matrix_angstrom: np.ndarray,
    det_rates_per_microsecond: np.ndarray,
    residue_types: list[str],
    *,
    random_seed: int | None = None,
) -> float:
    lifetime_chromophores_microseconds = np.asarray(
        [LIFETIMES_MICROSECONDS[residue_type] for residue_type in residue_types],
        dtype=float,
    )
    radiative_rates_per_microsecond = 1.0 / lifetime_chromophores_microseconds
    rates_per_microsecond = np.hstack(
        (
            det_rates_per_microsecond,
            radiative_rates_per_microsecond.reshape((-1, 1)),
        )
    )

    rng = np.random.default_rng(random_seed)
    distances_travelled_as_function_of_time = np.zeros(
        (EXCITATION_COUNT, MAX_TIME_STEPS + 1), dtype=float
    )
    where_it_moved_all_values = np.zeros(
        (EXCITATION_COUNT, MAX_TIME_STEPS + 1), dtype=int
    )
    starting_molecules = np.zeros(EXCITATION_COUNT, dtype=int)

    for excitation_index in range(EXCITATION_COUNT):
        starting_molecule = int(rng.integers(len(residue_types)))
        starting_molecules[excitation_index] = starting_molecule
        where_it_moved_all_values[excitation_index, 0] = starting_molecule + 1
        where_it_moved_location = starting_molecule

        for time_index in range(MAX_TIME_STEPS):
            weights = rates_per_microsecond[where_it_moved_location]
            total_weight = float(np.sum(weights))
            if total_weight <= 0 or not math.isfinite(total_weight):
                raise ValueError(
                    "Transition weights must have a positive finite sum."
                )
            threshold = rng.random() * total_weight
            acceptor_index = int(
                np.searchsorted(np.cumsum(weights), threshold, side="left")
            )

            if acceptor_index == len(weights) - 1:
                break

            where_it_moved_location = acceptor_index
            where_it_moved_all_values[
                excitation_index, time_index + 1
            ] = acceptor_index + 1
            distances_travelled_as_function_of_time[
                excitation_index, time_index + 1
            ] = distances_matrix_angstrom[
                starting_molecule, where_it_moved_location
            ]

    where_it_moved_all_values_last_value = np.max(
        np.where(
            where_it_moved_all_values != 0,
            np.arange(1, where_it_moved_all_values.shape[1] + 1),
            0,
        ),
        axis=1,
    )
    ending_molecules = np.empty(EXCITATION_COUNT, dtype=int)
    for excitation_index, last_column in enumerate(
        where_it_moved_all_values_last_value
    ):
        ending_molecules[excitation_index] = (
            where_it_moved_all_values[excitation_index, last_column - 1] - 1
        )

    distances_travelled_by_each_excitation = np.zeros(EXCITATION_COUNT, dtype=float)
    for excitation_index in range(EXCITATION_COUNT):
        distances_travelled_by_each_excitation[excitation_index] = (
            distances_matrix_angstrom[
                starting_molecules[excitation_index],
                ending_molecules[excitation_index],
            ]
        )

    return float(np.mean(distances_travelled_by_each_excitation))


def process_single_pdb(
    pdb_filename: str | Path,
    microtubule_diameter_angstrom: float,
    *,
    random_seed: int | None = None,
) -> float:
    """Return the DET-model mean endpoint displacement in Angstroms."""
    try:
        microtubule_diameter = float(microtubule_diameter_angstrom)
        microtubule_radius = microtubule_diameter / 2.0
        _ = microtubule_radius

        atoms = _read_pdb_atoms(pdb_filename)
        residue_types, groups = _group_residues(atoms)
        number_of_molecules = len(residue_types)
        if number_of_molecules < 2:
            return 0.0

        positions_xyz = np.asarray(
            [_mean_omit_nan(coordinates) for _, coordinates in groups],
            dtype=float,
        )
        distances_matrix_angstrom = cdist(
            positions_xyz, positions_xyz, metric="euclidean"
        )
        distances_matrix_m = distances_matrix_angstrom * 10**-10

        transfer_parameter_matrix = _build_transfer_parameter_matrix(residue_types)
        det_rates_per_microsecond = _calculate_transfer_rates(
            distances_matrix_m,
            transfer_parameter_matrix,
        )
        total_hops = int(np.sum(~np.isnan(det_rates_per_microsecond)))
        _ = total_hops

        return _simulate_diffusion_length(
            distances_matrix_angstrom,
            det_rates_per_microsecond,
            residue_types,
            random_seed=random_seed,
        )
    except Exception as error:
        print(f"Error processing {pdb_filename}: {error}", file=sys.stderr)
        return float("nan")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Calculate the DET diffusion length for one PDB structure."
    )
    parser.add_argument("pdb_filename", help="Input PDB file")
    parser.add_argument(
        "microtubule_diameter_angstrom",
        type=float,
        help="Microtubule diameter in Angstroms (unused by the MATLAB algorithm)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Optional Python RNG seed for repeatable Monte Carlo runs",
    )
    args = parser.parse_args()
    print(
        process_single_pdb(
            args.pdb_filename,
            args.microtubule_diameter_angstrom,
            random_seed=args.seed,
        )
    )


if __name__ == "__main__":
    main()
