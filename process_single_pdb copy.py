"""Python translation of ``process_single_pdb_v1.m``.

Requires NumPy and SciPy (install with ``pip install numpy scipy``). Monte Carlo
draws use NumPy's generator; transition probabilities and stopping rules match
the MATLAB implementation, but its random stream cannot be reproduced without
MATLAB's RNG state.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import svd
from scipy.spatial.distance import cdist


RING_ATOMS = {
    "TRP": {"CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
}


def _read_pdb_atoms(pdb_filename: str | Path) -> list[tuple[str, int, str, float, float, float]]:
    """Read the first PDB model in the fixed-column format used by pdbread."""
    atoms = []
    saw_model = False
    inside_first_model = True

    with Path(pdb_filename).open("r", encoding="ascii", errors="replace") as pdb_file:
        for line in pdb_file:
            record = line[0:6].strip().upper()
            if record == "MODEL":
                if saw_model:
                    break
                saw_model = True
                inside_first_model = True
                continue
            if record == "ENDMDL" and saw_model:
                break
            if record == "END":
                break
            if not inside_first_model or record not in {"ATOM", "HETATM"}:
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


def _matlab_residue_groups(atoms, residue_name):
    matching = [atom for atom in atoms if residue_name in atom[0].upper()]
    residue_numbers = sorted({atom[1] for atom in matching})
    groups = []
    for residue_number in residue_numbers:
        group = [atom for atom in matching if atom[1] == residue_number]
        groups.append(
            (
                residue_number,
                [atom[2] for atom in group],
                np.asarray([[atom[3], atom[4], atom[5]] for atom in group], dtype=float)
                .reshape((-1, 3)),
            )
        )
    return groups


def _build_parameter_matrices(types):
    count = len(types)
    j_values = {
        ("TRP", "TRP"): 5.5953e-06,
        ("TYR", "TYR"): 4.7933e-04,
        ("PHE", "PHE"): 6.1123e-04,
        ("TRP", "TYR"): 2.3247e-05,
        ("TYR", "TRP"): 9.2550e-04,
        ("PHE", "TYR"): 0.0043,
        ("PHE", "TRP"): 0.0039,
        ("TRP", "PHE"): 1.3034e-07,
        ("TYR", "PHE"): 1.8891e-05,
    }
    dipoles = {"TRP": 2.074, "TYR": 1.18, "PHE": 0.28}
    debye_units = 3.33564e-30

    j_matrix = np.empty((count, count), dtype=float)
    mu_squared_matrix = np.empty((count, count), dtype=float)
    for row, donor_type in enumerate(types):
        for column, acceptor_type in enumerate(types):
            j_matrix[row, column] = j_values[(donor_type, acceptor_type)]
            mu_squared_matrix[row, column] = (
                dipoles[donor_type] * dipoles[acceptor_type] * debye_units**2
            )
    return j_matrix, mu_squared_matrix


def _calculate_ring_normals(groups, types, positions_xyz):
    ring_normals = np.zeros((len(groups), 3), dtype=float)
    for index, (residue_type, group) in enumerate(zip(types, groups)):
        _, atom_names, coordinates = group
        ring_coordinates = np.asarray(
            [
                coordinate
                for atom_name, coordinate in zip(atom_names, coordinates)
                if atom_name in RING_ATOMS[residue_type]
            ],
            dtype=float,
        ).reshape((-1, 3))
        ring_coordinates = ring_coordinates[~np.isnan(ring_coordinates).any(axis=1)]

        if len(ring_coordinates) == 0:
            raise ValueError("Cannot calculate a ring normal without ring atoms.")

        ring_centroid = np.mean(ring_coordinates, axis=0)
        _, _, right_vectors = svd(ring_coordinates - ring_centroid, full_matrices=False)
        normal = right_vectors[-1]
        normal = normal / np.linalg.norm(normal)

        toward_residue_centroid = positions_xyz[index] - ring_centroid
        if np.dot(normal, toward_residue_centroid) < 0:
            normal = -normal
        ring_normals[index] = normal
    return ring_normals


def process_single_pdb(
    pdb_filename,
    microtubule_diameter_angstrom,
    *,
    random_seed: int | None = None,
):
    """Return diffusion length in Angstroms; ``random_seed`` enables repeatable runs."""
    try:
        # These values are intentionally computed but unused in the MATLAB code.
        microtubule_diameter = float(microtubule_diameter_angstrom)
        microtubule_radius = microtubule_diameter / 2.0
        del microtubule_radius

        atoms = _read_pdb_atoms(pdb_filename)
        if not atoms:
            return 0.0

        residue_names = ("TRP", "TYR", "PHE")
        groups_by_type = {
            residue_type: _matlab_residue_groups(atoms, residue_type)
            for residue_type in residue_names
        }
        type_groups = [
            (residue_type, group)
            for residue_type in residue_names
            for group in groups_by_type[residue_type]
        ]
        number_of_molecules = len(type_groups)
        if number_of_molecules < 2:
            return 0.0

        centroids = []
        for _, (_, _, coordinates) in type_groups:
            if len(coordinates):
                centroids.append(np.mean(coordinates, axis=0))
            else:
                centroids.append(np.full(3, np.nan))
        positions_xyz = np.asarray(centroids, dtype=float)

        types = [residue_type for residue_type, _ in type_groups]
        group_data = [group for _, group in type_groups]
        ring_normals = _calculate_ring_normals(group_data, types, positions_xyz)

        distances_matrix = cdist(positions_xyz, positions_xyz, metric="euclidean")
        kappa_matrix = np.zeros((number_of_molecules, number_of_molecules), dtype=float)
        for i in range(number_of_molecules):
            for j in range(number_of_molecules):
                if i == j:
                    continue
                r_hat = (positions_xyz[j] - positions_xyz[i]) / distances_matrix[i, j]
                kappa_matrix[i, j] = (
                    np.dot(ring_normals[i], ring_normals[j])
                    - 3.0
                    * np.dot(ring_normals[i], r_hat)
                    * np.dot(ring_normals[j], r_hat)
                )

        j_matrix, mu_squared_matrix = _build_parameter_matrices(types)
        epsilon_0 = 8.854e-12
        constants = 1.0 / (4.0 * math.pi * epsilon_0)
        distances_matrix_m = distances_matrix * 1e-10
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            coupling_joules = (
                constants
                * kappa_matrix
                * mu_squared_matrix
                / distances_matrix_m**3
            )
            coupling_cm_minus_1 = coupling_joules * (5.03e22)
            fermi_rates_per_second = (
                1.0
                / (3.0e10 * (5.29e-12) ** 2)
                * (1.0 / 1.4**4)
                * (j_matrix * coupling_cm_minus_1**2)
            )
        fermi_rates_per_ns = fermi_rates_per_second * 1e-9
        fermi_rates_per_ns[~np.isfinite(fermi_rates_per_ns)] = 0.0
        radiative_rates_per_ns = np.asarray(
            [1.0 / {"TRP": 3.6, "TYR": 7.5, "PHE": 3.39}[kind] for kind in types]
        )
        rates = np.column_stack((fermi_rates_per_ns, radiative_rates_per_ns))

        if not np.isfinite(rates).all() or np.any(rates < 0):
            raise ValueError("The transition-rate matrix contains invalid weights.")
        rng = np.random.default_rng(random_seed)
        distances_travelled = np.zeros(10_000, dtype=float)

        # Keep MATLAB's nested loop order: finish one excitation before
        # drawing the next starting molecule.
        for excitation in range(10_000):
            starting_molecule = int(rng.integers(number_of_molecules))
            current_molecule = starting_molecule

            for _ in range(1_000):
                weights = rates[current_molecule]
                total_weight = float(weights.sum())
                if total_weight <= 0.0 or not math.isfinite(total_weight):
                    raise ValueError(
                        "Transition weights must have a positive finite sum."
                    )

                threshold = rng.random() * total_weight
                acceptor_index = int(
                    np.searchsorted(np.cumsum(weights), threshold, side="left")
                )
                if acceptor_index == number_of_molecules:
                    break
                current_molecule = acceptor_index

            distances_travelled[excitation] = distances_matrix[
                starting_molecule, current_molecule
            ]

        return float(np.mean(distances_travelled))
    except Exception as error:
        print(f"Error processing {pdb_filename}: {error}", file=sys.stderr)
        return float("nan")


def main():
    parser = argparse.ArgumentParser(
        description="Calculate a PDB diffusion length using the MATLAB-translated pipeline."
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
