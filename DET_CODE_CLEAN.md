# `det_code_clean.m` Python translation

`det_code_clean.py` is a separate, step-by-step Python translation of
`det_code_clean.m`. It intentionally does **not** replace
`process_single_pdb.py`: that module contains the project's other
single-structure calculation, with different coupling constants, equations,
lifetimes, and time units.

## Run it

From the `Plexcal` directory:

```powershell
python -m pip install -r requirements.txt
python det_code_clean.py "1TUB.pdb" 0
```

The second positional argument is the microtubule diameter in Angstroms. The
MATLAB code calculates the corresponding radius but does not use either value
in its distance or transfer calculations. An optional seed makes a Python
run repeatable:

```powershell
python det_code_clean.py "1TUB.pdb" 0 --seed 123
```

The Python and MATLAB random-number generators differ, so this does not
reproduce MATLAB's random stream or bit-for-bit output.

## Workflow and equations

The Python function `process_single_pdb(pdb_filename,
microtubule_diameter_angstrom, random_seed=None)` follows the MATLAB stages:

1. **Read PDB records.** Read `ATOM` and `HETATM` atom names, residue names,
   residue sequence numbers, and Cartesian coordinates from fixed-width PDB
   fields. The first `MODEL` is used when models are present. Coordinates are
   in Angstroms.
2. **Select and group chromophores.** Select residue names containing `TRP`,
   `TYR`, or `PHE` (case-insensitive); group each type by residue sequence
   number and process types in TRP, TYR, PHE order. Chain ID and insertion
   code are not included in the grouping key, matching the MATLAB code's
   `resSeq` grouping.
3. **Count and handle small inputs.** If fewer than two selected residues
   exist, the returned dictionary reports a diffusion length of `0.0`, zero
   total and coherent pairs, and a `None` coherent-pair fraction.
4. **Calculate residue positions.** For each residue, take the coordinate-wise
   mean of all atom coordinates, omitting NaNs as MATLAB's `mean(...,
   'omitnan')` does:

   \[
   \mathbf r_i = \operatorname{mean}_{a \in i}(\mathbf x_{ia}).
   \]

5. **Calculate pair distances.** `pdist2` is translated with Euclidean
   distances between every pair of residue centroids:

   \[
   r_{ij,\mathrm{\AA}} = \lVert\mathbf r_i-\mathbf r_j\rVert_2,
   \qquad
   r_{ij,\mathrm m}=10^{-10}r_{ij,\mathrm{\AA}}.
   \]

   The positions and output distances stay in Angstroms; the coupling model
   uses the converted metre distances.
6. **Build the directional \(J\) matrix.** Rows are donors and columns are
   acceptors. These constants are transcribed from the MATLAB source:

   | Donor \ Acceptor | TRP | TYR | PHE |
   |---|---:|---:|---:|
   | TRP | 0.00058 | 0.000487 | 0.00031 |
   | TYR | 0.0011 | 0.0012 | 0.00098 |
   | PHE | 0.0039 | 0.0012 | 0.00098 |

7. **Calculate distance-dependent DET coupling.** For each distance, the
   MATLAB expression is

   \[
   V_{ij,\mathrm{eV}} =
   \frac{3.22\times10^{-2}}{1.4\times1.4}
   \exp\!\left[-(10^{10}\times2.5579)
   (r_{ij,\mathrm m}-4.15\times10^{-10})\right].
   \]

   It is converted using the exact source-code factors:

   \[
   V_{ij,\mathrm{cm}^{-1}} =
   V_{ij,\mathrm{eV}}(1.602\times10^{-19})(5.03\times10^{22}).
   \]

   The diagonal of the converted coupling matrix is then set to zero by
   subtracting its diagonal matrix.
8. **Calculate transfer rates.** With
   \(c=3\times10^{10}\) and
   \(\hbar_{\mathrm{cm}^{-1}}=5.29\times10^{-12}\), the code uses

   \[
   k_{\mathrm{DET},ij,\mathrm{s}^{-1}} =
   \frac{1}{c\,\hbar_{\mathrm{cm}^{-1}}^2}
   J_{ij}V_{ij,\mathrm{cm}^{-1}}^2.
   \]

   The distance cutoff is `Inf`, so no finite pair distance is excluded. The
   rates are converted to per-microsecond values by multiplying by
   \(10^{-6}\). NaN and infinite transfer rates are then replaced with zero.
9. **Add radiative decay rates.** The source uses lifetimes of 7.3
   microseconds (TRP), 5.6 microseconds (TYR), and 3.1 microseconds (PHE):

   \[
   k_{\mathrm{rad},i,\mathrm{\mu s}^{-1}} = \frac{1}{\tau_i}.
   \]

   Each row of the event-weight matrix contains that chromophore's transfer
   rates to all acceptors followed by its radiative rate.
10. **Run the Monte Carlo walk.** Run 10,000 excitations. Each starts at a
    uniformly sampled chromophore. For each of up to 1,000 time steps, sample
    one event using the row's nonnegative rate weights. Equivalently, the
    probability of an event with rate \(k_q\) is

    \[
    P(q\mid i)=\frac{k_q}{\sum_p k_p}.
    \]

    A transfer event changes the current chromophore. Selection of the final
    radiative-rate column ends the excitation. The implementation retains the
    MATLAB trajectory matrices and records the displacement from the original
    start to each visited chromophore.
11. **Calculate the reported result.** Find the last nonzero chromophore index
    for each trajectory, calculate the straight-line distance from its start
    to that endpoint, and average over all excitations:

    \[
    L_{\mathrm{\AA}} =
    \frac{1}{10^4}\sum_{m=1}^{10^4}
    \left\lVert\mathbf r_{\mathrm{end},m}
    -\mathbf r_{\mathrm{start},m}\right\rVert_2.
    \]

    Thus, the reported value is mean start-to-end displacement, not total hop
    path length or FRET efficiency.

The function returns a dictionary with `diffusion_length_angstrom`,
`total_pairs`, `coherent_pairs`, and `coherent_pair_fraction`. Coherent DET
pairs use the strict threshold `abs(coupling) > 53.052` cm⁻¹ and only the
strict lower triangle of the coupling matrix, counting each off-diagonal pair
once. The fraction is `coherent_pairs / total_pairs`; it is `None` when the
total pair count is zero. The command-line interface continues to print only
the diffusion length.

## Error handling and implementation details

- The only MATLAB function in the selected file,
  `process_single_pdb`, is exposed under the same name in
  `det_code_clean.py`. Python helper functions split out parsing, grouping,
  \(J\)-matrix construction, rate calculation, and simulation so each stage
  can be inspected independently.
- As in the MATLAB `try/catch`, any processing exception is printed to
  standard error and the `diffusion_length_angstrom` result is `NaN`. Fewer
  than two selected residues reports a diffusion length of `0.0`.
- The source computes `Total_hops`, the microtubule radius, and the
  per-time-step travelled-distance matrix, though they do not affect the
  returned average. The Python translation retains these calculations and
  intermediate trajectory data.
- Random-number sequences differ between MATLAB and Python. The algorithm's
  rate weights and stopping rules are carried over, but a seeded Python result
  is not a bit-for-bit MATLAB reproduction.
- The constants and unit conversions above are reproduced as written in the
  MATLAB source; this guide documents the implementation and does not
  independently validate its physical assumptions or units.

## Difference from the other Python processor

Use `det_code_clean.py` for the distance-exponential DET model documented
above. Use `process_single_pdb.py` for the existing calculation documented in
[PROCESS_SINGLE_PDB.md](PROCESS_SINGLE_PDB.md). They are separate calculations
and their outputs should not be treated as interchangeable.
