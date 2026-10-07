# Single-PDB Diffusion-Length Calculation

This guide explains what `process_single_pdb.py` reads, how it calculates
energy-transfer and radiative rates, and how it turns those rates into the
reported diffusion length. The implementation is a Python translation of
`process_single_pdb_v1.m`. It describes the model implemented in those files;
it is not an independent validation of the physical assumptions or parameter
values.

## What the program calculates

For one PDB structure, the program:

1. Reads atom coordinates for tryptophan (TRP), tyrosine (TYR), and
   phenylalanine (PHE) residues.
2. Treats each selected residue as a chromophore.
3. Calculates each chromophore's centroid, ring-plane normal, and distance to
   every other chromophore.
4. Calculates a geometry-dependent transfer rate for each ordered
   donor-to-acceptor pair, and a radiative-decay rate for each chromophore.
5. Simulates 10,000 excitation histories. Each history starts at a randomly
   selected chromophore and makes transfer hops until radiative decay or the
   1,000-step limit.
6. Averages the start-to-end straight-line distance over those histories.

The returned value is a **mean displacement in Angstroms**, not a path length
(the sum of hop distances), a diffusion coefficient, or a conventional FRET
efficiency.

## Requirements and running one structure

From the `Plexcal` directory, install the calculation dependencies if needed:

```powershell
python -m pip install -r requirements.txt
```

The requirements are NumPy and SciPy. Run a structure by passing its path and
the required microtubule diameter argument:

```powershell
python process_single_pdb.py "1TUB.pdb" 0
```

An optional seed makes that individual Monte Carlo run repeatable:

```powershell
python process_single_pdb.py "1TUB.pdb" 0 --seed 123
```

The diameter is converted to a number and its radius is computed, but neither
value is used by the calculation. Changing the diameter therefore does not
change the result. The Python and MATLAB implementations use different random
number generators, so the same seed does not make their results match.

## Input interpretation

The parser uses the fixed-width PDB fields for atom name, residue name, residue
sequence number, and x/y/z coordinates. It accepts `ATOM` and `HETATM` records
and uses the first model when a file contains `MODEL` records. PDB coordinates
are treated as Angstroms.

Only residues whose names contain `TRP`, `TYR`, or `PHE` are included. Their
atoms are grouped by residue name and residue sequence number; chain ID and
insertion code are not part of the grouping key. The code does not filter
alternate atom locations. These details can affect structures with reused
residue numbering across chains or alternate conformations.

If the file has no atom records, or has fewer than two selected residues, the
function returns `0.0`. If parsing or calculation raises an exception, it
prints an error to standard error and returns `NaN`.

## Geometry and physical model

### 1. Chromophore positions and pair distances

For chromophore \(i\), the position is the arithmetic mean of the coordinates
of all atoms in that residue:

\[
\mathbf{r}_i = \frac{1}{N_i}\sum_{a=1}^{N_i}\mathbf{x}_{ia}
\]

The pair distance is the Euclidean distance between residue centroids:

\[
r_{ij} = \left\|\mathbf{r}_j-\mathbf{r}_i\right\|
\]

These distances are initially in Angstroms. Before the coupling calculation,
the code converts them to metres using
\(r_{ij,\mathrm m}=10^{-10}r_{ij,\mathrm{\AA}}\). Distances used in the final
displacement remain in Angstroms.

### 2. Aromatic-ring orientation

For each chromophore, only the following atom names are used to fit its
aromatic ring plane:

| Residue | Ring atom names |
|---|---|
| TRP | CG, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3, CH2 |
| TYR | CG, CD1, CD2, CE1, CE2, CZ |
| PHE | CG, CD1, CD2, CE1, CE2, CZ |

The ring-atom coordinates are centered by subtracting their mean. Singular
value decomposition (SVD) finds the plane of best fit; the right-singular
vector associated with the smallest singular value is the ring-plane normal
\(\mathbf{n}_i\). The normal is normalized and its sign is chosen to point from
the ring centroid toward the centroid of the whole residue.

At least one recognized ring atom is required. A missing ring atom set that
prevents a normal from being calculated causes the calculation to fail for
that PDB.

### 3. Orientation factor

For each ordered pair \(i\ne j\), the unit vector from donor \(i\) to acceptor
\(j\) is

\[
\hat{\mathbf{r}}_{ij} =
\frac{\mathbf{r}_j-\mathbf{r}_i}{r_{ij}}.
\]

The code calculates the dipole-orientation factor

\[
\kappa_{ij} =
\mathbf{n}_i\cdot\mathbf{n}_j
-3(\mathbf{n}_i\cdot\hat{\mathbf{r}}_{ij})
(\mathbf{n}_j\cdot\hat{\mathbf{r}}_{ij}).
\]

The diagonal is set to zero, so a chromophore does not transfer to itself.

### 4. Dipole products and pair parameter \(J_{ij}\)

The dipole magnitudes in the code are 2.074 D for TRP, 1.18 D for TYR, and
0.28 D for PHE. A dipole value in Debye is converted to SI units with
\(1\ \mathrm{D}=3.33564\times10^{-30}\ \mathrm{C\,m}\). For each pair, the
code uses the product of the donor and acceptor dipole magnitudes:

\[
\mu_{ij} = \mu_i\mu_j.
\]

The ordered pair parameter matrix is directional; rows are donor residue type
and columns are acceptor residue type:

| Donor \ Acceptor | TRP | TYR | PHE |
|---|---:|---:|---:|
| TRP | \(5.5953\times10^{-6}\) | \(2.3247\times10^{-5}\) | \(1.3034\times10^{-7}\) |
| TYR | \(9.2550\times10^{-4}\) | \(4.7933\times10^{-4}\) | \(1.8891\times10^{-5}\) |
| PHE | \(3.9\times10^{-3}\) | \(4.3\times10^{-3}\) | \(6.1123\times10^{-4}\) |

The values are copied from the MATLAB translation. The calculation does not
derive them from the PDB coordinates.

### 5. Dipole coupling and transfer rate

The code first calculates a pairwise coupling in joules:

\[
V_{ij,\mathrm J} =
\frac{1}{4\pi\epsilon_0}
\frac{\kappa_{ij}\mu_{ij}}{r_{ij,\mathrm m}^{3}},
\qquad \epsilon_0=8.854\times10^{-12}.
\]

It converts this coupling to inverse centimetres by multiplying by
\(5.03\times10^{22}\):

\[
V_{ij,\mathrm{cm}^{-1}} =
V_{ij,\mathrm J}(5.03\times10^{22}).
\]

The transfer rate in inverse seconds is then calculated using the implemented
Fermi-rate expression:

\[
k_{\mathrm{ET},ij,\mathrm{s}^{-1}} =
\frac{1}{c\,\hbar_{\mathrm{cm}^{-1}}^2}
\frac{1}{n^4}
J_{ij}V_{ij,\mathrm{cm}^{-1}}^2,
\]

with \(c=3\times10^{10}\), \(\hbar_{\mathrm{cm}^{-1}}=5.29\times10^{-12}\),
and protein refractive index \(n=1.4\). The code converts the result to
inverse nanoseconds by multiplying by \(10^{-9}\). Non-finite transfer-rate
entries are set to zero.

### 6. Radiative decay and competing event probabilities

The radiative lifetimes used are 3.6 ns (TRP), 7.5 ns (TYR), and 3.39 ns
(PHE). For chromophore \(i\), the radiative rate is

\[
k_{\mathrm{rad},i} = \frac{1}{\tau_i}.
\]

The transfer rates from chromophore \(i\) to each possible acceptor and its
radiative rate form the event weights. The probability of choosing event \(q\)
is its rate divided by the sum of all event rates from the current
chromophore:

\[
P(i\to q) = \frac{k_{i\to q}}{\sum_j k_{\mathrm{ET},ij}
+ k_{\mathrm{rad},i}}.
\]

Here \(q\) can be an acceptor transfer or radiative decay. The program samples
from the equivalent unnormalized rates (normalizing implicitly by drawing a
uniform threshold between zero and their sum). A radiative-decay selection
ends that excitation; a transfer selection moves it to the selected acceptor.

### 7. Monte Carlo estimate of the reported length

The Python implementation performs 10,000 independent excitation histories.
Each begins at a uniformly chosen chromophore. At each of up to 1,000
iterations, it samples a transfer or radiative event using the current
chromophore's rates. The simulation ends on radiative decay or when the
iteration limit is reached.

For history \(m\), let \(s_m\) be its starting chromophore and \(e_m\) its
ending chromophore. Its recorded displacement is

\[
d_m = \left\|\mathbf{r}_{e_m}-\mathbf{r}_{s_m}\right\|.
\]

The function returns the arithmetic mean

\[
L = \frac{1}{10{,}000}\sum_{m=1}^{10{,}000} d_m
\quad\text{(Angstroms)}.
\]

This is the endpoint displacement, not the sum of all distances traveled
during the hops. Because the generator is unseeded by default, repeated runs
can produce slightly different Monte Carlo estimates.

## Processing a folder in parallel

To run the same function over every `.pdb` file in `all_proteins_pdb/`, use
the separate batch script:

```powershell
python process_all_pdbs.py
```

The batch script calls the existing `process_single_pdb()` function without
changing its calculation, uses up to four worker processes by default, and
writes results in sorted filename order to
`all_proteins_pdb/diffusion_lengths.csv`. The CSV columns are
`pdb_file` and `diffusion_length_angstrom`. A calculation error is represented
by `NaN` in that file, in line with the single-PDB function's return value.
The output file is overwritten when the batch is run again.

Choose another worker count with, for example:

```powershell
python process_all_pdbs.py --workers 2
```

To process a different folder, pass its path:

```powershell
python process_all_pdbs.py --pdb-dir "D:\structures\pdbs" --workers 2
```

## Implementation caveats

- The diameter and radius arguments are currently unused, so this model does
  not impose a microtubule cylinder or geometric cutoff.
- Residue grouping uses residue type and sequence number only, not chain ID or
  insertion code.
- Structures with missing aromatic ring atoms may fail during normal
  calculation; inspect `NaN` results and the printed error for the cause.
- Results are Monte Carlo estimates, and Python's random stream differs from
  MATLAB's.
- The code reports mean displacement from a rate-based hopping model. It does
  not report a FRET efficiency or independently validate the physical
  parameterization.
