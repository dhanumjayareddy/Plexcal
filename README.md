# Plexcal: excitation energy-transfer diffusion lengths

Plexcal is a structure-aware scientific workspace for estimating excitation
energy-transfer diffusion lengths with the existing FRET and DET models.
Uploaded, RCSB, and AlphaFold PDB structures can be explored in the full
interactive PDBe Mol* viewer. The local analysis server does not build or query
the former SQLite structure-record database.

## Run the local application

Install the calculation dependencies and start the server:

```text
python -m pip install -r requirements.txt
python server.py
```

Open <http://127.0.0.1:8000>. The structure workspace loads an uploaded PDB,
an RCSB PDB ID, or an AlphaFold structure ID into the PDBe Mol* viewer.
Uploads are previewed as soon as they are selected; ID-based previews can be
loaded independently of calculation. The viewer keeps its Mol* sequence,
selection, component, representation, and display controls available. The
visibility panel can independently toggle protein polymers, ligands and other
heteroatoms, water, carbohydrates, non-standard residues, and density maps;
native Mol* controls also provide structure rotation and representation options.

RCSB structures include deposited structure details, organism, experimental
method, resolution, and the primary citation when supplied by RCSB. AlphaFold
models include protein and organism annotations, model confidence and dates,
and a citation for the AlphaFold Protein Structure Database. Uploaded PDB
headers are parsed locally for available title, organism, method, and citation
details. Metadata retrieval does not modify the structure calculations.

Run FRET and DET from the analysis section. Both methods use their existing
calculation functions and the same local PDB file, whether it was uploaded or
downloaded. Uploads and downloads are limited to 25 MiB; temporary calculation
files are removed afterward. The results include excitation energy-transfer
diffusion length, total pairs, coherent pairs, and coherent-pair fraction for
each method, with separate comparison charts for lengths, pair counts, and
fractions.

Completed results are retained in the current browser tab across page refreshes.
Uploaded PDB files are not persisted by the browser, so select the file again
to restore its 3D preview after refreshing.

The uploaded-file FRET calculator is a Python translation of
`process_single_pdb_v1.m`. Its Monte Carlo calculation uses 10,000 excitation
trajectories with up to 1,000 transitions per trajectory. Because Python and
MATLAB use different random-number streams, the result is not guaranteed to
match MATLAB's result bit-for-bit.

The separate distance-exponential DET calculation translated from
`det_code_clean.m` is available as `det_code_clean.py`. Run it on one structure
with `python det_code_clean.py "1TUB.pdb" 0`; see
[`DET_CODE_CLEAN.md`](DET_CODE_CLEAN.md) for its workflow and equations.

## Batch processing

To calculate diffusion lengths and pair statistics for the structures in
`all_proteins_pdb/`, run:

```text
python process_all_pdbs.py
```

The script uses up to four worker processes by default and writes
`all_proteins_pdb/diffusion_lengths.csv`. To choose another worker count, pass
`--workers`, for example `python process_all_pdbs.py --workers 2`.

For details of the single-PDB FRET workflow, input interpretation, physical
model, equations, and limitations, see
[`PROCESS_SINGLE_PDB.md`](PROCESS_SINGLE_PDB.md).
