# PDB Diffusion Database

The project contains two ways to serve the same UI and CSV records:

- **GitHub Pages (static, no backend):** build and publish `docs/`. The site
  fetches a small index when it opens and downloads only the data shard for a
  searched PDB ID. IDs are grouped by their first two characters; the build
  checks that no shard exceeds 900 KB. The generated shards total about
  110 MiB to store in the repository, but a visitor does not download all of
  them: the index is about 2.3 KB and a search downloads just one shard.
- **Local API (optional):** run `python server.py` to serve lookups from a
  generated SQLite database.

## Build and test the GitHub Pages site

From this folder, run:

```text
python build_pages.py
```

This copies the current HTML, CSS, and static-site JavaScript into `docs/`,
then generates the index and JSON shards from
`RCSB_MASTER_DATA_EXPANDED_ANNOTATED.csv`. Re-run it whenever the UI or source
CSV changes. The original CSV must be present alongside `build_pages.py` while
building; GitHub Pages serves only the generated files in `docs/`.

To check the static site locally, run:

```text
python -m http.server 8001 --directory docs
```

Then open <http://127.0.0.1:8001>.

## Publish with GitHub Pages

1. Build the site and commit the generated `docs/` directory.
2. Push the repository to GitHub.
3. In **Settings → Pages**, choose **Deploy from a branch**, select the branch
   (for example, `main`) and the `/docs` folder, then save.
4. Open the Pages URL shown in the Settings page.

The older `pdb_database.json` files and the local SQLite database are generated
artifacts and are not needed for GitHub Pages; `.gitignore` excludes them from
new commits. If any of these files are already tracked by Git, remove them from
Git tracking before pushing: each generated database file exceeds GitHub's
100 MB per-file limit.

## Run the optional local API

Install the calculation dependencies and start the app:

```text
python -m pip install -r requirements.txt
python server.py
```

Open <http://127.0.0.1:8000>. The first launch builds a SQLite index in `data/`;
later launches reuse it and rebuild it if the source CSV changes. In the local
calculation form, provide a PDB upload, an RCSB four-character PDB ID, or an
AlphaFold/UniProt structure ID. Both models run on the resulting local PDB file.
Uploads and downloads are limited to 25 MiB and temporary files are removed
after calculation. GitHub Pages supports database lookups only; it cannot run
structure calculations.

The uploaded-file calculator is a Python translation of
`process_single_pdb_v1.m`. Its Monte Carlo calculation uses 10,000 excitation
trajectories with up to 1,000 transitions per trajectory. Because Python and
MATLAB use different random-number streams, the result is not guaranteed to
match MATLAB's result bit-for-bit.

## Process a folder of PDB files in parallel

To calculate diffusion lengths for the PDB files in `all_proteins_pdb/` using
the same calculation as the single-file processor, run:

```text
python process_all_pdbs.py
```

The script uses up to four worker processes by default and writes
`all_proteins_pdb/diffusion_lengths.csv` with each PDB filename, diffusion
length in Angstroms, total chromophore pairs, coherent pairs, and coherent-pair
fraction. To choose another worker count, pass
`--workers`, for example `python process_all_pdbs.py --workers 2`. The existing
diffusion-length calculation remains unchanged.

For a detailed explanation of the single-PDB workflow, input interpretation,
physical model, equations, and limitations, see
[PROCESS_SINGLE_PDB.md](PROCESS_SINGLE_PDB.md).

The separate distance-exponential DET calculation translated from
`det_code_clean.m` is available as `det_code_clean.py`. Run it on one structure
with `python det_code_clean.py "1TUB.pdb" 0`; see
[DET_CODE_CLEAN.md](DET_CODE_CLEAN.md) for its workflow and equations. It does
not replace `process_single_pdb.py`. The local analysis form accepts an uploaded
PDB, a four-character RCSB PDB ID, or a UniProt/AlphaFold structure ID. It runs
both existing calculation functions on the resulting local PDB and displays
the diffusion lengths, pair statistics, and unit-specific comparison charts.
