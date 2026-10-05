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

Run `python server.py` and open <http://127.0.0.1:8000>. The first launch builds
a SQLite index in `data/`; later launches reuse it and rebuild it if the source
CSV changes. No third-party Python packages are required.
