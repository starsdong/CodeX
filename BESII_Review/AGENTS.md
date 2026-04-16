# Repository Guidelines

## Project Structure & Module Organization
- `scripts/`: standalone Python data-processing and plotting scripts.
- `data/`: canonical CSV datasets and generated figures.
- `data/particle_figures/`: per-particle log-log plots plus `manifest.csv`.
- `BESPapers.docx`: source bibliography document used by extraction/enrichment scripts.
- `.mplconfig/`: local Matplotlib cache/config for reproducible rendering in this workspace.

Keep generated outputs in `data/` and avoid writing artifacts outside the repo.

## Build, Test, and Development Commands
Use Python 3 from the repository root:

```bash
python3 scripts/build_first_group_plot.py
python3 scripts/plot_first_group_dn_dy_only.py
python3 scripts/split_dn_dy_by_particle.py
python3 scripts/extract_star_links.py --input data/publications_from_docx.csv --output data/publications_enriched.csv
```

- `build_first_group_plot.py`: fetches external HEPData JSON and regenerates the core yield dataset/plot.
- `plot_first_group_dn_dy_only.py`: filters to `dN/dy` rows and writes `first_group_dn_dy_vs_energy.*`.
- `split_dn_dy_by_particle.py`: generates per-particle PNGs and a manifest.
- `extract_star_links.py`: enriches publication metadata from STAR pages.

Note: some scripts require network access; `add_table8_0808_2041.py` also requires `pdftotext`.

## Coding Style & Naming Conventions
- Follow existing style: 4-space indentation, top-level script constants, small helper functions.
- Use `snake_case` for functions/variables and descriptive script names (e.g., `repair_core_dataset.py`).
- Prefer `pathlib.Path`, `csv.DictReader/DictWriter`, and explicit UTF-8 file handling.
- Preserve stable CSV column names and ordering when updating datasets.

## Testing Guidelines
There is no formal test framework in this repository yet. Validate changes by:
- Running the affected script(s) end-to-end.
- Confirming regenerated CSVs load and row counts/keys are sensible.
- Verifying expected PNGs are produced in `data/` or `data/particle_figures/`.

For new logic, add small deterministic checks in-script or add a focused validation script under `scripts/`.

## Commit & Pull Request Guidelines
- Match current commit style: imperative, concise subject lines (e.g., `Keep only offset compare plots and update constrained fit window`).
- Keep commits focused: script logic, data refresh, and plot outputs should be grouped coherently.
- PRs should include:
  - what script(s) changed,
  - what data/figures were regenerated,
  - external sources touched (HEPData/STAR/arXiv),
  - before/after notes for any dataset schema or centrality/observable interpretation changes.
