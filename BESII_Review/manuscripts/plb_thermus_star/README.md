# PLB THERMUS STAR Draft

This directory contains a Physics Letters B-style manuscript draft based on the
THERMUS fits in this repository.

## Files

- `main.tex`: manuscript source using `elsarticle`.
- `references.bib`: bibliography used by the draft.
- `highlights.md`: submission highlights.
- `figures/freezeout_summary.pdf`: generated freeze-out parameter figure.

The manuscript also includes figures from `../../data/`:

- `strange_hadron_model_comparison_triptych.pdf`
- `thermus_sce_particle_set_comparison.pdf`

## Build

From this directory:

```bash
pdflatex main
bibtex main
pdflatex main
pdflatex main
```

Regenerate the freeze-out summary figure from the repository root:

```bash
python3 scripts/prepare_plb_thermus_article_figures.py
```

## Notes Before Submission

- Replace placeholder author names and affiliations.
- Verify all reference metadata and journal abbreviations.
- Recheck whether the 62.4 GeV point should be shown in the main table or moved
  to supplementary material because of its reduced input set.
- Replace the generative-AI declaration with the final author-approved statement.
