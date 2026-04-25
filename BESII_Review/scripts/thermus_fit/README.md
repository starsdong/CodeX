# THERMUS Fit (Per-Energy BES dN/dy)

This workflow mirrors the ThermalFIST setup, but runs fits with THERMUS and writes to separate outputs.

## Fitted Parameters Per Energy
- `T`
- `muB`
- `muS`
- `R` (reported also as `dV/dy = 4*pi*R^3/3`)

`gammaS` is fixed to `1`. `130 GeV` is skipped.

## Inputs and Particle Set
- Source data: `data/first_group_dn_dy_vs_energy.csv`
- Fitted particles: `pi+, pi-, K+, K-, p, pbar, Lambda, Lambda_bar, Xi, Xi_bar, phi`
- Proton feeddown correction: `p -> p - BR*Lambda`, `pbar -> pbar - BR*Lambda_bar`, with `BR=0.639`

## Run
From repo root:

```bash
scripts/thermus_fit/run_all.sh
```

Optional environment overrides:
- `THERMUS_DIR` (default: `THERMUS/source`)
- `PARTICLE_LIST` (default: `THERMUS/source/particles/PartList_PPB2014_CBHN.txt`)
- `DECAYS_DIR` (default: `THERMUS/source/particles`)

## Outputs
- `data/thermus_inputs/`
- `data/thermus_fit_results.csv`
- `data/thermus_fit_points/`
- `data/thermus_figures/`
  - `T_vs_energy.png`, `muB_vs_energy.png`, `muS_vs_energy.png`, `gammaS_vs_energy.png`, `dVdy_vs_energy.png`, `T_vs_muB.png`
  - `yield_comparisons/sqrts_*GeV_yield_compare.png`
