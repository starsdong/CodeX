# ThermalFIST Fit (Per-Energy BES dN/dy)

This workflow fits **five parameters per collision energy**:
- `T` (chemical freeze-out temperature)
- `mu_B` (baryochemical potential)
- `mu_S` (strangeness chemical potential)
- `gamma_S` (strangeness undersaturation factor)
- `dV/dy` (reported from fitted `R` via `dV/dy = 4*pi*R^3/3`)

## Files
- `build_fit_inputs.py`: converts `data/first_group_dn_dy_vs_energy.csv` into ThermalFIST input files per energy.
- `fit_by_energy.cpp`: performs one ThermalFIST minimization per energy and writes a CSV summary.
- `plot_fit_results.py`: builds final parameter figures from the fit CSV.
- `plot_yield_comparisons.py`: creates per-energy bar plots of data yields vs fitted ThermalFIST yields.
- `run_all.sh`: one-command wrapper to build inputs, compile, and run all fits.

## 1) Build fit inputs
Run from repository root:

```bash
python3 scripts/thermalfist_fit/build_fit_inputs.py \
  --input data/first_group_dn_dy_vs_energy.csv \
  --outdir data/thermalfist_inputs \
  --manifest data/thermalfist_inputs/manifest.csv \
  --lambda-to-p-br 0.639
```

## 2) Compile fitter
This expects ThermalFIST source/build in `../ThermalFIST/thermal-fist` and ROOT/Minuit2 available via `root-config`.

```bash
g++ -std=c++17 -O2 scripts/thermalfist_fit/fit_by_energy.cpp \
  -I../ThermalFIST/thermal-fist/include \
  -I../ThermalFIST/thermal-fist/build/include \
  -I../ThermalFIST/thermal-fist/thirdparty/MersenneTwister \
  -I../ThermalFIST/thermal-fist/thirdparty/eigen-3.4.0 \
  -L../ThermalFIST/thermal-fist/build/lib -lThermalFIST \
  -L"$(root-config --libdir)" -lMinuit2 \
  -Wl,-rpath,../ThermalFIST/thermal-fist/build/lib \
  -Wl,-rpath,"$(root-config --libdir)" \
  -o scripts/thermalfist_fit/fit_by_energy
```

## 3) Run fits

```bash
scripts/thermalfist_fit/fit_by_energy \
  data/thermalfist_inputs/manifest.csv \
  ../ThermalFIST/thermal-fist/input/list/PDG2014/list-withnuclei.dat \
  data/thermalfist_fit_results.csv \
  data/thermalfist_fit_points
```

Or run the full chain:

```bash
scripts/thermalfist_fit/run_all.sh
```

## Output
`data/thermalfist_fit_results.csv` columns include:
`energy_GeV, T_MeV, muB_MeV, muS_MeV, gammaS, R_fm, dVdy_fm3, chi2, ndf, chi2_ndf`.

`data/thermalfist_figures/` includes:
- `T_vs_energy.png`
- `muB_vs_energy.png`
- `muS_vs_energy.png`
- `gammaS_vs_energy.png`
- `dVdy_vs_energy.png`
- `T_vs_muB.png`

`data/thermalfist_figures/yield_comparisons/` includes one plot per collision energy:
- `sqrts_*GeV_yield_compare.png` (bar comparison: data vs fitted model for all fitted particle yields at that energy)

## Notes
- Only rows with observable containing `dN/dy` are used.
- Unsupported species (e.g. proxy rows and `Omega_sum`) are skipped.
- Total data uncertainty is combined as `sqrt(stat_err^2 + sys_err^2)`.
- Proton feeddown correction is applied before fitting:
  `p -> p - BR*Lambda` and `pbar -> pbar - BR*Lambda_bar`, with propagated uncertainty.
  Default `BR = 0.639` (configurable via `--lambda-to-p-br`).
- The C++ fitter parses the generated input text directly (instead of relying on ThermalFIST data-file loader).
