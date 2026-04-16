#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO_DIR"

THERMALFIST_DIR="${THERMALFIST_DIR:-../ThermalFIST/thermal-fist}"
PARTICLE_LIST="${PARTICLE_LIST:-$THERMALFIST_DIR/input/list/PDG2014/list-withnuclei.dat}"
mkdir -p "$REPO_DIR/.mplconfig"
export MPLCONFIGDIR="$REPO_DIR/.mplconfig"

rm -rf data/thermalfist_fit_points data/thermalfist_figures

python3 scripts/thermalfist_fit/build_fit_inputs.py \
  --input data/first_group_dn_dy_vs_energy.csv \
  --outdir data/thermalfist_inputs \
  --manifest data/thermalfist_inputs/manifest.csv \
  --fit-particles "pi+,pi-,K+,K-,p,pbar,Lambda,Lambda_bar,Xi,Xi_bar"

g++ -std=c++17 -O2 scripts/thermalfist_fit/fit_by_energy.cpp \
  -I"$THERMALFIST_DIR/include" \
  -I"$THERMALFIST_DIR/build/include" \
  -I"$THERMALFIST_DIR/thirdparty/MersenneTwister" \
  -I"$THERMALFIST_DIR/thirdparty/eigen-3.4.0" \
  -L"$THERMALFIST_DIR/build/lib" -lThermalFIST \
  -L"$(root-config --libdir)" -lMinuit2 \
  -Wl,-rpath,"$THERMALFIST_DIR/build/lib" \
  -Wl,-rpath,"$(root-config --libdir)" \
  -o scripts/thermalfist_fit/fit_by_energy

scripts/thermalfist_fit/fit_by_energy \
  data/thermalfist_inputs/manifest.csv \
  "$PARTICLE_LIST" \
  data/thermalfist_fit_results.csv \
  data/thermalfist_fit_points

python3 scripts/thermalfist_fit/plot_fit_results.py \
  --input data/thermalfist_fit_results.csv \
  --outdir data/thermalfist_figures

python3 scripts/thermalfist_fit/plot_yield_comparisons.py \
  --points-dir data/thermalfist_fit_points \
  --outdir data/thermalfist_figures/yield_comparisons

echo "Done. Results: data/thermalfist_fit_results.csv"
echo "Done. Figures: data/thermalfist_figures"
