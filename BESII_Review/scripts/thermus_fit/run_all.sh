#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO_DIR"

THERMUS_DIR="${THERMUS_DIR:-$REPO_DIR/THERMUS/source}"
PARTICLE_LIST="${PARTICLE_LIST:-$THERMUS_DIR/particles/PartList_PPB2014_CBHN.txt}"
DECAYS_DIR="${DECAYS_DIR:-$THERMUS_DIR/particles}"
mkdir -p "$REPO_DIR/.mplconfig"
export MPLCONFIGDIR="$REPO_DIR/.mplconfig"

rm -rf data/thermus_fit_points data/thermus_figures

python3 scripts/thermus_fit/build_fit_inputs.py \
  --input data/first_group_dn_dy_vs_energy.csv \
  --outdir data/thermus_inputs \
  --manifest data/thermus_inputs/manifest.csv \
  --fit-particles "pi+,pi-,K+,K-,p,pbar,Lambda,Lambda_bar,Xi,Xi_bar,phi" \
  --energies "7.7,11.5,19.6,27,39,62.4,130,200" \
  --exclude-energy-particles "62.4:Xi,Xi_bar"

g++ -std=c++17 -O2 scripts/thermus_fit/fit_by_energy.cpp \
  -I"$THERMUS_DIR/main" \
  -I"$THERMUS_DIR/include" \
  -I"$(root-config --incdir)" \
  $(root-config --cflags) \
  -L"$THERMUS_DIR/build/lib" -lTHERMUS -lFunctions \
  $(root-config --libs) \
  -Wl,-rpath,"$THERMUS_DIR/build/lib" \
  -Wl,-rpath,"$(root-config --libdir)" \
  -o scripts/thermus_fit/fit_by_energy

scripts/thermus_fit/fit_by_energy \
  data/thermus_inputs/manifest.csv \
  "$PARTICLE_LIST" \
  "$DECAYS_DIR" \
  data/thermus_fit_results.csv \
  data/thermus_fit_points \
  > data/thermus_fit.log 2>&1

python3 scripts/thermus_fit/plot_fit_results.py \
  --input data/thermus_fit_results.csv \
  --outdir data/thermus_figures

python3 scripts/thermus_fit/plot_yield_comparisons.py \
  --points-dir data/thermus_fit_points \
  --outdir data/thermus_figures/yield_comparisons \
  --fit-summary data/thermus_fit_results.csv

echo "Done. Results: data/thermus_fit_results.csv"
echo "Done. Figures: data/thermus_figures"
