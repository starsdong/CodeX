#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO_DIR"

THERMUS_DIR="${THERMUS_DIR:-$REPO_DIR/THERMUS/source}"
PARTICLE_LIST="${PARTICLE_LIST:-$THERMUS_DIR/particles/PartList_PPB2014_CBHN.txt}"
DECAYS_DIR="${DECAYS_DIR:-$THERMUS_DIR/particles}"
OUT_DIR="${OUT_DIR:-data/thermus_sce_3gev}"
INPUT_DIR="$OUT_DIR/input"
POINTS_DIR="$OUT_DIR/fit_points"
LAMBDA_TO_P_BR="${LAMBDA_TO_P_BR:-0.639}"

mkdir -p "$INPUT_DIR" "$POINTS_DIR"
mkdir -p "$REPO_DIR/.mplconfig"
export MPLCONFIGDIR="$REPO_DIR/.mplconfig"

python3 scripts/thermus_fit/build_fit_inputs.py \
  --input data/first_group_dn_dy_vs_energy.csv \
  --outdir "$INPUT_DIR" \
  --manifest "$INPUT_DIR/manifest.csv" \
  --min-points 1 \
  --energies "3.0" \
  --lambda-to-p-br "$LAMBDA_TO_P_BR" \
  --fit-particles "K-,Ks0,phi,p,Lambda,Xi"

g++ -std=c++17 -O2 scripts/thermus_fit/fit_3gev_sce.cpp \
  -I"$THERMUS_DIR/main" \
  -I"$THERMUS_DIR/include" \
  -I"$(root-config --incdir)" \
  $(root-config --cflags) \
  -L"$THERMUS_DIR/build/lib" -lTHERMUS -lFunctions \
  $(root-config --libs) \
  -Wl,-rpath,"$THERMUS_DIR/build/lib" \
  -Wl,-rpath,"$(root-config --libdir)" \
  -o scripts/thermus_fit/fit_3gev_sce

scripts/thermus_fit/fit_3gev_sce \
  "$INPUT_DIR/sqrts_3GeV.txt" \
  "$PARTICLE_LIST" \
  "$DECAYS_DIR" \
  "$OUT_DIR/fit_summary.csv" \
  "$POINTS_DIR/sqrts_3GeV_points.csv" \
  > "$OUT_DIR/fit.log" 2>&1

python3 scripts/thermus_fit/plot_yield_comparisons.py \
  --points-dir "$POINTS_DIR" \
  --outdir "$OUT_DIR/yield_comparisons" \
  --fit-summary "$OUT_DIR/fit_summary.csv"

echo "Done. Results: $OUT_DIR/fit_summary.csv"
echo "Done. Points:  $POINTS_DIR/sqrts_3GeV_points.csv"
echo "Done. Figure:  $OUT_DIR/yield_comparisons/sqrts_3GeV_yield_compare.png"
