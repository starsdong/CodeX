#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO_DIR"

THERMUS_DIR="${THERMUS_DIR:-$REPO_DIR/THERMUS/source}"
PARTICLE_LIST="${PARTICLE_LIST:-$THERMUS_DIR/particles/PartList_PPB2014_CBHN.txt}"
DECAYS_DIR="${DECAYS_DIR:-$THERMUS_DIR/particles}"
OUT_DIR="${OUT_DIR:-data/thermus_fit_predictions_full}"

mkdir -p "$OUT_DIR"

python3 scripts/thermus_fit/build_gce_fit_full_predictions.py

g++ -std=c++17 -O2 scripts/thermus_fit/predict_gce_grid.cpp \
  -I"$THERMUS_DIR/main" \
  -I"$THERMUS_DIR/include" \
  -I"$(root-config --incdir)" \
  $(root-config --cflags) \
  -L"$THERMUS_DIR/build/lib" -lTHERMUS -lFunctions \
  $(root-config --libs) \
  -Wl,-rpath,"$THERMUS_DIR/build/lib" \
  -Wl,-rpath,"$(root-config --libdir)" \
  -o scripts/thermus_fit/predict_gce_grid

scripts/thermus_fit/predict_gce_grid \
  "$OUT_DIR/prediction_manifest.csv" \
  "$PARTICLE_LIST" \
  "$DECAYS_DIR" \
  "$OUT_DIR/prediction_summary.csv" \
  "$OUT_DIR/prediction_points" \
  > "$OUT_DIR/predict.log" 2>&1

echo "Done. Parameters: $OUT_DIR/parameter_grid.csv"
echo "Done. Summary:    $OUT_DIR/prediction_summary.csv"
echo "Done. Points:     $OUT_DIR/prediction_points"
