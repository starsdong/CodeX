#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO_DIR"

THERMUS_DIR="${THERMUS_DIR:-$REPO_DIR/THERMUS/source}"
PARTICLE_LIST="${PARTICLE_LIST:-$THERMUS_DIR/particles/PartList_PPB2014_CBHN.txt}"
DECAYS_DIR="${DECAYS_DIR:-$THERMUS_DIR/particles}"
OUT_DIR="${OUT_DIR:-data/thermus_sce_blended}"
MODEL_LABEL="${MODEL_LABEL:-THERMUS SCE (endpoint-interpolated)}"
TITLE="${TITLE:-Strange-Hadron Yields vs Energy with THERMUS SCE (endpoint-interpolated)}"
WEIGHT_POWER="${WEIGHT_POWER:-2.0}"
MATCH_REF_DIR="${MATCH_REF_DIR:-data/thermus_gc_from_sce_effective/prediction_points}"
MATCH_PARTICLES="${MATCH_PARTICLES:-p phi}"

mkdir -p "$OUT_DIR"
mkdir -p "$REPO_DIR/.mplconfig"
export MPLCONFIGDIR="$REPO_DIR/.mplconfig"

python3 scripts/thermus_fit/build_sce_blended_endpoint.py

g++ -std=c++17 -O2 scripts/thermus_fit/predict_sce_grid.cpp \
  -I"$THERMUS_DIR/main" \
  -I"$THERMUS_DIR/include" \
  -I"$(root-config --incdir)" \
  $(root-config --cflags) \
  -L"$THERMUS_DIR/build/lib" -lTHERMUS -lFunctions \
  $(root-config --libs) \
  -Wl,-rpath,"$THERMUS_DIR/build/lib" \
  -Wl,-rpath,"$(root-config --libdir)" \
  -o scripts/thermus_fit/predict_sce_grid

scripts/thermus_fit/predict_sce_grid \
  "$OUT_DIR/endpoint_7p7/prediction_manifest.csv" \
  "$PARTICLE_LIST" \
  "$DECAYS_DIR" \
  "$OUT_DIR/endpoint_7p7/prediction_summary.csv" \
  "$OUT_DIR/endpoint_7p7/prediction_points" \
  > "$OUT_DIR/endpoint_7p7/predict.log" 2>&1

blend_args=(
  --outdir "$OUT_DIR"
  --weight-power "$WEIGHT_POWER"
)
if [[ -n "$MATCH_REF_DIR" && -n "$MATCH_PARTICLES" ]]; then
  blend_args+=(--match-ref-dir "$MATCH_REF_DIR" --match-particles $MATCH_PARTICLES)
fi

python3 scripts/thermus_fit/build_sce_blended_predictions.py "${blend_args[@]}"

python3 scripts/thermus_fit/plot_sce_low_energy_predictions.py \
  --model-dir "$OUT_DIR/prediction_points" \
  --outdir "$OUT_DIR" \
  --title "$TITLE" \
  --model-label "$MODEL_LABEL"

echo "Done. Endpoints: $OUT_DIR/endpoint_parameters.csv"
echo "Done. Weights:   $OUT_DIR/interpolation_grid.csv"
echo "Done. Summary:   $OUT_DIR/prediction_summary.csv"
echo "Done. Figure:    $OUT_DIR/strange_hadron_yields_with_sce.png"
