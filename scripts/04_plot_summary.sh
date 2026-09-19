#!/usr/bin/env bash

# =============================================================================
# STEP 4 SHELL WRAPPER: MULTI-GENOME SUMMARY FIGURES
# =============================================================================
# This wrapper activates the common environment/logging context, runs the R
# plotting script, and verifies that both SVG figures were produced and nonempty.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

mkdir -p results logs
LOG="${ROOT}/logs/04_plot_summary.log"
exec > >(tee "${LOG}") 2>&1

echo "============================================================"
echo "STEP 4: MULTI-GENOME SUMMARY FIGURES"
echo "============================================================"

# Run the base-R plotting implementation against Step 3's deterministic tables.
Rscript "${SCRIPT_DIR}/04_plot_summary.R"

# Treat a missing or zero-byte figure as a failed stage rather than silently
# continuing to later steps.
for file in results/Fol4287_multigenome_heatmap.svg results/core_vs_lineage_specific_by_query.svg; do
  [[ -s "${file}" ]] || { echo "ERROR: missing figure ${file}"; exit 1; }
done

echo "STEP 4 COMPLETE"
