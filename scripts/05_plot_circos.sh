#!/usr/bin/env bash

# =============================================================================
# STEP 5 SHELL WRAPPER: DETAILED FOL4287-vs-FO47 CIRCOS FIGURE
# =============================================================================
# The multi-genome result is summarized in Step 4.  This separate figure keeps
# the detailed chromosome-to-chromosome ribbons limited to one biologically
# interpretable pair so the plot does not become unreadably crowded.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"
source "${ROOT}/config/parameters.sh"

# Export the display thresholds from parameters.sh so the R plotting script can
# reproduce the filtering annotation in the figure subtitle.
export PLOT_MIN_ALIGNMENT_LENGTH
export PLOT_MIN_IDENTITY

mkdir -p results logs
LOG="${ROOT}/logs/05_plot_circos.log"
exec > >(tee "${LOG}") 2>&1

echo "============================================================"
echo "STEP 5: DETAILED FOL4287-vs-FO47 CIRCOS FIGURE"
echo "============================================================"

# Generate the vector SVG using circlize and the filtered link table from Step 3.
Rscript "${SCRIPT_DIR}/05_plot_circos.R"
[[ -s results/Fol4287_vs_Fo47_circos.svg ]] || { echo "ERROR: Circos SVG missing"; exit 1; }

echo "STEP 5 COMPLETE"
