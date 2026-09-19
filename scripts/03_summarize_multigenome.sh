#!/usr/bin/env bash

# =============================================================================
# STEP 3 SHELL WRAPPER: MULTI-GENOME SYNTENY SUMMARIZATION
# =============================================================================
# The substantive calculations are implemented in 03_summarize_multigenome.py.
# This wrapper supplies the shared Conda environment, exports the configured
# visualization parameters so Python can read them, captures a log, and verifies
# that all expected summary tables were created successfully.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"
source "${ROOT}/config/parameters.sh"

# Export shell variables because Python reads these values via os.environ rather
# than sourcing the shell configuration directly.
export PLOT_MIN_ALIGNMENT_LENGTH
export PLOT_MIN_IDENTITY
export CIRCOS_QUERY

mkdir -p results logs
LOG="${ROOT}/logs/03_summarize_multigenome.log"
exec > >(tee "${LOG}") 2>&1

echo "============================================================"
echo "STEP 3: MULTI-GENOME CHROMOSOME SYNTENY SUMMARY"
echo "============================================================"

# Execute the deterministic parser/summarizer against all .1coords files produced
# in Step 2.
python "${SCRIPT_DIR}/03_summarize_multigenome.py"

# Validate the central long-form table, heatmap matrix, class summary,
# across-query chromosome summary, and filtered Fo47 Circos-link table.
for file in \
  results/multigenome_chromosome_synteny.tsv \
  results/percent_aligned_matrix.tsv \
  results/core_vs_lineage_specific_by_query.tsv \
  results/chromosome_across_query_summary.tsv \
  results/fo47_circos_links.tsv; do
  [[ -s "${file}" ]] || { echo "ERROR: missing output ${file}"; exit 1; }
done

echo "STEP 3 COMPLETE"
