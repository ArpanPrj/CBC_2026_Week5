#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${ROOT}"
mkdir -p logs
LOG="${ROOT}/logs/run_all.log"
exec > >(tee "${LOG}") 2>&1
START="$(date +%s)"

echo "============================================================"
echo "FOL4287 FIVE-GENOME SYNTENY PIPELINE"
echo "============================================================"

# Reproducibly bootstrap the software environment if necessary.
bash setup.sh

# Remove all generated scientific data/results so every run starts cleanly.
rm -rf data tmp results
rm -f CHECKSUMS.txt
mkdir -p data/raw data/analysis tmp results

bash scripts/01_download_genomes.sh
bash scripts/02_run_mummer.sh
bash scripts/03_summarize_multigenome.sh
bash scripts/04_plot_summary.sh
bash scripts/05_plot_circos.sh
bash scripts/06_checksums.sh

REQUIRED=(
  "results/genome_metadata.tsv"
  "results/multigenome_chromosome_synteny.tsv"
  "results/percent_aligned_matrix.tsv"
  "results/core_vs_lineage_specific_by_query.tsv"
  "results/chromosome_across_query_summary.tsv"
  "results/Fol4287_multigenome_heatmap.svg"
  "results/core_vs_lineage_specific_by_query.svg"
  "results/Fol4287_vs_Fo47_circos.svg"
  "CHECKSUMS.txt"
)
for file in "${REQUIRED[@]}"; do
  [[ -s "${file}" ]] || { echo "ERROR: required output missing: ${file}"; exit 1; }
done

END="$(date +%s)"
ELAPSED="$((END-START))"
echo "============================================================"
echo "PIPELINE COMPLETE"
echo "Runtime: $((ELAPSED/60)) min $((ELAPSED%60)) sec"
echo "============================================================"
echo "Primary results:"
echo "  results/percent_aligned_matrix.tsv"
echo "  results/core_vs_lineage_specific_by_query.tsv"
echo "  results/Fol4287_multigenome_heatmap.svg"
echo "  results/Fol4287_vs_Fo47_circos.svg"
echo "  CHECKSUMS.txt"
