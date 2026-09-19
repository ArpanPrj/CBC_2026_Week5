#!/usr/bin/env bash

# =============================================================================
# STEP 2: SEQUENTIAL PAIRWISE MUMMER/DNADIFF COMPARISONS
# =============================================================================
# Reference
# ---------
# The 15-record Fol4287 chromosome-only FASTA generated in Step 1.
#
# Queries
# -------
# Every row in config/genomes.tsv whose role column is "query".  The query
# assemblies are the complete downloaded FASTAs, so the reference chromosome
# sequence can align to chromosome or scaffold sequence present in each query.
#
# Why dnadiff?
# ------------
# dnadiff is MUMmer's two-genome comparison workflow.  It runs nucmer with
# max-match anchoring and derives filtered one-to-one comparisons among its
# outputs.  This project retains the human-readable .report and .1coords files
# for every pair; the larger intermediate delta files remain in tmp/.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"
source "${ROOT}/config/parameters.sh"

mkdir -p tmp/mummer results/reports results/pairwise_coords logs
LOG="${ROOT}/logs/02_run_mummer.log"
exec > >(tee "${LOG}") 2>&1

# Construct the stable path to the prepared Fol4287 chromosome reference using
# the REFERENCE_NAME configured in parameters.sh.
REFERENCE="data/analysis/${REFERENCE_NAME}.chromosomes.fna"
if [[ ! -s "${REFERENCE}" ]]; then
  echo "ERROR: missing reference chromosomes: ${REFERENCE}"
  echo "Run Step 1 first."
  exit 1
fi

# Delete previous MUMmer intermediates before beginning the new set of pairwise
# comparisons.  This guarantees that a stale .delta file cannot be mistaken for
# output from the current run.
rm -rf tmp/mummer
mkdir -p tmp/mummer
rm -rf results/reports results/pairwise_coords
mkdir -p results/reports results/pairwise_coords

# Initialize an index table that records, for each query, its accession/context
# and the two retained output paths.
printf 'query\taccession\thost_context\treport\tcoords\n' > results/pairwise_outputs.tsv

echo "============================================================"
echo "STEP 2: SEQUENTIAL MUMMER DNADIFF COMPARISONS"
echo "Reference: ${REFERENCE_NAME}"
echo "============================================================"
echo "dnadiff uses nucmer --maxmatch internally and then creates a 1-to-1 filtered alignment."

# Iterate over genomes.tsv in its tracked order.  Only rows explicitly marked
# role=query proceed to MUMmer, so the Fol4287 reference row is skipped.
while IFS=$'\t' read -r name accession role host_context circos; do
  [[ "${name}" == "name" ]] && continue
  [[ "${role}" != "query" ]] && continue

  # Each query is aligned using the complete raw NCBI FASTA downloaded in Step 1.
  query="data/raw/${name}.fna"
  # dnadiff creates several related files sharing this temporary prefix.
  prefix="tmp/mummer/${REFERENCE_NAME}_vs_${name}"
  # Capture a per-query start time for informative runtime logging only.
  start="$(date +%s)"

  if [[ ! -s "${query}" ]]; then
    echo "ERROR: missing query FASTA: ${query}"
    exit 1
  fi

  echo
  echo "------------------------------------------------------------"
  echo "Aligning ${REFERENCE_NAME} vs ${name} (${accession})"
  echo "------------------------------------------------------------"

  # Main external computation.  -p supplies the output prefix, followed by the
  # reference FASTA and query FASTA in that order.
  dnadiff -p "${prefix}" "${REFERENCE}" "${query}"

  # Validate the three MUMmer products this project relies on.  `.1coords` is the
  # one-to-one coordinate table parsed in Step 3; `.report` is retained as a
  # standard MUMmer summary; `.1delta` confirms creation of the filtered delta.
  for required in "${prefix}.report" "${prefix}.1coords" "${prefix}.1delta"; do
    if [[ ! -s "${required}" ]]; then
      echo "ERROR: expected MUMmer output missing: ${required}"
      exit 1
    fi
  done

  # Stable result paths are separated from temporary MUMmer internals.
  report_out="results/reports/${REFERENCE_NAME}_vs_${name}.report"
  coords_out="results/pairwise_coords/${REFERENCE_NAME}_vs_${name}.1coords.tsv"
  # Copy the small report into results/ while leaving the temporary original.
  cp "${prefix}.report" "${report_out}"
  # Copy the one-to-one coordinate table that drives all quantitative summaries.
  cp "${prefix}.1coords" "${coords_out}"

  # Append one deterministic row to the pairwise output index for this completed
  # comparison.
  printf '%s\t%s\t%s\t%s\t%s\n' \
    "${name}" "${accession}" "${host_context}" "${report_out}" "${coords_out}" \
    >> results/pairwise_outputs.tsv

  # Log this individual comparison's elapsed seconds; this value is not used in
  # biological calculations or checksummed tables.
  end="$(date +%s)"
  echo "Completed ${name} in $((end-start)) seconds."
done < config/genomes.tsv

echo
echo "STEP 2 COMPLETE"
