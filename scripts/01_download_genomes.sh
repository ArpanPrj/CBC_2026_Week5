#!/usr/bin/env bash

# =============================================================================
# STEP 1: DOWNLOAD EXACT ASSEMBLIES AND PREPARE ANALYSIS INPUTS
# =============================================================================
# Inputs
# ------
# config/genomes.tsv
#   The tracked table of short names, versioned NCBI accessions, biological
#   context labels, and the flag identifying the detailed Circos query.
#
# Main external programs
# ----------------------
# * NCBI Datasets: downloads each exact assembly accession.
# * unzip: expands each NCBI data package.
# * SeqKit: records broad FASTA assembly statistics.
# * Python helper 01_parse_genomes.py: parses FASTA records and creates the
#   chromosome-level reference/Fo47 inputs and metadata tables.
#
# Key outputs
# -----------
# data/raw/<name>.fna                 full downloaded genomic FASTAs
# data/analysis/Fol4287.chromosomes.fna
# data/analysis/Fo47.chromosomes.fna  chromosome-only FASTAs prepared by Python
# results/genome_stats.tsv            SeqKit statistics for all raw genomes
# results/genome_metadata.tsv         parsed record counts and assembly sizes
# =============================================================================
# Fail immediately on errors, unset variables, or hidden failures in pipelines.
set -euo pipefail

# Resolve this file's directory so common.sh can be sourced correctly even when
# the user launches the script from somewhere other than scripts/.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# common.sh moves to the repository root, fixes locale/timezone, and activates
# the already-created project Conda environment.
source "${SCRIPT_DIR}/common.sh"

# Create all directories this stage may write.  `-p` makes the operation safe on
# reruns when some or all directories already exist.
mkdir -p data/raw data/analysis tmp/downloads results logs
LOG="${ROOT}/logs/01_download_genomes.log"
exec > >(tee "${LOG}") 2>&1

# The genome panel is data-driven: adding/removing rows in genomes.tsv changes the
# set downloaded without requiring hard-coded per-genome commands in this script.
CONFIG="config/genomes.tsv"

echo "============================================================"
echo "STEP 1: DOWNLOAD EXACT VERSIONED GENOME ASSEMBLIES"
echo "============================================================"

# Bash array used to collect every final raw FASTA path.  The complete array is
# passed to one SeqKit command after all downloads finish.
FASTA_FILES=()

# -----------------------------------------------------------------------------
# LOOP OVER CONFIGURED ASSEMBLIES
# -----------------------------------------------------------------------------
# IFS is set to a literal tab because genomes.tsv is tab-separated.  `read -r`
# avoids backslash interpretation.  Each column is assigned to a named variable
# matching the header fields.
# Read the versioned assembly accessions from the tracked config file.
while IFS=$'\t' read -r name accession role host_context circos; do
  # Skip the TSV header row rather than treating the word "name" as a genome.
  [[ "${name}" == "name" ]] && continue
  # Ignore an accidental blank record safely.
  [[ -z "${name}" ]] && continue

  # Use a genome-specific temporary directory so each downloaded NCBI package is
  # isolated from every other assembly.
  download_dir="tmp/downloads/${name}"
  # Store the NCBI ZIP under a filename containing the exact accession.
  archive="${download_dir}/${accession}.zip"
  # NCBI Datasets ZIP contents are expanded into this package directory.
  package_dir="${download_dir}/package"
  # The FASTA used by later analysis receives a stable short-name path independent
  # of NCBI's internal package directory structure.
  final_fasta="data/raw/${name}.fna"

  echo
  echo "Downloading ${name} (${accession})"
  # Remove any prior temporary package for this genome so a rerun cannot reuse
  # stale extracted contents.
  rm -rf "${download_dir}"
  mkdir -p "${download_dir}" "${package_dir}"

  # Ask NCBI Datasets for the exact accession version and include only the genome
  # sequence payload needed for this project.  The backslash simply continues
  # the same shell command onto the following lines.
  datasets download genome accession "${accession}" \
    --include genome \
    --filename "${archive}"

  # Quietly extract the downloaded NCBI archive into the isolated package folder.
  unzip -q "${archive}" -d "${package_dir}"

  # Count genomic FASTA files inside the package.  The workflow requires exactly
  # one .fna file so an unexpected package layout triggers an explicit error.
  count="$(find "${package_dir}/ncbi_dataset/data" -type f -name '*.fna' | wc -l | tr -d ' ')"
  if [[ "${count}" -ne 1 ]]; then
    echo "ERROR: expected exactly one genomic FASTA for ${accession}; found ${count}."
    find "${package_dir}/ncbi_dataset/data" -type f -name '*.fna' -print
    exit 1
  fi

  # Locate the single FASTA deterministically.  sort + head provides a stable path
  # selection after the preceding count check has already established count=1.
  source_fasta="$(find "${package_dir}/ncbi_dataset/data" -type f -name '*.fna' -print | sort | head -n 1)"
  # Copy, rather than modify, the NCBI-provided FASTA into the project's stable
  # data/raw naming scheme.
  cp "${source_fasta}" "${final_fasta}"

  if [[ ! -s "${final_fasta}" ]]; then
    echo "ERROR: downloaded FASTA is missing or empty: ${final_fasta}"
    exit 1
  fi

  # Append this validated path to the Bash array for the all-genome SeqKit summary.
  FASTA_FILES+=("${final_fasta}")
done < "${CONFIG}"

echo
echo "Generating assembly statistics..."
# `-a` requests extended statistics and `-T` emits tab-delimited output.  Passing
# all FASTAs in one invocation creates one consistent assembly-statistics table.
seqkit stats -a -T "${FASTA_FILES[@]}" > results/genome_stats.tsv

# Hand off FASTA-aware record parsing to Python.  This avoids fragile shell/AWK
# parsing of descriptive FASTA headers and creates the chromosome-specific files
# needed by MUMmer and Circos.
python "${SCRIPT_DIR}/01_parse_genomes.py"

echo
echo "Reference chromosome FASTA:"
# Print a quick human-readable sanity check for the prepared 15-chromosome
# Fol4287 reference FASTA into the stage log.
seqkit stats data/analysis/Fol4287.chromosomes.fna

echo
echo "Fo47 chromosome FASTA:"
# Print the analogous sanity check for Fo47's 12 chromosome records.
seqkit stats data/analysis/Fo47.chromosomes.fna

echo
echo "Downloaded genomes: ${#FASTA_FILES[@]}"
echo "STEP 1 COMPLETE"
