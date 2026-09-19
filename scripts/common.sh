#!/usr/bin/env bash

# =============================================================================
# COMMON RUNTIME INITIALIZATION FOR ALL PIPELINE SHELL SCRIPTS
# =============================================================================
# Every numbered shell stage sources this file first.  Its job is to establish
# exactly the same working directory, locale, timezone, and Conda environment
# regardless of where the user launched a stage from.
#
# This avoids a common reproducibility problem in which a script works only
# when run from a particular directory or only after the user has manually
# activated a shell environment.
# =============================================================================
# Strict Bash mode:
#   -e  aborts when an unhandled command returns a non-zero status;
#   -u  treats an unset variable as an error;
#   -o pipefail propagates failures from any command in a pipeline.
# Together these settings make failures explicit instead of allowing later
# stages to continue with missing or partial inputs.
set -euo pipefail

# BASH_SOURCE[0] is the path of common.sh itself.  Moving one directory upward
# from scripts/ gives the repository root.  cd + pwd canonicalizes the path so
# later relative paths are anchored to the project rather than the caller's
# current directory.
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Exporting ROOT makes it available not only to this shell but also to child
# processes launched by the numbered scripts.
export ROOT
# All subsequent relative paths (config/, data/, results/, etc.) are therefore
# interpreted from the repository root.
cd "${ROOT}"

# The C locale fixes text sorting and numeric/text formatting behavior across
# systems.  This is especially helpful for deterministic tables and checksums.
export LC_ALL=C
# LANG is also fixed to C for programs that consult LANG rather than LC_ALL.
export LANG=C
# UTC avoids machine-specific local timezone effects in programs that may emit
# timestamps or interpret dates.
export TZ=UTC

# Load the pinned Miniforge version and the two installation paths defined in
# config/bootstrap.sh.
source "${ROOT}/config/bootstrap.sh"

# Refuse to proceed if setup.sh has not created both the bootstrap Conda binary
# and the dedicated analysis environment.  The numbered stages intentionally do
# not attempt ad-hoc installation themselves; setup.sh is the single setup path.
if [[ ! -x "${MINIFORGE_PREFIX}/bin/conda" ]] || [[ ! -d "${ENV_PREFIX}" ]]; then
  echo "ERROR: project environment not found. Run: bash setup.sh"
  exit 1
fi

# Sourcing conda.sh defines the shell function required for `conda activate` in
# a non-interactive script.  This does not require `conda init` or modification
# of the user's ~/.bashrc.
source "${MINIFORGE_PREFIX}/etc/profile.d/conda.sh"
# Activate the exact project environment before calling datasets, SeqKit,
# MUMmer, Python, or R.
conda activate "${ENV_PREFIX}"

# Cross-platform SHA256 helper.  GNU/Linux commonly supplies `sha256sum`, while
# macOS commonly supplies `shasum -a 256`.  The function accepts one or more
# path arguments and prints standard checksum records for whichever utility is
# available.
sha256_files() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$@"
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$@"
  else
    echo "ERROR: no SHA256 utility found."
    exit 1
  fi
}
