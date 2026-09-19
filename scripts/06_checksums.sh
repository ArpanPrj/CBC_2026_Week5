#!/usr/bin/env bash

# =============================================================================
# STEP 6: SHA256 REPRODUCIBILITY MANIFEST
# =============================================================================
# The assignment requires checksums for the data files actually used and for the
# generated outputs.  This stage builds one lexically sorted file list covering
# data/ and results/, computes SHA256 for every listed file, writes CHECKSUMS.txt,
# and immediately verifies that manifest against the current files.
#
# Runtime logs and temporary MUMmer files are intentionally excluded: they are
# machine/run diagnostics or regenerable intermediates rather than the submitted
# scientific data/result set.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

mkdir -p tmp logs
LOG="${ROOT}/logs/06_checksums.log"
exec > >(tee "${LOG}") 2>&1

# The manifest itself is written at repository root for easy submission/use.
CHECKSUM_FILE="${ROOT}/CHECKSUMS.txt"
# A temporary sorted path list separates file discovery/order from hashing.
FILE_LIST="${ROOT}/tmp/checksum_file_list.txt"
rm -f "${CHECKSUM_FILE}"

# Discover every regular file under data/ and results/, combine the two streams,
# and sort in the C locale so manifest row order is deterministic.
{
  find data -type f -print
  find results -type f -print
} | LC_ALL=C sort > "${FILE_LIST}"

if [[ ! -s "${FILE_LIST}" ]]; then
  echo "ERROR: no data/results files found for checksum generation."
  exit 1
fi

# Hash files one at a time using common.sh's Linux/macOS-compatible helper.  The
# emitted records include both SHA256 digest and relative path.
while IFS= read -r file; do
  sha256_files "${file}"
done < "${FILE_LIST}" > "${CHECKSUM_FILE}"

echo "============================================================"
echo "VERIFYING SHA256 MANIFEST"
echo "============================================================"
# Immediately verify the just-created manifest.  A checksum mismatch makes this
# stage fail instead of allowing a corrupted/incompletely written file to pass.
if command -v sha256sum >/dev/null 2>&1; then
  sha256sum --check "${CHECKSUM_FILE}"
elif command -v shasum >/dev/null 2>&1; then
  shasum -a 256 --check "${CHECKSUM_FILE}"
else
  echo "ERROR: no SHA256 verification utility found."
  exit 1
fi

echo "Files checksummed: $(wc -l < "${CHECKSUM_FILE}" | tr -d ' ')"
echo "Wrote CHECKSUMS.txt"
