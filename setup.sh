#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${ROOT}"
source "${ROOT}/config/bootstrap.sh"

mkdir -p logs tmp/bootstrap
LOG="${ROOT}/logs/setup.log"
exec > >(tee "${LOG}") 2>&1

SYSTEM="$(uname -s)"
ARCH="$(uname -m)"

case "${SYSTEM}" in
  Linux)  MF_OS="Linux" ;;
  Darwin) MF_OS="MacOSX" ;;
  *)
    echo "ERROR: Unsupported OS: ${SYSTEM}"
    echo "Windows users should run this project inside WSL."
    exit 1
    ;;
esac

case "${ARCH}" in
  x86_64|amd64) MF_ARCH="x86_64" ;;
  arm64)        MF_ARCH="arm64" ;;
  aarch64)
    if [[ "${MF_OS}" == "MacOSX" ]]; then MF_ARCH="arm64"; else MF_ARCH="aarch64"; fi
    ;;
  *)
    echo "ERROR: Unsupported CPU architecture: ${ARCH}"
    exit 1
    ;;
esac

INSTALLER="Miniforge3-${MINIFORGE_VERSION}-${MF_OS}-${MF_ARCH}.sh"
BASE_URL="https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}"
INSTALLER_URL="${BASE_URL}/${INSTALLER}"
CHECKSUM_URL="${INSTALLER_URL}.sha256"
INSTALLER_PATH="${ROOT}/tmp/bootstrap/${INSTALLER}"
CHECKSUM_PATH="${INSTALLER_PATH}.sha256"

# Download with whichever common command is already present on the host.
download_file() {
  local url="$1"
  local output="$2"
  if command -v curl >/dev/null 2>&1; then
    curl --fail --location --retry 3 --output "${output}" "${url}"
  elif command -v wget >/dev/null 2>&1; then
    wget --tries=3 --output-document="${output}" "${url}"
  else
    echo "ERROR: neither curl nor wget is installed."
    exit 1
  fi
}

echo "============================================================"
echo "REPRODUCIBLE SOFTWARE SETUP"
echo "============================================================"
echo "Project: ${ROOT}"
echo "Platform: ${MF_OS}-${MF_ARCH}"
echo "Miniforge: ${MINIFORGE_VERSION}"

if [[ ! -x "${MINIFORGE_PREFIX}/bin/conda" ]]; then
  echo "Installing pinned Miniforge into ${MINIFORGE_PREFIX}"
  rm -f "${INSTALLER_PATH}" "${CHECKSUM_PATH}"
  download_file "${INSTALLER_URL}" "${INSTALLER_PATH}"
  download_file "${CHECKSUM_URL}" "${CHECKSUM_PATH}"

  cd "${ROOT}/tmp/bootstrap"
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum --check "$(basename "${CHECKSUM_PATH}")"
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 --check "$(basename "${CHECKSUM_PATH}")"
  else
    echo "ERROR: no SHA256 utility found."
    exit 1
  fi
  cd "${ROOT}"

  bash "${INSTALLER_PATH}" -b -p "${MINIFORGE_PREFIX}"
else
  echo "Pinned Miniforge already present; skipping installation."
fi

source "${MINIFORGE_PREFIX}/etc/profile.d/conda.sh"

echo "Conda: $(conda --version)"

if command -v sha256sum >/dev/null 2>&1; then
  ENV_HASH="$(sha256sum environment.yml | awk '{print $1}')"
else
  ENV_HASH="$(shasum -a 256 environment.yml | awk '{print $1}')"
fi
HASH_FILE="${ENV_PREFIX}/.environment_yml_sha256"

if [[ -d "${ENV_PREFIX}" ]]; then
  if [[ ! -f "${HASH_FILE}" ]] || [[ "$(cat "${HASH_FILE}")" != "${ENV_HASH}" ]]; then
    echo "environment.yml changed; rebuilding the analysis environment."
    rm -rf "${ENV_PREFIX}"
  fi
fi

if [[ ! -d "${ENV_PREFIX}" ]]; then
  conda env create --prefix "${ENV_PREFIX}" --file environment.yml --solver libmamba
  echo "${ENV_HASH}" > "${HASH_FILE}"
else
  echo "Analysis environment already matches environment.yml."
fi

conda activate "${ENV_PREFIX}"

{
  echo -e "component\tversion"
  echo -e "conda\t$(conda --version | awk '{print $2}')"
  echo -e "python\t$(python --version 2>&1 | awk '{print $2}')"
  echo -e "mummer\t$(nucmer --version 2>&1 | tail -n 1 | tr -d '\r')"
  echo -e "datasets\t$(datasets version 2>&1 | tail -n 1 | tr -d '\r')"
  echo -e "seqkit\t$(seqkit version 2>&1 | awk '{print $3}')"
  echo -e "R\t$(Rscript --version 2>&1 | sed 's/Rscript (R) version //')"
  echo -e "circlize\t$(Rscript -e 'cat(as.character(packageVersion("circlize")))')"
} > software_versions.tsv

cat software_versions.tsv

echo "============================================================"
echo "SETUP COMPLETE"
echo "============================================================"
