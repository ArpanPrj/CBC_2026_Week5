#!/usr/bin/env bash

# =============================================================================
# BOOTSTRAP CONFIGURATION: WHERE THE REPRODUCIBLE SOFTWARE STACK LIVES
# =============================================================================
# This file contains *configuration only*.  It does not install software by
# itself.  setup.sh sources this file and uses the values below to decide:
#   1. which exact Miniforge release to download;
#   2. where that Miniforge installation should be placed; and
#   3. where the project-specific Conda environment should be created.
#
# Keeping these values in one small file avoids scattering hard-coded paths or
# version numbers across multiple scripts.  The analysis scripts later source
# common.sh, and common.sh sources this file, so every stage resolves the same
# environment locations.
#
# IMPORTANT REPRODUCIBILITY POINT:
# The version is intentionally pinned rather than using a moving "latest"
# installer.  If a different Miniforge release became current later, a fresh
# reproducer would still request the exact release specified here.
# =============================================================================

# Exact Miniforge release used to bootstrap Conda.
# MINIFORGE_VERSION controls the bootstrap distribution itself, not the
# biological analysis.  setup.sh interpolates this literal version into the
# GitHub release URL and into the local installation-directory name.
MINIFORGE_VERSION="26.7.2-0"

# Keep the software installation outside the repository. This avoids slow
# environments when a Windows/WSL user stores the project under /mnt/c or /mnt/d.
# SOFTWARE_ROOT is deliberately outside the Git repository.  The shell syntax
# ${VARIABLE:-default} means: use the environment variable
# FOL4287_FIVEGENOME_SYNTENY_SOFTWARE_DIR if the user explicitly supplied one;
# otherwise use the cache directory under the user's home directory.  This
# gives advanced users an override without requiring any project-file edits.
SOFTWARE_ROOT="${FOL4287_FIVEGENOME_SYNTENY_SOFTWARE_DIR:-${HOME}/.cache/fol4287_five_genome_synteny}"
# MINIFORGE_PREFIX is the location of the bootstrap Conda distribution.  The
# release number is embedded in the directory name so that the path itself
# records which Miniforge version is being used.
MINIFORGE_PREFIX="${SOFTWARE_ROOT}/miniforge-${MINIFORGE_VERSION}"
# ENV_PREFIX is the dedicated analysis environment created from environment.yml.
# Separating it from the base Miniforge installation prevents project packages
# from being installed into the bootstrap/base environment.
ENV_PREFIX="${SOFTWARE_ROOT}/environment"
