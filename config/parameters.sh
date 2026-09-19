#!/usr/bin/env bash

# =============================================================================
# GLOBAL ANALYSIS PARAMETERS
# =============================================================================
# Every value in this file is sourced by one or more pipeline stages.  Central
# configuration makes the workflow auditable: a reviewer can see the reference
# assembly, the detailed Circos comparison, and the visualization thresholds
# without searching through the implementation scripts.
#
# No computation occurs in this file; it only assigns shell variables.
# =============================================================================

# REFERENCE_NAME is the short project label used when constructing filenames.
# It must correspond to the reference row in config/genomes.tsv.
REFERENCE_NAME="Fol4287"
# REFERENCE_ACCESSION records the exact, versioned NCBI assembly accession for
# the reference genome.  The ".1" suffix is important: it pins a particular
# assembly version rather than an unversioned accession that could resolve to a
# later revision.
REFERENCE_ACCESSION="GCF_000149955.1"

# Detailed Circos comparison retained for this query genome only.
# CIRCOS_QUERY selects exactly one of the query genomes for the detailed ribbon
# plot.  All configured query genomes are still included in the multi-genome
# quantitative analysis; this variable only controls the uncluttered pairwise
# Circos visualization.
CIRCOS_QUERY="Fo47"

# Circos is a visualization, so only large/high-identity blocks are shown.
# PLOT_MIN_ALIGNMENT_LENGTH is a *display* threshold measured in base pairs.
# Alignment blocks shorter than this value are omitted from the detailed
# Circos plot so the visual is readable.  The multi-genome chromosome coverage
# summaries are calculated from the full one-to-one MUMmer coordinate output,
# not from this plotting subset.
PLOT_MIN_ALIGNMENT_LENGTH=50000
# PLOT_MIN_IDENTITY is the second display-only filter.  A value of 90 means a
# block must have at least 90 percent nucleotide identity to be drawn in the
# Fo47 Circos panel.
PLOT_MIN_IDENTITY=90

# PREFIX provides a stable short label that other scripts may use when naming
# analysis products.  It is intentionally independent of absolute file paths.
PREFIX="Fol4287"
