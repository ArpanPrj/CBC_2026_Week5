#!/usr/bin/env python3

# =============================================================================
# STEP 3: PARSE MUMMER COORDINATES AND BUILD MULTI-GENOME SYNTENY SUMMARIES
# =============================================================================
# This is the quantitative core of the workflow.  It reads the one-to-one
# `.1coords` output from every Fol4287-vs-query dnadiff comparison and measures
# how much of each Fol4287 chromosome is represented by those alignments.
#
# Important analytical choices
# ----------------------------
# * Coverage is computed on the Fol4287/reference coordinate system.
# * Overlapping alignment intervals are merged before counting aligned bases,
#   preventing the same reference base from being counted twice.
# * Percent aligned = unique aligned Fol4287 bases / Fol4287 chromosome length.
# * Mean nucleotide identity is weighted by alignment-block length.
# * "same" and "opposite" describe coordinate orientation only; the script does
#   not infer that every opposite-orientation block is a biological inversion.
# * Circos thresholds are applied only when making the Fo47 visualization-link
#   table.  They do not determine chromosome coverage in the multi-genome tables.
# =============================================================================
from __future__ import annotations

# csv reads/writes stable tab-delimited data products.
import csv
# math.nan is used to represent a weighted identity when no alignment blocks exist.
import math
# os.environ supplies plotting thresholds exported by the shell wrapper.
import os
# defaultdict(list) groups parsed alignment blocks efficiently by reference ID.
from collections import defaultdict
from pathlib import Path

# Resolve repository-root-relative paths from this script's own location.
ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
REF_META = RESULTS / "reference_chromosome_metadata.tsv"
FO47_META = RESULTS / "Fo47_chromosome_metadata.tsv"
GENOMES = ROOT / "config" / "genomes.tsv"

# Visualization-only minimum block length.  Defaults mirror parameters.sh in case
# the Python script is ever invoked directly without the wrapper.
PLOT_MIN_LENGTH = int(os.environ.get("PLOT_MIN_ALIGNMENT_LENGTH", "50000"))
# Visualization-only nucleotide identity threshold.
PLOT_MIN_IDENTITY = float(os.environ.get("PLOT_MIN_IDENTITY", "90"))
# Query selected for the detailed pairwise Circos plot; currently Fo47.
CIRCOS_QUERY = os.environ.get("CIRCOS_QUERY", "Fo47")


# -----------------------------------------------------------------------------
# GENERIC TSV READER
# -----------------------------------------------------------------------------
# DictReader preserves named columns and returns strings exactly as represented
# in the file until the analysis explicitly converts numeric fields.
def read_tsv(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


# -----------------------------------------------------------------------------
# GENERIC TSV WRITER
# -----------------------------------------------------------------------------
# The caller supplies an explicit field order.  lineterminator="\\n" prevents
# platform-specific line endings from altering result bytes.
def write_tsv(path: Path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


# -----------------------------------------------------------------------------
# MUMMER .1COORDS PARSER
# -----------------------------------------------------------------------------
# dnadiff emits a tabular one-to-one alignment coordinate file.  For every row
# this parser extracts reference/query start/end positions, aligned lengths,
# percent identity, sequence lengths, coverage values, and terminal sequence
# identifiers.  Numeric conversion is validated immediately so malformed output
# fails at the exact line where parsing became impossible.
def parse_coords(path: Path):
    """Parse dnadiff .1coords output (show-coords -THrcl)."""
    rows = []
    with path.open(encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n\r")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 13:
                raise ValueError(f"Unexpected .1coords format at {path}:{line_no}: {line}")
            try:
                s1, e1, s2, e2 = map(int, parts[0:4])
                len1, len2 = map(int, parts[4:6])
                identity = float(parts[6])
                len_r, len_q = map(int, parts[7:9])
                cov_r, cov_q = map(float, parts[9:11])
            except ValueError as exc:
                raise ValueError(f"Could not parse numeric fields at {path}:{line_no}") from exc
            # Sequence IDs are taken from the final two columns rather than fixed
            # middle positions, matching the `show-coords -THrcl` style layout.
            ref_id, query_id = parts[-2], parts[-1]
            rows.append(
                {
                    "ref_start": s1,
                    "ref_end": e1,
                    "query_start": s2,
                    "query_end": e2,
                    "ref_aln_len": len1,
                    "query_aln_len": len2,
                    "pct_identity": identity,
                    "ref_seq_len": len_r,
                    "query_seq_len": len_q,
                    "ref_cov_pct": cov_r,
                    "query_cov_pct": cov_q,
                    "ref_id": ref_id,
                    "query_id": query_id,
                    # Query coordinates increasing from start to end are labeled "same";
                    # decreasing coordinates are labeled "opposite" relative to the
                    # reference interval orientation.
                    "orientation": "same" if s2 <= e2 else "opposite",
                }
            )
    return rows


# -----------------------------------------------------------------------------
# UNIQUE COVERAGE OF A SET OF INCLUSIVE GENOMIC INTERVALS
# -----------------------------------------------------------------------------
# MUMmer blocks can overlap on the reference.  Summing their lengths directly
# could therefore overestimate chromosome coverage.  This function normalizes
# each interval to low<=high, sorts intervals, merges overlapping/adjacent spans,
# and returns the exact number of reference bases covered at least once.
def union_length(intervals):
    intervals = [(min(a, b), max(a, b)) for a, b in intervals]
    if not intervals:
        return 0
    intervals.sort()
    total = 0
    cur_start, cur_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= cur_end + 1:
            cur_end = max(cur_end, end)
        else:
            total += cur_end - cur_start + 1
            cur_start, cur_end = start, end
    total += cur_end - cur_start + 1
    return total


# -----------------------------------------------------------------------------
# ALIGNMENT-LENGTH-WEIGHTED MEAN
# -----------------------------------------------------------------------------
# A 1 Mb alignment should contribute more to the chromosome-wide identity than
# a 1 kb alignment.  The weighted mean implements that proportional contribution.
# When total weight is zero, NaN signals that no meaningful identity exists.
def weighted_mean(values, weights):
    denom = sum(weights)
    return sum(v * w for v, w in zip(values, weights)) / denom if denom else math.nan


# =============================================================================
# LOAD CONFIGURATION AND REFERENCE METADATA
# =============================================================================
genomes = read_tsv(GENOMES)
# Preserve the query order from genomes.tsv; that same order propagates into the
# matrix columns, grouped summaries, and figures.
queries = [row for row in genomes if row["role"] == "query"]
ref_meta = read_tsv(REF_META)
# Sort Fol4287 chromosomes numerically for deterministic rows 1 through 15.
ref_meta.sort(key=lambda row: int(row["chromosome"]))
# Direct RefSeq ID -> Fol4287 metadata lookup used repeatedly while classifying
# alignment blocks.
ref_by_id = {row["sequence_id"]: row for row in ref_meta}
fo47_meta = read_tsv(FO47_META)
# Direct Fo47 RefSeq ID -> chromosome metadata lookup used for Circos labels.
fo47_by_id = {row["sequence_id"]: row for row in fo47_meta}

if len(ref_meta) != 15:
    raise RuntimeError(f"Expected 15 reference chromosomes, found {len(ref_meta)}")

# =============================================================================
# PARSE EVERY PAIRWISE ONE-TO-ONE COORDINATE FILE
# =============================================================================
# all_blocks becomes a dictionary keyed by query short name.  Each value is the
# list of parsed MUMmer blocks for Fol4287 versus that query.
all_blocks = {}
for query in queries:
    coords = RESULTS / "pairwise_coords" / f"Fol4287_vs_{query['name']}.1coords.tsv"
    if not coords.exists():
        raise FileNotFoundError(coords)
    all_blocks[query["name"]] = parse_coords(coords)

# =============================================================================
# LONG-FORM CHROMOSOME SYNTENY TABLE
# =============================================================================
# One output row is generated for every combination of query genome and one of
# the 15 Fol4287 chromosomes.
# Per-query, per-reference-chromosome coverage.
long_rows = []
for query in queries:
    blocks = all_blocks[query["name"]]
    # Group this query's blocks by Fol4287 sequence ID so every reference
    # chromosome can be summarized independently.
    by_ref = defaultdict(list)
    for block in blocks:
        by_ref[block["ref_id"]].append(block)

    for meta in ref_meta:
        chr_blocks = by_ref.get(meta["sequence_id"], [])
        chr_len = int(meta["length_bp"])
        # Unique coverage across all one-to-one blocks on this Fol4287 chromosome.
        aligned_bp = union_length([(b["ref_start"], b["ref_end"]) for b in chr_blocks])
        # Unique reference bases covered specifically by same-orientation blocks.
        same_bp = union_length(
            [(b["ref_start"], b["ref_end"]) for b in chr_blocks if b["orientation"] == "same"]
        )
        # Unique reference bases covered specifically by opposite-orientation blocks.
        opposite_bp = union_length(
            [(b["ref_start"], b["ref_end"]) for b in chr_blocks if b["orientation"] == "opposite"]
        )
        # Identity weighted by reference aligned length; this is descriptive of the
        # aligned sequence and does not treat each block as an equal replicate.
        identity = weighted_mean(
            [b["pct_identity"] for b in chr_blocks],
            [b["ref_aln_len"] for b in chr_blocks],
        )
        # Largest individual one-to-one block, or zero if this chromosome has no blocks.
        largest = max((b["ref_aln_len"] for b in chr_blocks), default=0)
        long_rows.append(
            {
                "query": query["name"],
                "accession": query["accession"],
                "host_context": query["host_context"],
                "chromosome": meta["chromosome"],
                "ref_id": meta["sequence_id"],
                "chromosome_class": meta["chromosome_class"],
                "chromosome_length_bp": chr_len,
                "aligned_bp": aligned_bp,
                "percent_aligned": f"{100 * aligned_bp / chr_len:.4f}",
                "alignment_blocks": len(chr_blocks),
                "weighted_mean_identity": "NA" if math.isnan(identity) else f"{identity:.4f}",
                "largest_block_bp": largest,
                "same_orientation_aligned_bp": same_bp,
                "opposite_orientation_aligned_bp": opposite_bp,
            }
        )

long_fields = [
    "query", "accession", "host_context", "chromosome", "ref_id", "chromosome_class",
    "chromosome_length_bp", "aligned_bp", "percent_aligned", "alignment_blocks",
    "weighted_mean_identity", "largest_block_bp", "same_orientation_aligned_bp",
    "opposite_orientation_aligned_bp",
]
# Save the detailed long table before deriving matrix/class-level summaries.
write_tsv(RESULTS / "multigenome_chromosome_synteny.tsv", long_rows, long_fields)

# =============================================================================
# CHROMOSOME × QUERY MATRIX FOR THE HEATMAP
# =============================================================================
# Reformat the long table so rows are Fol4287 chromosomes and columns are query
# genomes.  Each cell is the percent of that Fol4287 chromosome covered by the
# query's one-to-one MUMmer alignments.
# Percent-aligned matrix, preserving query order and chromosome order.
lookup = {(r["query"], r["chromosome"]): r for r in long_rows}
matrix_rows = []
for meta in ref_meta:
    row = {
        "chromosome": meta["chromosome"],
        "chromosome_class": meta["chromosome_class"],
    }
    for query in queries:
        row[query["name"]] = lookup[(query["name"], meta["chromosome"])]["percent_aligned"]
    matrix_rows.append(row)
write_tsv(
    RESULTS / "percent_aligned_matrix.tsv",
    matrix_rows,
    ["chromosome", "chromosome_class"] + [q["name"] for q in queries],
)

# =============================================================================
# CORE VS LINEAGE-SPECIFIC COMPARTMENT SUMMARY
# =============================================================================
# For each query, aggregate Fol4287 chromosomes into the two biological classes
# tracked in reference_chromosome_metadata.tsv.  The total percent aligned uses
# summed unique aligned bases divided by summed chromosome lengths.
# Core-vs-lineage-specific summary for each query.
class_rows = []
for query in queries:
    qrows = [r for r in long_rows if r["query"] == query["name"]]
    for class_name in ["core", "lineage_specific"]:
        subset = [r for r in qrows if r["chromosome_class"] == class_name]
        total_len = sum(int(r["chromosome_length_bp"]) for r in subset)
        total_aligned = sum(int(r["aligned_bp"]) for r in subset)
        percents = [float(r["percent_aligned"]) for r in subset]
        class_blocks = [
            b for b in all_blocks[query["name"]]
            if ref_by_id.get(b["ref_id"], {}).get("chromosome_class") == class_name
        ]
        identity = weighted_mean(
            [b["pct_identity"] for b in class_blocks],
            [b["ref_aln_len"] for b in class_blocks],
        )
        # Calculate the median explicitly so behavior is transparent and does not
        # depend on a statistics-library implementation.
        ordered = sorted(percents)
        n = len(ordered)
        median = ordered[n // 2] if n % 2 else (ordered[n // 2 - 1] + ordered[n // 2]) / 2
        class_rows.append(
            {
                "query": query["name"],
                "accession": query["accession"],
                "host_context": query["host_context"],
                "chromosome_class": class_name,
                "chromosome_count": len(subset),
                "total_length_bp": total_len,
                "aligned_bp": total_aligned,
                "percent_aligned": f"{100 * total_aligned / total_len:.4f}",
                "median_chromosome_percent_aligned": f"{median:.4f}",
                "alignment_blocks": len(class_blocks),
                "weighted_mean_identity": "NA" if math.isnan(identity) else f"{identity:.4f}",
            }
        )

class_fields = [
    "query", "accession", "host_context", "chromosome_class", "chromosome_count",
    "total_length_bp", "aligned_bp", "percent_aligned", "median_chromosome_percent_aligned",
    "alignment_blocks", "weighted_mean_identity",
]
write_tsv(RESULTS / "core_vs_lineage_specific_by_query.tsv", class_rows, class_fields)

# =============================================================================
# ACROSS-QUERY SUMMARY FOR EACH FOL4287 CHROMOSOME
# =============================================================================
# Collapse the query dimension to mean/median/min/max conservation for each
# reference chromosome.  This is useful for identifying chromosomes that are
# consistently conserved or consistently poorly represented across the panel.
# Across-query summary by reference chromosome.
across_rows = []
for meta in ref_meta:
    values = [float(lookup[(q["name"], meta["chromosome"])]["percent_aligned"]) for q in queries]
    ordered = sorted(values)
    n = len(ordered)
    median = ordered[n // 2] if n % 2 else (ordered[n // 2 - 1] + ordered[n // 2]) / 2
    across_rows.append(
        {
            "chromosome": meta["chromosome"],
            "chromosome_class": meta["chromosome_class"],
            "mean_percent_aligned_across_queries": f"{sum(values) / len(values):.4f}",
            "median_percent_aligned_across_queries": f"{median:.4f}",
            "min_percent_aligned": f"{min(values):.4f}",
            "max_percent_aligned": f"{max(values):.4f}",
        }
    )
write_tsv(
    RESULTS / "chromosome_across_query_summary.tsv",
    across_rows,
    ["chromosome", "chromosome_class", "mean_percent_aligned_across_queries",
     "median_percent_aligned_across_queries", "min_percent_aligned", "max_percent_aligned"],
)

# =============================================================================
# FILTERED LINK TABLE FOR THE DETAILED FOL4287-vs-FO47 CIRCOS PLOT
# =============================================================================
# This is deliberately separated from the quantitative coverage analysis above.
# Only the visualization applies PLOT_MIN_LENGTH and PLOT_MIN_IDENTITY, reducing
# ribbon clutter without changing any multi-genome percentage calculations.
# Circos links only for the designated detailed comparison (Fo47).
if CIRCOS_QUERY not in all_blocks:
    raise RuntimeError(f"Circos query {CIRCOS_QUERY} is not among the configured query genomes")

circos_rows = []
for block in all_blocks[CIRCOS_QUERY]:
    # Skip blocks failing either display threshold.
    if block["ref_aln_len"] < PLOT_MIN_LENGTH or block["pct_identity"] < PLOT_MIN_IDENTITY:
        continue
    # Map the MUMmer sequence IDs back to biological chromosome metadata.  A block
    # lacking either mapping is skipped rather than plotted with an ambiguous label.
    ref = ref_by_id.get(block["ref_id"])
    qry = fo47_by_id.get(block["query_id"])
    if ref is None or qry is None:
        continue
    circos_rows.append(
        {
            "fol_id": block["ref_id"],
            "fol_chromosome": ref["chromosome"],
            "fol_chromosome_class": ref["chromosome_class"],
            "fol_start": min(block["ref_start"], block["ref_end"]),
            "fol_end": max(block["ref_start"], block["ref_end"]),
            "fo47_id": block["query_id"],
            "fo47_chromosome": qry["chromosome"],
            "fo47_chromosome_class": qry["chromosome_class"],
            "fo47_start": min(block["query_start"], block["query_end"]),
            "fo47_end": max(block["query_start"], block["query_end"]),
            "alignment_length_bp": block["ref_aln_len"],
            "percent_identity": f"{block['pct_identity']:.4f}",
            "orientation": block["orientation"],
        }
    )

# Sort links deterministically by Fol4287 chromosome/coordinate and then Fo47 ID
# and coordinate.  Stable row order improves reproducibility of both the TSV and
# subsequent SVG drawing order.
circos_rows.sort(key=lambda r: (int(r["fol_chromosome"]), int(r["fol_start"]), r["fo47_id"], int(r["fo47_start"])))
circos_fields = [
    "fol_id", "fol_chromosome", "fol_chromosome_class", "fol_start", "fol_end",
    "fo47_id", "fo47_chromosome", "fo47_chromosome_class", "fo47_start", "fo47_end",
    "alignment_length_bp", "percent_identity", "orientation",
]
write_tsv(RESULTS / "fo47_circos_links.tsv", circos_rows, circos_fields)

# Record exactly how many one-to-one blocks existed and how many survived the
# visualization thresholds, including same/opposite orientation counts.
with (RESULTS / "fo47_circos_filter_summary.tsv").open("w", newline="", encoding="utf-8") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(["metric", "value"])
    writer.writerow(["all_one_to_one_blocks", len(all_blocks[CIRCOS_QUERY])])
    writer.writerow(["plot_min_alignment_length_bp", PLOT_MIN_LENGTH])
    writer.writerow(["plot_min_identity_percent", f"{PLOT_MIN_IDENTITY:g}"])
    writer.writerow(["blocks_retained_for_plot", len(circos_rows)])
    writer.writerow(["same_orientation_blocks_for_plot", sum(r["orientation"] == "same" for r in circos_rows)])
    writer.writerow(["opposite_orientation_blocks_for_plot", sum(r["orientation"] == "opposite" for r in circos_rows)])

print(f"Summarized {len(queries)} query genomes across {len(ref_meta)} Fol4287 chromosomes")
print(f"Circos links retained for {CIRCOS_QUERY}: {len(circos_rows)}")
