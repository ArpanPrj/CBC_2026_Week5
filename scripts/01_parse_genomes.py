#!/usr/bin/env python3

# =============================================================================
# FASTA PARSING AND CHROMOSOME-LEVEL INPUT PREPARATION
# =============================================================================
# This Python helper is called by Step 1 after all complete NCBI genomic FASTAs
# have been downloaded.  It performs two related tasks:
#
# 1. Build a general metadata summary for every configured assembly by counting
#    FASTA records, total bases, and RefSeq-style NC_/NW_ record prefixes.
# 2. For Fol4287 and Fo47 specifically, identify chromosome records from the
#    actual FASTA headers, write chromosome-only FASTAs, and write metadata that
#    connects RefSeq sequence IDs to human-readable chromosome labels/classes.
#
# The full query genomes other than Fo47 remain as downloaded.  Step 2 aligns
# the chromosome-only Fol4287 reference against each full query assembly.
# =============================================================================
# Postpone evaluation of type annotations.  This keeps annotations lightweight
# and avoids some version-dependent runtime evaluation behavior.
from __future__ import annotations

# csv handles deterministic tab-delimited configuration/result tables.
import csv
# re extracts chromosome labels from descriptive FASTA headers.
import re
# pathlib provides explicit, cross-platform path composition instead of manual
# string concatenation.
from pathlib import Path

# __file__ is this Python script; parent is scripts/, and parent.parent is the
# repository root.  Every later path is derived from that root.
ROOT = Path(__file__).resolve().parent.parent
# Tracked genome-panel configuration.
GENOMES_TSV = ROOT / "config" / "genomes.tsv"
# Tracked Fol4287 chromosome classification table (core vs lineage-specific).
CLASS_TSV = ROOT / "config" / "fol4287_chromosome_classes.tsv"
# Directory populated by 01_download_genomes.sh with unmodified full FASTAs.
RAW_DIR = ROOT / "data" / "raw"
# Directory for chromosome-only derived FASTAs used by specialized downstream
# comparisons/visualizations.
ANALYSIS_DIR = ROOT / "data" / "analysis"
# Directory for compact, tracked scientific output tables.
RESULTS_DIR = ROOT / "results"

ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


# -----------------------------------------------------------------------------
# TABLE READER
# -----------------------------------------------------------------------------
# Return the entire TSV as a list of dictionaries keyed by its header names.
# newline="" is the csv module's recommended way to avoid newline translation.
def read_tsv(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


# -----------------------------------------------------------------------------
# STREAMING FASTA READER
# -----------------------------------------------------------------------------
# Yield one (header, sequence) pair at a time.  Sequence lines are accumulated
# only for the current record, whitespace is removed, and malformed FASTA with
# sequence before any header raises an explicit error.
def fasta_records(path: Path):
    # `None` means no FASTA record has started yet.
    header = None
    # Sequence pieces are buffered as a list and joined only when a record ends.
    seq_chunks = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n\r")
            if not line:
                continue
            # A new `>` header closes/yields the previous record before resetting
            # state for the new sequence.
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:]
                seq_chunks = []
            else:
                if header is None:
                    raise ValueError(f"Sequence encountered before FASTA header in {path}")
                seq_chunks.append("".join(line.split()))
    if header is not None:
        yield header, "".join(seq_chunks)


# -----------------------------------------------------------------------------
# FASTA HEADER SPLITTER
# -----------------------------------------------------------------------------
# NCBI FASTA headers begin with a sequence accession followed by descriptive
# text.  Splitting only once preserves the entire remaining description.
def split_header(header: str):
    parts = header.split(maxsplit=1)
    seq_id = parts[0]
    desc = parts[1] if len(parts) == 2 else ""
    return seq_id, desc


# -----------------------------------------------------------------------------
# CHROMOSOME LABEL EXTRACTION
# -----------------------------------------------------------------------------
# Search the description for text beginning with the word "chromosome" and
# capture characters up to the next comma.  This supports Fol4287 Arabic labels
# (e.g. 1, 14) and Fo47 Roman labels (e.g. VII, XII) without hard-coding IDs.
def chromosome_label(description: str):
    match = re.search(r"\bchromosome\s+([^,]+)", description, flags=re.IGNORECASE)
    return match.group(1).strip() if match else ""


# -----------------------------------------------------------------------------
# DETERMINISTIC FASTA WRITER
# -----------------------------------------------------------------------------
# Preserve the original record header/sequence while wrapping output sequence at
# a fixed 60 bases per line.  Explicit newline="\\n" makes line endings stable.
def write_fasta(records, path: Path, width: int = 60):
    with path.open("w", encoding="utf-8", newline="\n") as out:
        for header, seq in records:
            out.write(f">{header}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i : i + width] + "\n")


# Read the complete configured genome panel once for assembly-level summaries.
genomes = read_tsv(GENOMES_TSV)
# Convert the Fol4287 chromosome-class table into a direct label -> class lookup
# dictionary used when preparing the reference chromosome metadata.
classes = {row["chromosome"]: row["chromosome_class"] for row in read_tsv(CLASS_TSV)}

# =============================================================================
# GENERAL ASSEMBLY SUMMARY FOR EVERY CONFIGURED GENOME
# =============================================================================
# For each raw FASTA, count records and bases.  NC_ is typically used for RefSeq
# chromosome/complete sequence records and NW_ for scaffold-level RefSeq records;
# the counts are descriptive metadata and do not alter downstream sequences.
# General assembly summary.
summary_rows = []
for genome in genomes:
    path = RAW_DIR / f"{genome['name']}.fna"
    if not path.exists():
        raise FileNotFoundError(path)
    nseq = 0
    total = 0
    nc = 0
    nw = 0
    other = 0
    for header, seq in fasta_records(path):
        seq_id, _ = split_header(header)
        nseq += 1
        total += len(seq)
        if seq_id.startswith("NC_"):
            nc += 1
        elif seq_id.startswith("NW_"):
            nw += 1
        else:
            other += 1
    summary_rows.append(
        {
            **genome,
            "num_sequences": nseq,
            "total_bp": total,
            "NC_records": nc,
            "NW_records": nw,
            "other_records": other,
        }
    )

# Write one row per configured genome, retaining the original configuration
# columns alongside the newly calculated sequence statistics.
with (RESULTS_DIR / "genome_metadata.tsv").open("w", newline="", encoding="utf-8") as out:
    fields = [
        "name", "accession", "role", "host_context", "circos",
        "num_sequences", "total_bp", "NC_records", "NW_records", "other_records",
    ]
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows)


# =============================================================================
# PREPARE A CHROMOSOME-ONLY FASTA FOR A NAMED GENOME
# =============================================================================
# Selection criterion:
#   * sequence ID begins with NC_, AND
#   * a chromosome label can be parsed from the FASTA description.
#
# `expected_count` acts as a strong sanity check against unexpected assembly or
# header changes.  `classify` is an optional function that assigns a biological
# class label to each parsed chromosome name.
def prepare_chromosome_genome(name: str, expected_count: int, classify=None):
    raw = RAW_DIR / f"{name}.fna"
    selected = []
    metadata = []
    for header, seq in fasta_records(raw):
        seq_id, desc = split_header(header)
        chrom = chromosome_label(desc)
        if seq_id.startswith("NC_") and chrom:
            selected.append((header, seq))
            metadata.append(
                {
                    "genome": name,
                    "sequence_id": seq_id,
                    "chromosome": chrom,
                    "length_bp": len(seq),
                    "chromosome_class": classify(chrom) if classify else "",
                    "description": desc,
                }
            )
    if len(selected) != expected_count:
        raise RuntimeError(
            f"Expected {expected_count} chromosome records for {name}, found {len(selected)}"
        )
    write_fasta(selected, ANALYSIS_DIR / f"{name}.chromosomes.fna")
    return metadata


# Fol4287 is expected to contribute 15 chromosome records.  The class lookup
# marks the whole lineage-specific chromosomes according to the tracked config
# table and labels the remaining chromosomes as configured.
fol_meta = prepare_chromosome_genome(
    "Fol4287", 15, classify=lambda chrom: classes.get(chrom, "unclassified")
)
# Sort the Fol4287 metadata numerically so chromosome 10 follows chromosome 9
# rather than chromosome 1 in lexical string order.
# Numeric reference order.
fol_meta.sort(key=lambda row: int(row["chromosome"]))

# Fo47 is expected to contribute 12 chromosome records.  Chromosome VII receives
# the accessory label used only for the detailed Fo47 visualization; the other
# Fo47 chromosomes are labeled core.
fo47_meta = prepare_chromosome_genome(
    "Fo47", 12, classify=lambda chrom: "accessory" if chrom.upper() == "VII" else "core"
)
# Explicit Roman-numeral order guarantees deterministic Fo47 chromosome ordering
# in tables and Circos sectors.
roman_order = {r: i for i, r in enumerate(
    ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII"], start=1
)}
fo47_meta.sort(key=lambda row: roman_order.get(row["chromosome"].upper(), 999))

# Write separate metadata files for the reference and detailed Fo47 comparison.
# Downstream scripts use the RefSeq sequence_id field to map MUMmer coordinates
# back to chromosome labels, sizes, and biological classes.
for filename, rows in [
    ("reference_chromosome_metadata.tsv", fol_meta),
    ("Fo47_chromosome_metadata.tsv", fo47_meta),
]:
    with (RESULTS_DIR / filename).open("w", newline="", encoding="utf-8") as out:
        fields = [
            "genome", "sequence_id", "chromosome", "length_bp", "chromosome_class", "description"
        ]
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

print(f"Prepared {len(fol_meta)} Fol4287 chromosome records")
print(f"Prepared {len(fo47_meta)} Fo47 chromosome records")
