# Methods

## Study design

The chromosome-level RefSeq assembly of *Fusarium oxysporum* f. sp. *lycopersici* Fol4287 (`GCF_000149955.1`) is used as the reference. Four versioned query assemblies are downloaded from NCBI: Fo47 (`GCF_013085055.1`), FoCL57 (`GCA_000260155.3`), FoCotton (`GCA_000260175.2`), and FolD11 (`GCA_003977725.1`).

The analysis asks whether Fol4287 chromosomes classified as whole lineage-specific chromosomes (3, 6, 14, and 15) show lower sequence conservation across the query panel than Fol4287 whole core chromosomes.

The final four-query panel was selected after runtime benchmarking of a larger candidate panel. Two particularly slow comparisons, Fo5176 and FoHDV247, were excluded to keep the complete analysis appropriate for a laptop-scale reproducibility exercise while retaining biological diversity among the query assemblies.

## Software environment

`setup.sh` downloads a fixed Miniforge release, verifies the official installer SHA256, installs Miniforge outside the repository, and creates the analysis environment from `environment.yml`. The pipeline uses MUMmer4, NCBI Datasets, SeqKit, Python, R, and the R package `circlize`.

No user-level `conda init` is performed and no existing user Conda installation is modified. Windows users are instructed to run the workflow in WSL rather than PowerShell.

## Genome acquisition

Exact versioned assembly accessions are specified in `config/genomes.tsv`. NCBI Datasets downloads the genomic FASTA for each accession. The full downloaded FASTAs are retained as the input data used in the analysis and are included in `CHECKSUMS.txt`, but they are excluded from GitHub because they are public and reproducibly downloadable.

Fol4287 chromosome records are identified from `NC_` RefSeq records whose FASTA descriptions contain chromosome labels. Exactly 15 chromosome records are required. Fo47 is similarly required to contain exactly 12 chromosome records for the detailed circular comparison. Derived chromosome-only FASTAs are written deterministically with 60 bases per sequence line.

The complete assemblies of FoCL57, FoCotton, and FolD11 are used as query FASTAs, allowing chromosome-level and scaffold/contig-level query assemblies to be compared to the chromosome-level Fol4287 reference.

## Whole-genome alignment

Each query assembly is compared sequentially against the 15-chromosome Fol4287 reference FASTA with MUMmer4 `dnadiff`. `dnadiff` generates raw alignment files, one-to-one filtered alignments, many-to-many alignments, coordinate tables, and summary reports. The analysis uses the one-to-one coordinate output (`.1coords`) for quantitative coverage calculations.

The four pairwise comparisons are run sequentially. Runtime is therefore determined by actual comparative-genomics computation rather than artificial delays.

## Chromosome coverage

For each query genome and each Fol4287 chromosome, all one-to-one MUMmer reference intervals are collected. Overlapping reference intervals are merged before aligned base pairs are calculated, preventing double counting.

Percent aligned is calculated as:

`100 × unique aligned Fol4287 reference bp / Fol4287 chromosome length`.

The workflow also records alignment-block count, length-weighted nucleotide identity, largest alignment-block length, and aligned reference bp in same versus opposite orientation.

Opposite orientation is not automatically interpreted as a biological inversion because chromosome/scaffold orientation can differ between genome assemblies.

## Core versus lineage-specific comparison

Fol4287 chromosomes 3, 6, 14, and 15 are classified as whole lineage-specific chromosomes. Other Fol4287 chromosomes are classified as core at the whole-chromosome level. This whole-chromosome classification does not imply that every base of a core-classified chromosome is core; chromosomes 1 and 2 contain lineage-specific regions.

For each query genome, total Fol4287 reference length, unique aligned bp, total percent aligned, median chromosome percent aligned, alignment-block count, and length-weighted identity are summarized separately for core and whole lineage-specific chromosomes.

A second summary calculates the mean, median, minimum, and maximum percent alignment of each Fol4287 chromosome across the four query genomes.

## Visualization

The primary multi-genome visualization is a heatmap of Fol4287 chromosome percent aligned across Fo47, FoCL57, FoCotton, and FolD11. A grouped bar plot summarizes core versus lineage-specific coverage for each query.

A detailed Circos-style plot is produced only for Fol4287 versus Fo47. This avoids the severe clutter that would result from putting several query assemblies into a single circular plot.

The detailed Fo47 plot displays only one-to-one blocks at least 50 kb long and at least 90% identical. These thresholds are visualization filters only and are not used for the quantitative chromosome-coverage calculations.

Fo47 chromosome VII is displayed as an accessory chromosome. Accessory status is not equated with pathogenicity.

## Checksums and reproducibility

`CHECKSUMS.txt` contains SHA256 hashes for every file under `data/` and `results/`. Intermediate MUMmer files under `tmp/` and logs under `logs/` are excluded.

A complete run can be repeated and the two checksum manifests compared with `diff` to test byte-level reproducibility on the same machine. Text and numeric outputs are the primary cross-platform reproducibility targets because SVG rendering can vary slightly with platform graphics and fonts.

## Runtime reporting

The workflow records total wall-clock runtime at the end of `logs/run_all.log`. The measured value from the final submission laptop should be reported in `README.md`. Runtime should not be padded with `sleep` statements or redundant identical analyses.
