# Fol4287 five-genome synteny

A reproducible laptop-scale comparative-genomics workflow that asks whether the four whole lineage-specific chromosomes of *Fusarium oxysporum* f. sp. *lycopersici* Fol4287 are consistently less conserved than its core chromosomes across four additional *F. oxysporum* genomes.

The workflow was deliberately reduced from a larger ten-query panel after benchmarking showed that some genome comparisons were substantially slower. The final panel keeps a biologically useful mix of nonpathogenic/endophytic and pathogenic strains while keeping the project closer to the intended laptop-scale runtime.

## Biological question

**Are Fol4287 lineage-specific chromosomes consistently less conserved than core chromosomes across a small, diverse panel of *Fusarium oxysporum* genomes?**

Fol4287 chromosomes 3, 6, 14, and 15 are treated as whole lineage-specific chromosomes. Chromosomes 1, 2, 4, 5, and 7–13 are treated as core at the whole-chromosome level. This is intentionally a whole-chromosome classification: chromosomes 1 and 2 are known to contain lineage-specific regions even though they are not classified here as wholly lineage-specific.

## Genome panel

| Name | Versioned assembly | Context | Role |
|---|---|---|---|
| Fol4287 | GCF_000149955.1 | tomato wilt | reference |
| Fo47 | GCF_013085055.1 | endophyte/biocontrol | query; detailed Circos comparison |
| FoCL57 | GCA_000260155.3 | tomato crown/root rot | query |
| FoCotton | GCA_000260175.2 | cotton wilt | query |
| FolD11 | GCA_003977725.1 | tomato wilt | query |

The exact accessions used by the workflow are stored in `config/genomes.tsv`.

## Why these four query genomes?

The analysis keeps four query genomes rather than the earlier larger panel. `Fo5176` and `FoHDV247` were removed after local benchmarking because each required more than five minutes for the MUMmer comparison on the development laptop. The retained panel still includes a nonpathogenic/endophytic strain, pathogens associated with different hosts/disease contexts, and a second tomato-wilt isolate.

This is a scientific reduction of the panel rather than an artificial slowdown or speedup: every retained query is analyzed exactly once with the same MUMmer workflow.

## What the pipeline does

1. Bootstraps a pinned Miniforge installation if Conda is absent.
2. Creates the pinned analysis environment from `environment.yml`.
3. Downloads five exact versioned public genome assemblies with NCBI Datasets.
4. Extracts the 15 chromosome records from Fol4287 and the 12 chromosome records from Fo47.
5. Runs four sequential Fol4287-vs-query comparisons with MUMmer4 `dnadiff`.
6. Calculates, for every Fol4287 chromosome and every query genome, unique aligned reference bp, percent aligned, alignment-block count, weighted identity, largest block, and orientation-specific aligned bp.
7. Summarizes core versus lineage-specific conservation for each query genome.
8. Produces a chromosome-by-genome heatmap and a grouped core-vs-lineage-specific summary figure.
9. Retains a detailed Fol4287-vs-Fo47 Circos plot rather than putting all genomes into one crowded circular plot.
10. Writes `CHECKSUMS.txt` with SHA256 hashes for every downloaded/derived data file and every result file.


## Data

All genomic input data are public NCBI genome assemblies and are downloaded automatically by the pipeline. The repository therefore does not store large FASTA files in Git. Instead, the exact versioned assembly accessions are recorded in `config/genomes.tsv`, and `scripts/01_download_genomes.sh` retrieves them with NCBI Datasets every time the analysis is started from a clean clone.

The pipeline creates two generated data directories:

- `data/raw/` — contains the complete genomic FASTA downloaded for each of the five assemblies. These files are the original sequence inputs used by the workflow. The expected files are `Fol4287.fna`, `Fo47.fna`, `FoCL57.fna`, `FoCotton.fna`, and `FolD11.fna`.
- `data/analysis/` — contains chromosome-only FASTA files generated for analyses that require chromosome-scale records. `Fol4287.chromosomes.fna` contains the 15 Fol4287 chromosome records used as the MUMmer reference. `Fo47.chromosomes.fna` contains the 12 Fo47 chromosome records used for the detailed Circos visualization.

The other query genomes are compared to the Fol4287 chromosome-only reference using their complete downloaded assemblies in `data/raw/`. This preserves all query assembly sequence while keeping the Fol4287 reference restricted to its chromosome-scale `NC_` records.

Genome metadata and sequence structure are parsed by `scripts/01_parse_genomes.py`. The script records the number of sequences, total assembly size, and counts of `NC_`, `NW_`, and other sequence records in `results/genome_metadata.tsv`. It also writes:

- `results/reference_chromosome_metadata.tsv` — Fol4287 chromosome identifiers, chromosome numbers, lengths, chromosome class, and FASTA descriptions.
- `results/Fo47_chromosome_metadata.tsv` — Fo47 chromosome identifiers, chromosome labels, lengths, and core/accessory classification.

Fol4287 chromosome classes are not inferred during the run. They are explicitly defined in the tracked file `config/fol4287_chromosome_classes.tsv`, where chromosomes 3, 6, 14, and 15 are designated lineage-specific and the remaining whole chromosomes are designated core.

Because `data/` is generated automatically, it is excluded from Git by `.gitignore`. The data are nevertheless included in `CHECKSUMS.txt`, so the exact FASTA files used in a successful run can be verified with SHA256 hashes.

## Scripts

The analysis is organized as numbered scripts so each computational stage can be inspected and, if necessary, run separately. `run_all.sh` executes them in order after creating the software environment.

### `setup.sh`

Bootstraps the software environment. It detects the operating system and CPU architecture, downloads the pinned Miniforge installer if necessary, verifies the Miniforge installer with its official SHA256 checksum, creates the Conda environment from `environment.yml`, activates the environment, and records the installed software versions in `software_versions.tsv`.

### `scripts/common.sh`

Provides shared setup used by all numbered shell scripts. It identifies the project root, sets deterministic locale/time-zone variables, loads the project-local Miniforge installation, activates the analysis environment, and defines a cross-platform SHA256 helper that works with either `sha256sum` or `shasum -a 256`.

### `scripts/01_download_genomes.sh`

Reads the exact genome accessions from `config/genomes.tsv`, downloads each assembly with NCBI Datasets, extracts the genomic FASTA, and saves it under a stable genome name in `data/raw/`. It then runs SeqKit to generate `results/genome_stats.tsv` and calls `scripts/01_parse_genomes.py` to prepare chromosome-scale files and metadata.

### `scripts/01_parse_genomes.py`

Parses FASTA records and NCBI FASTA descriptions. It summarizes assembly structure, selects the 15 chromosome records from Fol4287 and the 12 chromosome records from Fo47, writes chromosome-only FASTA files to `data/analysis/`, and creates the chromosome metadata tables used by later analyses and figures.

### `scripts/02_run_mummer.sh`

Runs the main external comparative-genomics analysis. Fol4287 chromosome sequences are used as the reference, and each of the four query genomes is analyzed sequentially with MUMmer4 `dnadiff`. The script saves the direct `dnadiff` reports in `results/reports/`, one-to-one coordinate tables in `results/pairwise_coords/`, and an index of pairwise outputs in `results/pairwise_outputs.tsv`. Raw/intermediate MUMmer files remain in `tmp/mummer/`.

### `scripts/03_summarize_multigenome.sh`

Shell wrapper for the multi-genome summary step. It loads the reproducible environment, checks that all expected MUMmer coordinate tables are present, records the run in `logs/`, and executes `scripts/03_summarize_multigenome.py`.

### `scripts/03_summarize_multigenome.py`

Converts the pairwise MUMmer coordinates into biological summaries. For every Fol4287 chromosome and every query genome, it calculates unique aligned Fol4287 bases, percentage aligned, alignment-block count, weighted nucleotide identity, largest alignment block, and alignment orientation. It then generates the multi-genome result tables used by the figures, including `multigenome_chromosome_synteny.tsv`, `percent_aligned_matrix.tsv`, `core_vs_lineage_specific_by_query.tsv`, and `chromosome_across_query_summary.tsv`. It also prepares the filtered Fo47 link table used by the Circos plot.

### `scripts/04_plot_summary.sh`

Shell wrapper for the main multi-genome figures. It activates the environment, writes a plotting log, and runs `scripts/04_plot_summary.R`.

### `scripts/04_plot_summary.R`

Produces the two main summary visualizations: `Fol4287_multigenome_heatmap.svg`, showing the percentage of each Fol4287 chromosome aligned to each query genome, and `core_vs_lineage_specific_by_query.svg`, comparing total conservation of core versus lineage-specific Fol4287 chromosomes across the four query assemblies.

### `scripts/05_plot_circos.sh`

Shell wrapper for the detailed Fol4287-versus-Fo47 chromosome-scale visualization. It exports the plotting thresholds from `config/parameters.sh`, records a log, and runs `scripts/05_plot_circos.R`.

### `scripts/05_plot_circos.R`

Creates `Fol4287_vs_Fo47_circos.svg`. The figure displays the 15 Fol4287 chromosomes and 12 Fo47 chromosomes, distinguishes Fol4287 core and lineage-specific chromosomes, marks Fo47 chromosome VII as accessory, and draws filtered one-to-one MUMmer links between the two assemblies. The Circos thresholds are visualization filters only and do not change the quantitative chromosome-conservation summaries.

### `scripts/06_checksums.sh`

Generates the final reproducibility manifest. It finds every file under `data/` and `results/`, computes a SHA256 hash for each file, writes the hashes to `CHECKSUMS.txt`, and immediately verifies the completed manifest.

### `run_all.sh`

Master workflow entry point. It runs `setup.sh`, removes previously generated `data/`, `tmp/`, and `results/` directories so the analysis starts cleanly, executes Steps 1–6 in order, verifies that all required final outputs exist, and reports the total wall-clock runtime. A reproducer should normally run this file rather than invoking individual scripts manually.

## Requirements

### Windows

Use WSL Ubuntu, not PowerShell. For best filesystem performance, place the project under the Linux home directory rather than `/mnt/c` or `/mnt/d`.

### Linux/macOS

A normal Bash shell, Git, internet access, and either `curl` or `wget` are sufficient. Conda does **not** need to be preinstalled.

## Run from a clean clone

```bash
git clone https://github.com/ArpanPrj/CBC_2026_Week5.git
cd CBC_2026_Week5
bash run_all.sh
```

That is the complete reproduction command. `run_all.sh` invokes `setup.sh`, so a user who sees `conda: command not found` does not need to solve that manually.

To run setup separately:

```bash
bash setup.sh
```

Then run the analysis:

```bash
bash run_all.sh
```

## Expected runtime

The workflow performs four sequential whole-genome MUMmer comparisons rather than adding artificial waits or deliberately inefficient computation. Runtime depends on processor, disk speed, network speed, and whether the Conda environment already exists.

Record the measured runtime from the end of `logs/run_all.log` after the final test.

**Measured runtime for submission:** `FILL_IN_AFTER_FINAL_TEST`

If the final measured runtime is somewhat below 30 minutes, report the actual runtime rather than artificially slowing the workflow. The assignment target is an approximate scale, not a requirement to waste computation.

## Main outputs

The most important scientific outputs are:

- `results/multigenome_chromosome_synteny.tsv` — long-format chromosome × query summary.
- `results/percent_aligned_matrix.tsv` — 15 Fol4287 chromosomes × four query genomes.
- `results/core_vs_lineage_specific_by_query.tsv` — direct core-vs-lineage-specific comparison for every query.
- `results/chromosome_across_query_summary.tsv` — each Fol4287 chromosome summarized across all queries.
- `results/Fol4287_multigenome_heatmap.svg` — main multi-genome visualization.
- `results/core_vs_lineage_specific_by_query.svg` — class-level summary visualization.
- `results/Fol4287_vs_Fo47_circos.svg` — detailed pairwise synteny example.
- `results/reports/` — direct MUMmer `dnadiff` reports.
- `results/pairwise_coords/` — direct one-to-one MUMmer coordinate tables.
- `CHECKSUMS.txt` — SHA256 manifest for all data and result files.

Downloaded FASTA files are intentionally excluded from Git with `.gitignore`; they are re-downloaded from exact versioned accessions.

## Interpreting the figures

### Multi-genome heatmap

Each cell is the percentage of a Fol4287 chromosome covered by one-to-one MUMmer alignments to a query assembly. Rows marked `(LS)` correspond to Fol4287 chromosomes 3, 6, 14, and 15.

This heatmap is the main multi-genome figure because combining all query genomes in a single Circos plot would create severe visual clutter.

### Core versus lineage-specific plot

The grouped bar plot compares the total percentage of Fol4287 core-chromosome sequence and whole lineage-specific-chromosome sequence aligned to each query genome.

### Detailed Fo47 Circos plot

Fol4287 core chromosomes and whole lineage-specific chromosomes are shown separately. Fo47 chromosome VII is marked as an accessory chromosome. Accessory status does not imply that it is a pathogenicity chromosome.

Alignment orientation is labeled `same` or `opposite`; an opposite assembly orientation is not automatically interpreted as a biological inversion.

Only large/high-identity links are displayed in the Circos figure. These thresholds affect visualization only, not the quantitative chromosome-coverage calculations.

## Reproducibility test

After a complete run:

```bash
cp CHECKSUMS.txt /tmp/CHECKSUMS.run1.txt
bash run_all.sh
diff -u /tmp/CHECKSUMS.run1.txt CHECKSUMS.txt
```

No `diff` output means the checksummed files were byte-identical between the two runs on that machine. Across operating systems, text/numeric results are the primary reproducibility targets; SVG bytes can potentially differ because of graphics/font rendering even when the biological result is identical.

## Project layout

```text
fol4287_five_genome_synteny/
├── README.md
├── METHODS.md
├── environment.yml
├── setup.sh
├── run_all.sh
├── config/
│   ├── bootstrap.sh
│   ├── genomes.tsv
│   ├── parameters.sh
│   └── fol4287_chromosome_classes.tsv
├── scripts/
│   ├── common.sh
│   ├── 01_download_genomes.sh
│   ├── 01_parse_genomes.py
│   ├── 02_run_mummer.sh
│   ├── 03_summarize_multigenome.sh
│   ├── 03_summarize_multigenome.py
│   ├── 04_plot_summary.sh
│   ├── 04_plot_summary.R
│   ├── 05_plot_circos.sh
│   ├── 05_plot_circos.R
│   └── 06_checksums.sh
├── data/       # generated; ignored by Git
├── tmp/        # generated; ignored by Git
├── logs/       # generated; ignored by Git
└── results/    # generated and committed after the final run
```

## GitHub publication after validation

Only after the complete pipeline has been tested twice:

```bash
git init
git add .
git commit -m "Reproducible Fol4287 five-genome synteny analysis"
git branch -M main
git remote add origin YOUR_REPOSITORY_URL
git push -u origin main
```

Do not manually add `data/`, `tmp/`, or `logs/`; `.gitignore` excludes them.

## Key references

- Ma LJ et al. 2010. Comparative genomics reveals mobile pathogenicity chromosomes in *Fusarium*. *Nature* 464:367–373.
- Wang B et al. 2020. Chromosome-scale genome assembly of *Fusarium oxysporum* strain Fo47, a fungal endophyte and biocontrol agent. *Molecular Plant-Microbe Interactions*.
- MUMmer4 documentation: `dnadiff` performs whole-genome comparison and provides one-to-one coordinate and summary outputs used by this workflow.
