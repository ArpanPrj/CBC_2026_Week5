#!/usr/bin/env Rscript

# =============================================================================
# STEP 4: MULTI-GENOME SUMMARY VISUALIZATIONS
# =============================================================================
# This script turns the deterministic Step 3 TSV tables into two vector figures.
#
# Figure 1 is a chromosome-by-query heatmap.  Each cell reports the percentage
# of one Fol4287 chromosome covered by one-to-one MUMmer alignments to a query
# genome.  The text value and the cell shade encode the same percentage.
#
# Figure 2 is a grouped bar plot.  For each query genome it contrasts the
# aggregate percentage of Fol4287 core-chromosome sequence aligned with the
# aggregate percentage of Fol4287 lineage-specific-chromosome sequence aligned.
#
# Only comments have been added to this annotated copy; plotting calculations,
# dimensions, labels, colors, and filenames are unchanged.
# =============================================================================

# Keep imported text as character data rather than automatically creating
# factors, and discourage scientific notation for ordinary percentage values.
options(stringsAsFactors = FALSE, scipen = 999)

# Rscript supplies the path of the running script in a --file= argument.  The
# following block uses that path to recover the repository root, so the script
# does not depend on the user's current working directory.
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (length(script_arg) == 0) stop("Could not determine script location")
script_file <- sub("^--file=", "", script_arg[1])
root <- normalizePath(file.path(dirname(normalizePath(script_file)), ".."))
setwd(root)

# Step 3 produced exactly two tables required here: the chromosome-by-query
# matrix and the per-query core-vs-lineage-specific summary.
matrix_file <- "results/percent_aligned_matrix.tsv"
class_file <- "results/core_vs_lineage_specific_by_query.tsv"
if (!file.exists(matrix_file) || !file.exists(class_file)) stop("Missing Step 3 outputs")

# Read the heatmap matrix without allowing R to alter column names.  Preserving
# query names exactly keeps plot labels consistent with config/genomes.tsv.
mat_df <- read.delim(matrix_file, check.names = FALSE)
# The first two columns describe the Fol4287 chromosome and its class.  Every
# remaining column is therefore a query genome that should become one heatmap
# column.
query_names <- names(mat_df)[-(1:2)]
values <- as.matrix(mat_df[, query_names, drop = FALSE])
# Explicitly coerce the extracted matrix to numeric storage before using it in
# arithmetic, palette indexing, and numeric text formatting.
storage.mode(values) <- "numeric"

# ============================================================
# =============================================================================
# HEATMAP: FOL4287 CHROMOSOME CONSERVATION ACROSS QUERY GENOMES
# =============================================================================
# The code below uses base R primitives rather than a high-level heatmap package.
# That makes the mapping from percentage values to rectangles and text explicit.
# Heatmap
# ============================================================
# Open the SVG device first.  Every plotting command until dev.off() is written
# to this vector file.  The physical width/height are fixed for stable layout.
svg("results/Fol4287_multigenome_heatmap.svg", width = 12, height = 9, pointsize = 10, family = "sans")
par(mar = c(8, 8, 4, 4), xpd = FALSE)

# Cache matrix dimensions because they are used repeatedly for axes, cell loops,
# and the vertical color key.
n_chr <- nrow(values)
n_q <- ncol(values)
# Generate 101 discrete shades so rounded percentage values 0..100 map directly
# to palette entries 1..101.  Pale cells indicate little aligned sequence and
# darker cells indicate greater Fol4287 coverage.
colors <- colorRampPalette(c("#fff7fb", "#7fcdbb", "#2c7fb8", "#253494"))(101)

plot(
  NA,
  xlim = c(0.5, n_q + 1.35),
  ylim = c(0.5, n_chr + 0.5),
  xlab = "",
  ylab = "",
  xaxt = "n",
  yaxt = "n",
  bty = "n",
  main = "Conservation of Fol4287 chromosomes across Fusarium oxysporum genomes"
)

# Base R increases y upward from the bottom.  Reversing the y coordinates places
# chromosome 1 at the top of the heatmap and chromosome 15 at the bottom.
row_y <- n_chr:1
# Nested loops draw each chromosome-by-query cell individually.  This makes it
# possible to place the numeric percentage directly over its colored rectangle.
for (i in seq_len(n_chr)) {
  for (j in seq_len(n_q)) {
    # Clamp to the biologically meaningful plotting range 0..100 before choosing
    # a palette index.  Under the expected inputs the values already lie in this
    # range; the clamp is a defensive plotting safeguard.
    value <- max(0, min(100, values[i, j]))
    # R indexes vectors from 1, hence the +1 after rounding 0..100.
    col <- colors[round(value) + 1]
    # Draw one unit-width/unit-height heatmap cell centered at the current query
    # and chromosome coordinate.
    rect(j - 0.5, row_y[i] - 0.5, j + 0.5, row_y[i] + 0.5, col = col, border = "white")
    # Add one-decimal numeric values so exact comparisons do not depend solely on
    # judging color intensity.  The continuation line chooses contrasting text
    # color for dark versus light cells.
    text(j, row_y[i], labels = sprintf("%.1f", value), cex = 0.62,
         col = if (value >= 58) "white" else "black")
  }
}

# Query labels are rotated perpendicular to the axis (las=2) to fit multiple
# genome names without overlapping horizontally.
axis(1, at = seq_len(n_q), labels = query_names, las = 2, tick = FALSE, cex.axis = 0.85)
# Build chromosome-axis labels and visibly append “(LS)” only to chromosomes
# whose tracked class is lineage_specific.
row_labels <- paste0("Chr", mat_df$chromosome,
                     ifelse(mat_df$chromosome_class == "lineage_specific", " (LS)", ""))
axis(2, at = row_y, labels = row_labels, las = 1, tick = FALSE, cex.axis = 0.85)

# Draw a custom vertical 0–100% color key beside the heatmap using the identical
# palette used for the cells.
# Vertical color key.
legend_x1 <- n_q + 0.72
legend_x2 <- n_q + 0.93
for (k in 0:99) {
  rect(legend_x1, 1 + k * (n_chr - 2) / 100,
       legend_x2, 1 + (k + 1) * (n_chr - 2) / 100,
       col = colors[k + 1], border = NA)
}
text(n_q + 1.03, 1, "0%", adj = 0, cex = 0.7)
text(n_q + 1.03, n_chr - 1, "100%", adj = 0, cex = 0.7)
text(n_q + 0.82, n_chr, "Fol4287 bp\naligned", cex = 0.65)

# The bottom caption states exactly what a cell represents: coverage of the
# Fol4287 reference chromosome by one-to-one MUMmer alignments.
mtext("Cells show the percentage of each Fol4287 chromosome covered by one-to-one MUMmer alignments.",
      side = 1, line = 6.2, cex = 0.78)
box()
# Close and flush the heatmap SVG before beginning the second figure.
dev.off()

# ============================================================
# =============================================================================
# GROUPED BAR PLOT: CORE VS LINEAGE-SPECIFIC CONSERVATION
# =============================================================================
# This section reads the class-level Step 3 table and extracts one core value and
# one lineage-specific value per query genome.
# Core vs lineage-specific grouped bar plot
# ============================================================
# Preserve the table's literal column names exactly as written by Step 3.
class_df <- read.delim(class_file, check.names = FALSE)
# unique() retains first appearance, so bars follow the configured genome order
# rather than alphabetical order.
query_order <- unique(class_df$query)
# Extract the percent_aligned value for the core chromosome compartment of every
# query genome.
core <- sapply(query_order, function(q) {
  as.numeric(class_df$percent_aligned[class_df$query == q & class_df$chromosome_class == "core"])
})
# Extract the corresponding percent_aligned value for the lineage-specific
# chromosome compartment.
ls <- sapply(query_order, function(q) {
  as.numeric(class_df$percent_aligned[class_df$query == q & class_df$chromosome_class == "lineage_specific"])
})
# barplot(..., beside=TRUE) interprets these two rows as the two side-by-side
# series displayed for each query.
plot_matrix <- rbind(Core = core, Lineage_specific = ls)

# Open the second fixed-size vector graphics device.
svg("results/core_vs_lineage_specific_by_query.svg", width = 11, height = 7, pointsize = 10, family = "sans")
par(mar = c(8, 5, 4, 2))
barplot(
  plot_matrix,
  beside = TRUE,
  names.arg = query_order,
  las = 2,
  ylim = c(0, 100),
  ylab = "Fol4287 chromosome sequence aligned (%)",
  main = "Core versus lineage-specific Fol4287 sequence conservation",
  col = c("#377EB8", "#FF7F00"),
  border = NA,
  legend.text = c("Core chromosomes", "Lineage-specific chromosomes"),
  args.legend = list(x = "topright", bty = "n", cex = 0.85)
)
# Add faint horizontal reference lines every 20 percentage points.  They are
# visual guides only and do not transform the data.
abline(h = seq(0, 100, 20), lty = 3, col = "grey80")
# Close and write the completed grouped-bar SVG.
dev.off()

cat("Generated heatmap and core-vs-lineage-specific bar plot.\n")
