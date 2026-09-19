#!/usr/bin/env Rscript

# =============================================================================
# STEP 5: DETAILED CHROMOSOME-SCALE CIRCOS PLOT, FOL4287 vs FO47
# =============================================================================
# The multi-genome comparison is summarized by Step 4.  This separate figure
# deliberately keeps ribbon-level detail to one pair so chromosome relationships
# remain visually interpretable.
#
# Sector colors:
#   Fol4287 core chromosomes             blue
#   Fol4287 lineage-specific chromosomes orange
#   Fo47 core chromosomes                green
#   Fo47 accessory chromosome VII        purple
#
# Ribbon colors:
#   same coordinate orientation          blue
#   opposite coordinate orientation      red
#
# “Opposite” describes the MUMmer coordinate orientation and is not automatically
# interpreted as a biological inversion, because assembly sequence orientation
# can itself differ.
#
# Only comments have been added to this annotated copy; executable R statements,
# thresholds, colors, coordinates, ordering, and output filenames are unchanged.
# =============================================================================

# Keep text metadata as character strings and avoid unnecessary scientific
# notation in ordinary genomic coordinates/threshold labels.
options(stringsAsFactors = FALSE, scipen = 999)
# Load circlize while suppressing package startup messages so the pipeline log is
# focused on project messages and errors.
suppressPackageStartupMessages(library(circlize))

# Recover the repository root from Rscript's --file path, making the script
# independent of the shell's current working directory.
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (length(script_arg) == 0) stop("Could not determine script location")
script_file <- sub("^--file=", "", script_arg[1])
root <- normalizePath(file.path(dirname(normalizePath(script_file)), ".."))
setwd(root)

# Read the visualization thresholds exported by the Step 5 shell wrapper.  The
# defaults match config/parameters.sh if this R script is run directly.
plot_min_length <- as.numeric(Sys.getenv("PLOT_MIN_ALIGNMENT_LENGTH", unset = "50000"))
plot_min_identity <- as.numeric(Sys.getenv("PLOT_MIN_IDENTITY", unset = "90"))

# Fol4287 metadata provides RefSeq IDs, chromosome labels, lengths, and the
# core/lineage-specific class assigned in Step 1.
fol <- read.delim("results/reference_chromosome_metadata.tsv", check.names = FALSE)
# Fo47 metadata provides the equivalent information for 12 chromosomes,
# including the core/accessory classification used for sector colors.
fo47 <- read.delim("results/Fo47_chromosome_metadata.tsv", check.names = FALSE)
# The link table is the plotting-only subset constructed in Step 3 after applying
# the configured minimum alignment length and percent identity.
links <- read.delim("results/fo47_circos_links.tsv", check.names = FALSE)

# Strong record-count checks detect a changed/unexpected assembly or malformed
# metadata before any sectors are drawn.
if (nrow(fol) != 15) stop("Expected 15 Fol4287 chromosomes")
if (nrow(fo47) != 12) stop("Expected 12 Fo47 chromosomes")
if (nrow(links) == 0) stop("No Fo47 links passed the plotting filters")

# Fol4287 chromosome labels are numeric text, so conversion followed by sorting
# ensures the expected 1..15 sector order.
fol$chr_num <- as.integer(fol$chromosome)
fol <- fol[order(fol$chr_num), ]
# Fo47 uses Roman chromosome labels.  An explicit order vector avoids lexical
# ordering such as I, II, III, IX, IV, ... .
roman <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII")
fo47$order <- match(fo47$chromosome, roman)
if (any(is.na(fo47$order))) stop("Unexpected Fo47 chromosome label")
fo47 <- fo47[order(fo47$order), ]

# Internal sector IDs include genome prefixes, preventing naming collisions
# between genomes that both have a chromosome I/1-like label.
fol$sector <- paste0("Fol_", fol$chromosome)
fo47$sector <- paste0("Fo47_", fo47$chromosome)
# Compact display labels are stored independently from internal sector IDs.
fol$label <- paste0("Fol", fol$chromosome)
fo47$label <- paste0("Fo47-", fo47$chromosome)

# Place all Fol4287 sectors first and all Fo47 sectors second.  This intentional
# grouping is preserved later through factor levels passed to circlize.
sector_order <- c(fol$sector, fo47$sector)
# Each sector's angular span is proportional to its actual chromosome length.
sector_lengths <- c(as.numeric(fol$length_bp), as.numeric(fo47$length_bp))
names(sector_lengths) <- sector_order
sector_labels <- c(fol$label, fo47$label)
names(sector_labels) <- sector_order

# Fixed literal colors make biological categories stable across runs and ensure
# the legend has an unambiguous mapping.
color_fol_core <- "#377EB8"
color_fol_ls <- "#FF7F00"
color_fo47_core <- "#4DAF4A"
color_fo47_accessory <- "#984EA3"
color_same <- "#377EB8"
color_opposite <- "#E41A1C"

# Select the Fol4287 color from core/LS class and the Fo47 color from
# core/accessory class, yielding one color value per sector.
sector_colors <- c(
  ifelse(fol$chromosome_class == "core", color_fol_core, color_fol_ls),
  ifelse(fo47$chromosome_class == "accessory", color_fo47_accessory, color_fo47_core)
)
names(sector_colors) <- sector_order

# Build RefSeq-ID-to-sector lookup vectors because MUMmer coordinates identify
# sequences by accession rather than the shorter plotted chromosome labels.
fol_map <- setNames(fol$sector, fol$sequence_id)
fo47_map <- setNames(fo47$sector, fo47$sequence_id)
links$fol_sector <- unname(fol_map[links$fol_id])
links$fo47_sector <- unname(fo47_map[links$fo47_id])
if (any(is.na(links$fol_sector)) || any(is.na(links$fo47_sector))) {
  stop("One or more links could not be mapped to chromosome sectors")
}

# Explicitly coerce imported coordinate/length/identity fields to numeric before
# sorting or passing coordinates to circlize.
numeric_columns <- c("fol_start","fol_end","fo47_start","fo47_end","alignment_length_bp","percent_identity")
for (column in numeric_columns) links[[column]] <- as.numeric(links[[column]])
# Establish deterministic ribbon drawing order: longest first, followed by
# stable chromosome/coordinate tie-breakers.  Smaller ribbons drawn later remain
# visible above larger overlapping ribbons.
links <- links[order(-links$alignment_length_bp, links$fol_sector, links$fol_start,
                     links$fo47_sector, links$fo47_start), ]

# Open a square vector graphics device.  Every subsequent drawing call is written
# into this SVG until dev.off() closes the device.
svg("results/Fol4287_vs_Fo47_circos.svg", width = 12, height = 12, pointsize = 10, family = "sans")
par(mar = c(1, 1, 5, 1))
# Clear any prior circlize global state before initializing this plot.
circos.clear()

# Begin with a small 1.5-degree gap after every chromosome sector.
gaps <- rep(1.5, length(sector_order))
# Enlarge the gap after the last Fol4287 chromosome to visually separate the two
# genome groups.
gaps[length(fol$sector)] <- 10
# Enlarge the terminal gap as well, creating a matching group boundary where the
# circle wraps back from Fo47 to Fol4287.
gaps[length(sector_order)] <- 10
# Configure global circlize geometry before sector initialization: starting angle,
# inter-sector gaps, track spacing, cell padding, and overflow-warning behavior.
circos.par(
  start.degree = 90,
  gap.after = gaps,
  track.margin = c(0.005, 0.005),
  cell.padding = c(0, 0, 0, 0),
  points.overflow.warning = FALSE
)

# Give each chromosome its natural genomic coordinate system from 0 to full
# chromosome length.
xlim_matrix <- cbind(rep(0, length(sector_order)), sector_lengths)
# Factor levels explicitly preserve the designed sector order rather than letting
# R/circlize alphabetize names.
circos.initialize(factors = factor(sector_order, levels = sector_order), xlim = xlim_matrix)

# Draw the colored outer chromosome track.  panel.fun executes once for each
# sector and uses CELL_META to access the active sector's identity and limits.
circos.trackPlotRegion(
  ylim = c(0, 1), track.height = 0.085, bg.border = NA,
  panel.fun = function(x, y) {
    sector <- CELL_META$sector.index
    sector_xlim <- CELL_META$xlim
    circos.rect(sector_xlim[1], 0, sector_xlim[2], 1,
                col = sector_colors[sector], border = "white", lwd = 0.5)
    circos.text(CELL_META$xcenter, 0.5, labels = sector_labels[sector],
                facing = "bending.inside", niceFacing = TRUE,
                cex = 0.50, col = "white", font = 2)
  }
)

# Iterate over every retained Fo47 synteny block and draw one ribbon.  The link
# table was already deterministically sorted above.
for (i in seq_len(nrow(links))) {
  # Use orientation-specific translucent ribbon colors.  Transparency helps
  # dense overlaps remain visible rather than producing an opaque center.
  link_color <- if (links$orientation[i] == "opposite") {
    grDevices::adjustcolor(color_opposite, alpha.f = 0.24)
  } else {
    grDevices::adjustcolor(color_same, alpha.f = 0.18)
  }
  # Connect the Fol4287 interval (point1) to the corresponding Fo47 interval
  # (point2).  These low/high coordinates were normalized by Step 3.
  circos.link(
    sector.index1 = links$fol_sector[i],
    point1 = c(links$fol_start[i], links$fol_end[i]),
    sector.index2 = links$fo47_sector[i],
    point2 = c(links$fo47_start[i], links$fo47_end[i]),
    col = link_color,
    border = NA
  )
}

# Add the figure title outside the chromosome circle.
title("Detailed chromosome-scale synteny: Fol4287 vs Fo47", line = 3, cex.main = 1.2)
# Add a subtitle that explicitly reports the plotting-only length and identity
# thresholds used to decide which ribbons are shown.
mtext(
  paste0("Displayed one-to-one MUMmer links: ≥", format(plot_min_length, big.mark = ",", scientific = FALSE),
         " bp and ≥", plot_min_identity, "% identity"),
  side = 3, line = 1.6, cex = 0.72
)

# Document the four chromosome-sector categories and the two ribbon-orientation
# categories in a single legend.
legend(
  "topleft",
  legend = c(
    "Fol4287 core chromosome",
    "Fol4287 lineage-specific chromosome",
    "Fo47 core chromosome",
    "Fo47 accessory chromosome VII",
    "Same-orientation alignment",
    "Opposite-orientation alignment"
  ),
  fill = c(
    color_fol_core, color_fol_ls, color_fo47_core, color_fo47_accessory,
    grDevices::adjustcolor(color_same, alpha.f = 0.50),
    grDevices::adjustcolor(color_opposite, alpha.f = 0.50)
  ),
  border = NA, bty = "n", cex = 0.73, inset = 0.01
)

# Clear circlize's global state after drawing so no sector configuration leaks
# into any later plot in the same R session.
circos.clear()
# Close the SVG graphics device and flush the finished vector file to disk.
dev.off()
cat("Generated results/Fol4287_vs_Fo47_circos.svg with", nrow(links), "links.\n")
