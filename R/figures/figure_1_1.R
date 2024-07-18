# Loading packages and functions
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)

rstudioapi::getSourceEditorContext()$path %>%
    str_split("/") %>%
    unlist() %>%
    head(-3) %>%
    str_c(collapse = "/") %>%
    str_c("/") %>%
    setwd()

source("R/custom_fct.R")
# Loading data (path to change later)
rawcounts <- readcounts("/home/jules/Documents/phd/Data/Article_veranika/bulk/counts.csv", sep = ",", header = TRUE)
meta <- read.table("/home/jules/Documents/phd/Data/Article_veranika/bulk/metadata.csv", sep = ",", header = TRUE)

# Keeping only necessary samples
meta <- filter(meta, type %in% c("dorsal", "ventral"), CRISPR == "control")

counts <- rawcounts[, meta$sample]
rownames(counts) <- gene_converter(rownames(counts), "ENSEMBL", "SYMBOL")
counts <- counts[which(!is.na(rownames(counts))), ]


markers <- read.table("/home/jules/Documents/phd/Data/Article_veranika/GeneListFig1_18_07_24.csv", sep = ",", header = TRUE)
markers
markers <- markers$gene
setdiff(markers, rownames(counts))
markers <- intersect(markers, rownames(counts))
# Get rows corresponding to markers
retained_row <- counts[rownames(counts) %in% markers, ]

# filtering out lowly expressed genes
countsf <- counts[rowSums(counts) > 50, ]

# putting back potential fitltered out markers
if (length(which(!rownames(retained_row) %in% rownames(countsf))) == 1) {
    countsf <- rbind(countsf, retained_row[which(!rownames(retained_row) %in% rownames(countsf)), ])
    rownames(countsf)[nrow(countsf)] <- rownames(retained_row)[which(!rownames(retained_row) %in% rownames(countsf))]
} else if (length(which(!rownames(retained_row) %in% rownames(countsf))) > 1) {
    countsf <- rbind(countsf, retained_row[which(!rownames(retained_row) %in% rownames(countsf)), ])
}

# With dorsal / ventral distinction
norm <- countsf %>% DESeq2::varianceStabilizingTransformation()
norm <- limma::removeBatchEffect(norm, meta$line)
# scaled_norm <- t(apply(norm, 1, scale))
# colnames(scaled_norm) <- colnames(norm)
# rownames(scaled_norm) <- rownames(norm)
markers %>% length()
sample_ha <- columnAnnotation(
    line = meta[order(meta$type), ]$line,
    type = meta[order(meta$type), ]$type,
    col = list(
        line = c("LON" = "#00f7ff", "WTC" = "#ff0000"),
        type = c("dorsal" = "#A1A1DE", "ventral" = "#80AD3C")
    )
)
png(filename = "results/images/Figure_1/F1_1_marker_HM.png", width = 2000, height = 1800, res = 250)
Heatmap(
    norm[markers, meta[order(meta$type), ]$sample],
    name = "Normalized expression",
    row_title_gp = gpar(fontsize = 16, fontface = "bold"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_side = "left",
    show_column_names = TRUE,
    bottom_annotation = sample_ha,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    width = ncol(norm[markers, ]) * unit(4, "mm"),
    height = nrow(norm[markers, ]) * unit(4, "mm"),
    col = colorRampPalette(c(
        "black",
        "purple",
        "orange",
        "yellow"
    ))(1000),
)
dev.off()
