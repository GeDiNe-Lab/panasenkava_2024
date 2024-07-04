# Loading packages and functions
library(ggplot2)
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)

library(tibble)

rstudioapi::getSourceEditorContext()$path %>%
    str_split("/") %>%
    unlist() %>%
    head(-3) %>%
    str_c(collapse = "/") %>%
    str_c("/") %>%
    setwd()

source("R/custom_fct.R")
# Loading data (path to change later)
counts <- readcounts("/home/jules/Documents/phd/Data/Article_veranika/bulk/counts.csv", sep = ",", header = TRUE)
meta <- read.table("/home/jules/Documents/phd/Data/Article_veranika/bulk/metadata.csv", sep = ",", header = TRUE)

# Keeping only necessary samples
meta <- filter(meta, type %in% c("dorsal", "ventral"), CRISPR == "control")

counts <- counts[, meta$sample]
rownames(counts) <- gene_converter(rownames(counts), "ENSEMBL", "SYMBOL")
counts <- counts[which(!is.na(rownames(counts))), ]

full_marker <- read.table("/home/jules/Documents/phd/Data/Article_veranika/full_marker.csv", sep = ",", header = TRUE)
full_marker <- full_marker$x

wrong_genes <- setdiff(full_marker, rownames(counts))

# FGF15 ???
wrong_genes_name <- c("BCL11B", "NR2F2", "CUL3", "PPP1R1B", "none", "none", "GAD1", "TPRA1", "GSX2", "MNX1", "IKZF1", "EGR2", "ZNF503", "POU3F2", "RAX", "EOMES", "VIM")

full_marker[which(full_marker %in% wrong_genes)] <- wrong_genes_name

# Get rows corresponding to markers
retained_row <- counts[rownames(counts) %in% full_marker, ]

# filtering out lowly expressed genes
countsf <- counts[rowSums(counts) > ncol(counts), ]

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
rank <- norm %>% apply(2, percent_rank)
rownames(rank) <- rownames(norm)
colnames(rank) <- colnames(norm)

markersf <- full_marker[which(full_marker %in% rownames(norm))]

mat <- norm[markersf, ]

mat <- mat[, meta[order(meta$type), ]$sample]


Heatmap(
    mat[sort(rownames(mat))[1:60], ],
    name = "Normalized expression",
    row_title_gp = gpar(fontsize = 16, fontface = "bold"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_side = "left",
    show_column_names = TRUE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    width = ncol(mat) * unit(5, "mm"),
    height = nrow(mat) * unit(2, "mm"),
    col = colorRampPalette(c(
        "darkblue",
        "lightblue",
        "white",
        "red",
        "darkred"
    ))(1000),
)

Heatmap(
    mat[sort(rownames(mat))[61:121], ],
    name = "Normalized expression",
    row_title_gp = gpar(fontsize = 16, fontface = "bold"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_side = "left",
    show_column_names = TRUE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    width = ncol(mat) * unit(5, "mm"),
    height = nrow(mat) * unit(2, "mm"),
    col = colorRampPalette(c(
        "darkblue",
        "lightblue",
        "white",
        "red",
        "darkred"
    ))(1000),
)
Heatmap(
    mat[sort(rownames(mat))[122:182], ],
    name = "Normalized expression",
    row_title_gp = gpar(fontsize = 16, fontface = "bold"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_side = "left",
    show_column_names = TRUE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    width = ncol(mat) * unit(5, "mm"),
    height = nrow(mat) * unit(2, "mm"),
    col = colorRampPalette(c(
        "darkblue",
        "lightblue",
        "white",
        "red",
        "darkred"
    ))(1000),
)
