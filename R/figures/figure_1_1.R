# Loading packages and functions
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

markers <- c(
    "OTX2",
    "SIX3",
    "BMP7",
    "ARX",
    "LHX2",
    "SHH",
    "NKX2-1",
    "SIX6",
    "FGF10",
    "RAX",
    "NKX2-2",
    "FOXA2",
    "NTN1",
    "FOXO1",
    "SCUBE1",
    "PAX6",
    "SOX1",
    "SP8",
    "EMX2",
    "LHX5",
    "SEMA3A",
    "PAX3",
    "CNTN2",
    "MSX1",
    "OTX1",
    "EN1",
    "EN2",
    "GBX2",
    "HOXB1",
    "HOXA1"
)

# Get rows corresponding to markers
retained_row <- counts[rownames(counts) %in% markers, ]
retained_row
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

markers
mat <- norm[markers, ]

mat <- mat[, meta[order(meta$type), ]$sample]


Heatmap(
    mat[markers, ],
    name = "Normalized expression",
    row_title_gp = gpar(fontsize = 16, fontface = "bold"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_side = "left",
    show_column_names = TRUE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    width = ncol(mat) * unit(4, "mm"),
    height = nrow(mat) * unit(5, "mm"),
    col = colorRampPalette(c(
        "darkblue",
        "lightblue",
        "white",
        "red",
        "darkred"
    ))(1000),
)
