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
rawcounts <- readcounts("/home/jules/Documents/phd/Data/lab_RNAseq/diff13/diff13_counts.csv", sep = ",", header = TRUE)
meta <- read.table("/home/jules/Documents/phd/Data/lab_RNAseq/diff13/diff13_meta.csv", sep = ",", header = TRUE)

# Keeping only necessary samples
meta <- filter(meta, type %in% c("dorsal", "ventral"), day == "day12", sample != "LON71_D12_2")

counts <- rawcounts[, meta$sample]
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
markers_bis <- read.table("/home/jules/Documents/phd/Data/Article_veranika/GeneListFig1_15_07_24.csv", sep = ";", header = FALSE)
markers_bis
markers_bis <- markers_bis$V1
setdiff(markers_bis, rownames(counts))
markers_bis <- intersect(markers_bis, rownames(counts))
# Get rows corresponding to markers
retained_row <- counts[rownames(counts) %in% markers_bis, ]

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
# norm <- limma::removeBatchEffect(norm, meta$line)
# scaled_norm <- t(apply(norm, 1, scale))
# colnames(scaled_norm) <- colnames(norm)
# rownames(scaled_norm) <- rownames(norm)
markers_bis %>% length()
png(filename = "results/images/F1_1_marker_HM_diff13.png", width = 2000, height = 1800, res = 250)
Heatmap(
    norm[markers_bis, meta[order(meta$type), ]$sample],
    name = "Normalized expression",
    row_title_gp = gpar(fontsize = 16, fontface = "bold"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_side = "left",
    show_column_names = TRUE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    width = ncol(norm[markers_bis, ]) * unit(4, "mm"),
    height = nrow(norm[markers_bis, ]) * unit(4, "mm"),
    col = colorRampPalette(c(
        "blue",
        "white",
        "red"
    ))(1000),
)
dev.off()
