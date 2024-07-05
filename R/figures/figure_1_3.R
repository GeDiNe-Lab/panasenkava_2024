# Loading packages and functions
library(Matrix)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(org.Hs.eg.db)
library(ComplexHeatmap)

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
rawcounts <- readcounts("/home/jules/Documents/phd/Data/Article_veranika/bulk/counts.csv", sep = ",", header = TRUE)
meta <- read.table("/home/jules/Documents/phd/Data/Article_veranika/bulk/metadata.csv", sep = ",", header = TRUE)
# Keeping only necessary samples
meta <- filter(meta, type %in% c("dorsal", "ventral", "ipsc"), CRISPR %in% c("ipsc", "control"))

counts <- rawcounts[, meta$sample]
counts <- counts[which(rowSums(counts) > 50), ]
norm <- varianceStabilizingTransformation(counts)
norm <- limma::removeBatchEffect(norm, meta$line)

# Compute dorso_ventral DEGs
DEGs_DV <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ line + type
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "dorsal", "ventral")) %>%
    as.data.frame() %>%
    na.omit()

DEGs_DV_f <- filter(DEGs_DV, padj < 0.01, abs(log2FoldChange) > 2)
DEGs_DV_f$gene <- gene_converter(rownames(DEGs_DV_f), "ENSEMBL", "SYMBOL")
DEGs_DV_f <- filter(DEGs_DV_f, !is.na(gene))

write.csv(DEGs_DV_f, "results/tables/dorsal_VS_ventral_DEGs.csv")

norm_DEGs <- norm[rownames(DEGs_DV_f), ]
scaled_mat <- t(apply(norm_DEGs, 1, scale))
colnames(scaled_mat) <- colnames(norm_DEGs)


png(filename = "results/images/F1_3_DE_HM.png", width = 1600, height = 1600, res = 250)
Heatmap(
    scaled_mat,
    name = "Normalized expression",
    row_title_gp = gpar(fontsize = 16, fontface = "bold"),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    row_names_side = "left",
    show_column_names = TRUE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    width = ncol(scaled_mat) * unit(4, "mm"),
    # height = nrow(mat) * unit(5, "mm"),
    col = colorRampPalette(c(
        "blue",
        "white",
        "red"
    ))(1000),
)
dev.off()
