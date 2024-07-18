# Loading packages and functions
library(Matrix)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(VennDiagram)
library(tibble)

rstudioapi::getSourceEditorContext()$path %>%
    str_split("/") %>%
    unlist() %>%
    head(-3) %>%
    str_c(collapse = "/") %>%
    str_c("/") %>%
    setwd()

source("R/custom_fct.R")
rawcounts <- readcounts("/home/jules/Documents/phd/Data/lab_RNAseq/manip4/manip4_counts.csv")
meta <- read.table("/home/jules/Documents/phd/Data/lab_RNAseq/manip4/manip4_metadata.csv", sep = ",", header = T)

meta <- filter(meta, type %in% c("cyclo", "ventral"))
counts <- rawcounts[which(rowSums(rawcounts) >= 50), meta$samples]
norm <- varianceStabilizingTransformation(counts)
meta %>% View()
meta$cyclo_dose_qual <- meta$cyclo_dose %>% sapply(function(x) {
    if (x %in% c(0.125, 0.25)) {
        return("low")
    } else if (x %in% c(0.5, 1)) {
        return("high")
    } else {
        return("no_cyclo")
    }
})

DE_no_high <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~cyclo_dose_qual
) %>%
    DESeq() %>%
    results(alpha = 0.01, contrast = c("cyclo_dose_qual", "no_cyclo", "high")) %>%
    as.data.frame()
DE_no_high$gene <- gene_converter(rownames(DE_no_high), "ENSEMBL", "SYMBOL")
DE_no_high <- filter(DE_no_high, padj < 0.05 & abs(log2FoldChange) > 1 & !is.na(gene))

DE_low_high <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~cyclo_dose_qual
) %>%
    DESeq() %>%
    results(alpha = 0.01, contrast = c("cyclo_dose_qual", "low", "high")) %>%
    as.data.frame()
DE_low_high$gene <- gene_converter(rownames(DE_low_high), "ENSEMBL", "SYMBOL")
DE_low_high <- filter(DE_low_high, padj < 0.05 & abs(log2FoldChange) > 1 & !is.na(gene))

DE_no_low <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~cyclo_dose_qual
) %>%
    DESeq() %>%
    results(alpha = 0.01, contrast = c("cyclo_dose_qual", "no_cyclo", "low")) %>%
    as.data.frame()
DE_no_low$gene <- gene_converter(rownames(DE_no_low), "ENSEMBL", "SYMBOL")
DE_no_low <- filter(DE_no_low, padj < 0.05 & abs(log2FoldChange) > 1 & !is.na(gene))

write.csv(DE_no_high, "results/tables/DEG_cyclo_vAN_VS_high.csv")
write.csv(DE_low_high, "results/tables/DEG_cyclo_low_VS_high.csv")
write.csv(DE_no_low, "results/tables/DEG_cyclo_vAN_VS_low.csv")

scaled_mat <- t(apply(norm[rownames(DE_no_high)], 1, scale))
colnames(scaled_mat) <- colnames(mat)
DE_no_high %>% nrow()
# png(filename = "results/images/Figure_2A/F2A_DE_HM.png", width = 2400, height = 1600, res = 250)
Heatmap(
    scaled_mat,
    name = "Normalized expression",
    column_names_gp = gpar(fontsize = 6),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    row_names_side = "left",
    show_column_names = TRUE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    width = ncol(scaled_mat) * unit(2, "mm"),
    # height = nrow(mat) * unit(5, "mm"),
    col = colorRampPalette(c(
        "darkblue",
        "blue",
        "white",
        "red",
        "darkred"
    ))(1000),
)
dev.off()
