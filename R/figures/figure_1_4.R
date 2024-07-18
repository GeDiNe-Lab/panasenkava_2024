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
# Loading data (path to change later)
rawcounts <- readcounts("/home/jules/Documents/phd/Data/Article_veranika/bulk/counts.csv", sep = ",", header = TRUE)
meta <- read.table("/home/jules/Documents/phd/Data/Article_veranika/bulk/metadata.csv", sep = ",", header = TRUE)
# Keeping only necessary samples
meta <- filter(meta, type %in% c("dorsal", "ventral"), CRISPR %in% c("control"))

counts <- rawcounts[, meta$sample]

counts_LON <- counts[, filter(meta, line == "LON")$sample]
counts_LON <- counts_LON[which(rowSums(counts_LON) > 50), ]

counts_WTC <- counts[, filter(meta, line == "WTC")$sample]
counts_WTC <- counts_WTC[which(rowSums(counts_WTC) > 50), ]

# Compute dorso_ventral DEGs
DEGs_DV_LON <- DESeqDataSetFromMatrix(
    countData = counts_LON,
    colData = filter(meta, line == "LON"),
    design = ~type
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "dorsal", "ventral")) %>%
    as.data.frame() %>%
    na.omit()

DEGs_DV_WTC <- DESeqDataSetFromMatrix(
    countData = counts_WTC,
    colData = filter(meta, line == "WTC"),
    design = ~type
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "dorsal", "ventral")) %>%
    as.data.frame() %>%
    na.omit()

DEGs_DV_LON$gene <- gene_converter(rownames(DEGs_DV_LON), "ENSEMBL", "SYMBOL")
DEGs_DV_LON <- filter(DEGs_DV_LON, !is.na(gene))
DEGs_DV_LON_f <- filter(DEGs_DV_LON, padj < 0.05, abs(log2FoldChange) > 0.5)

DEGs_DV_WTC$gene <- gene_converter(rownames(DEGs_DV_WTC), "ENSEMBL", "SYMBOL")
DEGs_DV_WTC <- filter(DEGs_DV_WTC, !is.na(gene))
DEGs_DV_WTC_f <- filter(DEGs_DV_WTC, padj < 0.05, abs(log2FoldChange) > 0.5)


write.csv(DEGs_DV_LON_f, "results/tables/Figure_1/dorsal_VS_ventral_DEGs_LON.csv")
write.csv(DEGs_DV_WTC_f, "results/tables/Figure_1/dorsal_VS_ventral_DEGs_WTC.csv")

# 0.5 FC Venn diagram

VennDiagram::venn.diagram(
    x = list(DEGs_DV_LON_f$gene, DEGs_DV_WTC_f$gene),
    main = paste0(
        round((length(intersect(DEGs_DV_LON_f$gene, DEGs_DV_WTC_f$gene)) / length(union(DEGs_DV_LON_f$gene, DEGs_DV_WTC_f$gene))) * 100, 2),
        "% of common Genes above 0.5 absolute log2FC"
    ),
    category.names = c("LON", "WTC"),
    filename = "results/images/Figure_1/F1_4_LON_WTC_DE_FC0_5_venn.png",
    output = TRUE,
    disable.logging = TRUE
)

# 1 FC Venn diagram
VennDiagram::venn.diagram(
    x = list(filter(DEGs_DV_LON_f, abs(log2FoldChange) >= 1)$gene, filter(DEGs_DV_WTC_f, abs(log2FoldChange) >= 1)$gene),
    main = paste0(
        round((length(intersect(filter(DEGs_DV_LON_f, abs(log2FoldChange) >= 1)$gene, filter(DEGs_DV_WTC_f, abs(log2FoldChange) >= 1)$gene)) / length(union(filter(DEGs_DV_LON_f, abs(log2FoldChange) >= 1)$gene, filter(DEGs_DV_WTC_f, abs(log2FoldChange) >= 1)$gene))) * 100, 2),
        "% of common Genes above 1 absolute log2FC"
    ),
    category.names = c("LON", "WTC"),
    filename = "results/images/Figure_1/F1_4_LON_WTC_DE_FC1_venn.png",
    output = TRUE,
    disable.logging = TRUE,
)

# 1.5 FC Venn diagram
VennDiagram::venn.diagram(
    x = list(filter(DEGs_DV_LON_f, abs(log2FoldChange) >= 1.5)$gene, filter(DEGs_DV_WTC_f, abs(log2FoldChange) >= 1.5)$gene),
    main = paste0(
        round((length(intersect(filter(DEGs_DV_LON_f, abs(log2FoldChange) >= 1.5)$gene, filter(DEGs_DV_WTC_f, abs(log2FoldChange) >= 1.5)$gene)) / length(union(filter(DEGs_DV_LON_f, abs(log2FoldChange) >= 1.5)$gene, filter(DEGs_DV_WTC_f, abs(log2FoldChange) >= 1.5)$gene))) * 100, 2),
        "% of common Genes above 1.5 absolute log2FC"
    ),
    category.names = c("LON", "WTC"),
    filename = "results/images/Figure_1/F1_4_LON_WTC_DE_FC1_5_venn.png",
    output = TRUE,
    disable.logging = TRUE,
)

# 2 FC Venn diagram
VennDiagram::venn.diagram(
    x = list(filter(DEGs_DV_LON_f, abs(log2FoldChange) >= 2)$gene, filter(DEGs_DV_WTC_f, abs(log2FoldChange) >= 2)$gene),
    main = paste0(
        round((length(intersect(filter(DEGs_DV_LON_f, abs(log2FoldChange) >= 2)$gene, filter(DEGs_DV_WTC_f, abs(log2FoldChange) >= 2)$gene)) / length(union(filter(DEGs_DV_LON_f, abs(log2FoldChange) >= 2)$gene, filter(DEGs_DV_WTC_f, abs(log2FoldChange) >= 2)$gene))) * 100, 2),
        "% of common Genes above 2 absolute log2FC"
    ),
    category.names = c("LON", "WTC"),
    filename = "results/images/Figure_1/F1_4_LON_WTC_DE_FC2_venn.png",
    output = TRUE,
    disable.logging = TRUE,
)
