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

counts <- rawcounts[which(rowSums(rawcounts) > 50), meta$sample]

counts_ventral <- rawcounts[, filter(meta, type == "ventral")$sample]
counts_ventral <- counts_ventral[which(rowSums(counts_ventral) > 50), ]

counts_dorsal <- rawcounts[, filter(meta, type == "dorsal")$sample]
counts_dorsal <- counts_dorsal[which(rowSums(counts_dorsal) > 50), ]

counts_LON <- rawcounts[, filter(meta, line == "LON")$sample]
counts_LON <- counts_LON[which(rowSums(counts_LON) > 50), ]

counts_WTC <- rawcounts[, filter(meta, line == "WTC")$sample]
counts_WTC <- counts_WTC[which(rowSums(counts_WTC) > 50), ]

DEGs_line_ventral <- DESeqDataSetFromMatrix(
    countData = counts_ventral,
    colData = filter(meta, type == "ventral"),
    design = ~line
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("line", "LON", "WTC")) %>%
    as.data.frame() %>%
    na.omit()

DEGs_line_dorsal <- DESeqDataSetFromMatrix(
    countData = counts_dorsal,
    colData = filter(meta, type == "dorsal"),
    design = ~line
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("line", "LON", "WTC")) %>%
    as.data.frame() %>%
    na.omit()

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

DEGs_line_ventral$gene <- gene_converter(rownames(DEGs_line_ventral), "ENSEMBL", "SYMBOL")
DEGs_line_ventral <- filter(DEGs_line_ventral, !is.na(gene))
DEGs_line_ventral_f <- filter(DEGs_line_ventral, padj < 0.05, abs(log2FoldChange) > 0.5)

DEGs_line_dorsal$gene <- gene_converter(rownames(DEGs_line_dorsal), "ENSEMBL", "SYMBOL")
DEGs_line_dorsal <- filter(DEGs_line_dorsal, !is.na(gene))
DEGs_line_dorsal_f <- filter(DEGs_line_dorsal, padj < 0.05, abs(log2FoldChange) > 0.5)

DEGs_DV_LON$gene <- gene_converter(rownames(DEGs_DV_LON), "ENSEMBL", "SYMBOL")
DEGs_DV_LON <- filter(DEGs_DV_LON, !is.na(gene))
DEGs_DV_LON_f <- filter(DEGs_DV_LON, padj < 0.05, abs(log2FoldChange) > 0.5)

DEGs_DV_WTC$gene <- gene_converter(rownames(DEGs_DV_WTC), "ENSEMBL", "SYMBOL")
DEGs_DV_WTC <- filter(DEGs_DV_WTC, !is.na(gene))
DEGs_DV_WTC_f <- filter(DEGs_DV_WTC, padj < 0.05, abs(log2FoldChange) > 0.5)

write.csv(DEGs_line_ventral, "results/tables/Figure_1/LON_VS_WTC_DEGs_ventral.csv")
write.csv(DEGs_line_dorsal_f, "results/tables/Figure_1/LON_VS_WTC_DEGs_dorsal.csv")
write.csv(DEGs_DV_LON_f, "results/tables/Figure_1/dorsal_VS_ventral_DEGs_LON.csv")
write.csv(DEGs_DV_WTC_f, "results/tables/Figure_1/dorsal_VS_ventral_DEGs_WTC.csv")

# 0.5 FC Venn diagram
VennDiagram::venn.diagram(
    x = list(DEGs_DV_LON_f$gene, DEGs_DV_WTC_f$gene),
    main = paste0(
        "Dorsal VS ventral ",
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
        "Dorsal VS ventral ",
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
        "Dorsal VS ventral ",
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
        "Dorsal VS ventral ",
        round((length(intersect(filter(DEGs_DV_LON_f, abs(log2FoldChange) >= 2)$gene, filter(DEGs_DV_WTC_f, abs(log2FoldChange) >= 2)$gene)) / length(union(filter(DEGs_DV_LON_f, abs(log2FoldChange) >= 2)$gene, filter(DEGs_DV_WTC_f, abs(log2FoldChange) >= 2)$gene))) * 100, 2),
        "% of common Genes above 2 absolute log2FC"
    ),
    category.names = c("LON", "WTC"),
    filename = "results/images/Figure_1/F1_4_LON_WTC_DE_FC2_venn.png",
    output = TRUE,
    disable.logging = TRUE,
)
