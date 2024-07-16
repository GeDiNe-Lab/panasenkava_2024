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
rawcounts <- readcounts("/home/jules/Documents/phd/Data/lab_RNAseq/diff13/diff13_counts.csv", sep = ",", header = TRUE)
meta <- read.table("/home/jules/Documents/phd/Data/lab_RNAseq/diff13/diff13_meta.csv", sep = ",", header = TRUE)

# Keeping only necessary samples
meta <- filter(meta, type %in% c("dorsal", "ventral"), day == "day12", sample != "LON7171_D12_2")
meta %>% View()
counts <- rawcounts[, meta$sample]
counts <- counts[which(rowSums(counts) > 50), ]

counts_LON71 <- counts[, filter(meta, line == "LON71")$sample] + 1
counts_LON80 <- counts[, filter(meta, line == "LON80")$sample] + 1
counts_WTC <- counts[, filter(meta, line == "WTC")$sample] + 1

# Compute dorso_ventral DEGs
DEGs_DV_LON71 <- DESeqDataSetFromMatrix(
    countData = counts_LON71,
    colData = filter(meta, line == "LON71"),
    design = ~type
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "dorsal", "ventral")) %>%
    as.data.frame() %>%
    na.omit()

DEGs_DV_LON80 <- DESeqDataSetFromMatrix(
    countData = counts_LON80,
    colData = filter(meta, line == "LON80"),
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

DEGs_DV_LON71_f <- filter(DEGs_DV_LON71, padj < 0.05, abs(log2FoldChange) > 1)
DEGs_DV_LON71_f$gene <- gene_converter(rownames(DEGs_DV_LON71_f), "ENSEMBL", "SYMBOL")
DEGs_DV_LON71_f <- filter(DEGs_DV_LON71_f, !is.na(gene))

DEGs_DV_LON80_f <- filter(DEGs_DV_LON80, padj < 0.05, abs(log2FoldChange) > 1)
DEGs_DV_LON80_f$gene <- gene_converter(rownames(DEGs_DV_LON80_f), "ENSEMBL", "SYMBOL")
DEGs_DV_LON80_f <- filter(DEGs_DV_LON80_f, !is.na(gene))

DEGs_DV_WTC_f <- filter(DEGs_DV_WTC, padj < 0.05, abs(log2FoldChange) > 1)
DEGs_DV_WTC_f$gene <- gene_converter(rownames(DEGs_DV_WTC_f), "ENSEMBL", "SYMBOL")
DEGs_DV_WTC_f <- filter(DEGs_DV_WTC_f, !is.na(gene))

VennDiagram::venn.diagram(
    x = list(DEGs_DV_LON71_f$gene, DEGs_DV_LON80_f$gene, DEGs_DV_WTC_f$gene),
    category.names = c("LON71", "LON80", "WTC"),
    filename = "results/images/F1_4_venn_test.png",
    output = TRUE
)

Reduce(intersect, list(DEGs_DV_LON71_f$gene, DEGs_DV_LON80_f$gene, DEGs_DV_WTC_f$gene)) %>% sort()
View(DEGs_DV_WTC_f)
View(DEGs_DV_LON80_f)
