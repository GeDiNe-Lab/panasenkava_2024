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
meta <- filter(meta, type %in% c("dorsal", "ventral", "ipsc"), CRISPR %in% c("ipsc", "control"))

counts <- rawcounts[, meta$sample]
counts <- counts[which(rowSums(counts) > 50), ]

counts_LON <- counts[, filter(meta, line == "LON")$sample]
counts_WTC <- counts[, filter(meta, line == "WTC")$sample]

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

DEGs_DV_LON_f <- filter(DEGs_DV_LON, padj < 0.05, abs(log2FoldChange) > 1)
DEGs_DV_LON_f$gene <- gene_converter(rownames(DEGs_DV_LON_f), "ENSEMBL", "SYMBOL")
DEGs_DV_LON_f <- filter(DEGs_DV_LON_f, !is.na(gene))

DEGs_DV_WTC_f <- filter(DEGs_DV_WTC, padj < 0.05, abs(log2FoldChange) > 1)
DEGs_DV_WTC_f$gene <- gene_converter(rownames(DEGs_DV_WTC_f), "ENSEMBL", "SYMBOL")
DEGs_DV_WTC_f <- filter(DEGs_DV_WTC_f, !is.na(gene))
