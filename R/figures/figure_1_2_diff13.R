# Loading packages and functions
library(Matrix)
library(tidyverse)
library(DESeq2)
library(ggplot2)
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
rawcounts <- readcounts("/home/jules/Documents/phd/Data/lab_RNAseq/diff13/diff13_counts.csv", sep = ",", header = TRUE)
meta <- read.table("/home/jules/Documents/phd/Data/lab_RNAseq/diff13/diff13_meta.csv", sep = ",", header = TRUE)

# Keeping only necessary samples
meta <- filter(meta, type %in% c("dorsal", "ventral"), day == "day12", sample != "LON71_D12_2")

counts <- rawcounts[, meta$sample]
counts <- counts[which(rowSums(counts) >= 50), ]
View(meta)
dds_PCA <- DESeqDataSetFromMatrix(
    countData = counts[, meta$sample],
    colData = meta,
    design = ~ sexe + type
)
# Variance stabilizing transformation
vsd <- vst(dds_PCA, blind = FALSE)

# PCA plot
pca.data <- plotPCA(vsd, intgroup = c("type", "line"), returnData = TRUE)
percentVar <- round(100 * attr(pca.data, "percentVar"))

png(filename = "results/images/Figure_1/F1_2_PCA_diff13.png", width = 1600, height = 1200, res = 250)
ggplot(pca.data, aes(PC1, PC2, color = type, shape = line)) +
    geom_point(size = 2, stroke = 1) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_color_manual(values = c("#A1A1DE", "#80AD3C")) +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6)) +
    custom_theme() +
    ggtitle("PCA of dorsal and ventral kinetics")
dev.off()
