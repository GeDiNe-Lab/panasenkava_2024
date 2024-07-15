# Loading packages and functions
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(DESeq2)

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

# LON71_D12_2 is a biiiiiig outlier
meta <- filter(meta, sample != "LON71_D12_2", type %in% c("ventral", "dorsal"))
counts <- rawcounts[which(rowSums(rawcounts) >= 50), meta$sample]

dds <- DESeqDataSetFromMatrix(
    countData = counts[, meta$sample],
    colData = meta,
    design = ~ type + day
)

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# PCA plot
pca.data <- plotPCA(vsd, intgroup = c("type", "day", "line"), returnData = TRUE)
percentVar <- round(100 * attr(pca.data, "percentVar"))

png(filename = "results/images/Figure_2A/F2A_1_PCA_type.png", width = 1600, height = 1200, res = 250)
ggplot(pca.data, aes(PC1, PC2, color = type, shape = day)) +
    geom_point(size = 2, stroke = 1) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_color_manual(values = c("#A1A1DE", "#80AD3C")) +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6)) +
    custom_theme() +
    ggtitle("PCA of dorsal and ventral kinetics")
dev.off()

png(filename = "results/images/Figure_2A/F2A_1_PCA_line.png", width = 1600, height = 1200, res = 250)
ggplot(pca.data, aes(PC1, PC2, color = line, shape = day)) +
    geom_point(size = 2, stroke = 1) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6)) +
    custom_theme() +
    ggtitle("PCA of dorsal and ventral kinetics")
dev.off()
