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
rawcounts <- readcounts("/home/jules/Documents/phd/Data/Article_veranika/bulk/counts.csv", sep = ",", header = TRUE)
meta <- read.table("/home/jules/Documents/phd/Data/Article_veranika/bulk/metadata.csv", sep = ",", header = TRUE)

# Keeping only necessary samples
meta <- filter(meta, type %in% c("dorsal", "ventral", "ipsc"), CRISPR %in% c("ipsc", "control"))

counts <- rawcounts[, meta$sample]
counts <- counts[which(rowSums(counts) > ncol(counts)), ]
norm <- varianceStabilizingTransformation(counts)
norm <- limma::removeBatchEffect(norm, meta$line)

pca_results <- ggPCA(t(norm), ncp = 5, graph = FALSE, scale.unit = TRUE)

meta_pca <- cbind(meta, pca_results$gg.ind)

png(filename = "results/images/F1_2_PCA.png", width = 2100, height = 1600, res = 250)
ggplot(data = meta_pca, aes(x = PC1, y = PC2, color = type, shape = line)) +
    geom_point() +
    custom_theme()
dev.off()
