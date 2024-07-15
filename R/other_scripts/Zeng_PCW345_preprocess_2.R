library(Matrix)
library(tidyverse)
library(Seurat)
source("/home/jules/Documents/phd/custom_fct.R")

meta_w3 <- read.table("/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/meta_week3.csv", header = TRUE, sep = ",")
genes_w3 <- read.table("/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/genes_w3.csv", header = TRUE, sep = ",")$x
week3 <- readMM("/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/week3_filtered.mtx")
colnames(week3) <- meta_w3$index
rownames(week3) <- genes_w3

meta_w4 <- read.table("/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/meta_week4.csv", header = TRUE, sep = ",")
genes_w4 <- read.table("/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/genes_w4.csv", header = TRUE, sep = ",")$x
week4 <- readMM("/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/week4_filtered.mtx")
colnames(week4) <- meta_w4$index
rownames(week4) <- genes_w4

meta_w5 <- read.table("/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/meta_week5.csv", header = TRUE, sep = ",")
genes_w5 <- read.table("/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/genes_w5.csv", header = TRUE, sep = ",")$x
week5 <- readMM("/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/week5_filtered.mtx")
colnames(week5) <- meta_w5$index
rownames(week5) <- genes_w5

common_genes <- Reduce(intersect, list(genes_w3, genes_w4, genes_w5))

counts <- Reduce(cbind, list(week5[common_genes, ], week4[common_genes, ], week3[common_genes, ]))
meta <- Reduce(rbind, list(meta_w5, meta_w4, meta_w3))
rm(week3, week4, week5)
rm(meta_w3, meta_w4, meta_w5)
gc()

write.csv(rownames(counts), "/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/genes_final.csv")
write.csv(meta, "/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/meta_final.csv")
writeMM(counts, "/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/counts_w345.mtx")
