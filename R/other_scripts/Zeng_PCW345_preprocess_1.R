library(Matrix)
library(tidyverse)
library(Seurat)
source("/home/jules/Documents/phd/custom_fct.R")

meta <- read.table("/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/indices.csv", header = TRUE, sep = ",")
genes <- read.table("/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/genes.tsv", header = FALSE)$V1
meta$week_stage %>% unique()
# preprocess week3
week3 <- readMM("/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/week3.mtx") %>% t()

meta_week3 <- filter(meta, week_stage == "W3-1")
rownames(week3) <- genes
colnames(week3) <- meta_week3$index


seurat_week3 <- CreateSeuratObject(counts = week3[which(rowSums(week3) > 0), ], min.cells = 3, min.features = 200)
# seurat_week3[["percent.mt"]] <- PercentageFeatureSet(seurat_week3, pattern = "^MT-")
# VlnPlot(seurat_week3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seurat_week3 <- subset(seurat_week3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000)
week3 <- seurat_week3@assays$RNA$counts
meta_week3 <- filter(meta_week3, index %in% colnames(week3))
rm(seurat_week3)
gc()
write.csv(rownames(week3), "/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/genes_w3.csv")
write.csv(meta_week3, "/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/meta_week3.csv")
writeMM(week3, "/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/week3_filtered.mtx")
rm(week3)
rm(meta_week3)
gc()

# preprocess week3
week4 <- readMM("/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/week4.mtx") %>% t()

meta_week4 <- filter(meta, week_stage %in% c("W4-2", "W4-3"))
rownames(week4) <- genes
colnames(week4) <- meta_week4$index

seurat_week4 <- CreateSeuratObject(counts = week4[which(rowSums(week4) > 0), ], min.cells = 3, min.features = 200)

seurat_week4 <- subset(seurat_week4, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000)
week4 <- seurat_week4@assays$RNA$counts
ncol(week4)
meta_week4 %>% nrow()
meta_week4 <- filter(meta_week4, index %in% colnames(week4))
rm(seurat_week4)
gc()
write.csv(rownames(week4), "/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/genes_w4.csv")
write.csv(meta_week4, "/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/meta_week4.csv")
writeMM(week4, "/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/week4_filtered.mtx")
rm(week4)
rm(meta_week4)
gc()

# preprocess week3
week5 <- readMM("/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/week5.mtx") %>% t()

meta_week5 <- filter(meta, week_stage %in% c("W5-2", "W5-3"))
rownames(week5) <- genes
colnames(week5) <- meta_week5$index


seurat_week5 <- CreateSeuratObject(counts = week5[which(rowSums(week5) > 0), ], min.cells = 3, min.features = 200)

seurat_week5 <- subset(seurat_week5, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000)
week5 <- seurat_week5@assays$RNA$counts
meta_week5 <- filter(meta_week5, index %in% colnames(week5))
rm(seurat_week5)
gc()
write.csv(rownames(week5), "/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/genes_w5.csv")
write.csv(meta_week5, "/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/meta_week5.csv")
writeMM(week5, "/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/week5_filtered.mtx")
rm(week5)
rm(meta_week5)
gc()
