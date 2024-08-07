#  Loading packages and functions
library(Matrix)
library(Seurat)
library(data.table)
library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(DESeq2)
library(ggrepel)

rstudioapi::getSourceEditorContext()$path %>%
    str_split("/") %>%
    unlist() %>%
    head(-3) %>%
    str_c(collapse = "/") %>%
    str_c("/") %>%
    setwd()

source("R/custom_fct.R")

rawcounts <- readcounts("/home/jules/Documents/phd/Data/lab_RNAseq/diff13/diff13_counts.csv", sep = ",", header = TRUE)
meta <- read.table("/home/jules/Documents/phd/Data/lab_RNAseq/diff13/diff13_meta.csv", sep = ",", header = TRUE)

# LON71_D12_2 is a biiiiiig outlier
meta <- filter(meta, sample != "LON71_D12_2", type %in% c("ventral", "dorsal"), line %in% c("LON71", "WTC"))
counts <- rawcounts[which(rowSums(rawcounts) >= 50), meta$sample]

norm <- varianceStabilizingTransformation(counts)

# Get the path to the downloaded data
data_dir <- "/home/jules/Documents/phd/Data/literature/Eze_et_al/download"

# Initialize lists to store data and metadata
data_list <- list()
metadata_list <- list()

# Iterate through the folders
folders <- list.files(data_dir, full.names = TRUE)

for (folder in folders) {
    # Extract information from the folder name
    brain_region <- str_split(basename(folder), "_")[[1]][2]
    stage <- str_split(basename(folder), "_")[[1]][1]
    print(brain_region)
    # Read the matrix.mtx file
    matrix_file <- file.path(folder, "matrix.mtx")
    mtx <- readMM(matrix_file)

    genes_file <- file.path(folder, "genes.tsv")
    genes <- read.table(genes_file, header = FALSE, stringsAsFactors = FALSE)

    # Read the barcodes.tsv
    barcodes_file <- file.path(folder, "barcodes.tsv")
    barcodes <- read.table(barcodes_file, header = FALSE, stringsAsFactors = FALSE)

    # set colnames and rownames
    rownames(mtx) <- genes$V1
    colnames(mtx) <- paste(barcodes$V1, brain_region, sep = "_")

    # Build metadata dataframe for this folder
    curr_meta <- data.frame(
        cell = colnames(mtx),
        brain_region = rep(brain_region, ncol(mtx)),
        stage = rep(stage, ncol(mtx))
    )
    metadata_list[[basename(folder)]] <- curr_meta
    data_list[[basename(folder)]] <- mtx
    gc()
}
# Aggregate data into one metadata dataframe and one counts matrix
metadata <- do.call(rbind, metadata_list)
rm(metadata_list)
raw_sc_counts <- do.call(cbind, data_list)
rm(data_list)
gc()
# Diencephalon forebrain prosencephalon Telencephalon
# check for cells order:
ant_metadata <- filter(metadata, brain_region %in% c("Diencephalon", "forebrain", "forebrain2", "prosencephalon", "Telencephalon"))
ant_sc_counts <- raw_sc_counts[, ant_metadata$cell]
table(ant_metadata$stage, ant_metadata$brain_region)

# filter out dropout genes
ant_sc_counts <- ant_sc_counts[which(rowSums(ant_sc_counts) > 0), ]
# Clean single cell data using Seurat workflow
# throwing out cells with less than 200 unique features count
seurat_obj <- CreateSeuratObject(counts = CreateAssayObject(ant_sc_counts,
    assay = "RNA",
    min.features = 200,
    min.cells = 20
))

# Get mitochondrial genes expression %
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
# there is 0 expression of mithocondrial genes

# vizualizing QCs
# VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
seurat_obj <- subset(seurat_obj, nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
seurat_obj <- NormalizeData(seurat_obj)
ant_sc_norm <- seurat_obj@assays$RNA$data %>% as("dgCMatrix")

ant_metadata <- ant_metadata %>% dplyr::filter(cell %in% colnames(ant_sc_norm))

common_genes <- intersect(rownames(ant_sc_norm), rownames(norm))
CS12_cor <- WGCNA::cor(ant_sc_norm[common_genes, filter(ant_metadata, stage == "CS12")$cell], norm[common_genes, ])
CS12_cor[which(is.na(CS12_cor))] <- 0
meta$CS12_mean <- CS12_cor %>% colMeans()
meta$CS12_sd <- CS12_cor %>% colSds()

CS13_cor <- WGCNA::cor(ant_sc_norm[common_genes, filter(ant_metadata, stage == "CS13")$cell], norm[common_genes, ])
CS13_cor[which(is.na(CS13_cor))] <- 0
meta$CS13_mean <- CS13_cor %>% colMeans()
meta$CS13_sd <- CS13_cor %>% colSds()

CS15_cor <- WGCNA::cor(ant_sc_norm[common_genes, filter(ant_metadata, stage == "CS15")$cell], norm[common_genes, ])
CS15_cor[which(is.na(CS15_cor))] <- 0
meta$CS15_mean <- CS15_cor %>% colMeans()
meta$CS15_sd <- CS15_cor %>% colSds()

ggplot(meta, aes(x = day, y = CS12_mean, color = type, group = type)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = CS12_mean - CS12_sd, ymax = CS12_mean + CS12_sd), width = 0.2) +
    ylim(c(-0.1, 1)) +
    custom_theme()

ggplot(meta, aes(x = day, y = CS13_mean, color = type, group = type)) +
    geom_point() +
    geom_line() +
    ylim(c(-0.1, 1)) +
    geom_errorbar(aes(ymin = CS13_mean - CS13_sd, ymax = CS13_mean + CS13_sd), width = 0.2) +
    custom_theme()

ggplot(meta, aes(x = day, y = CS15_mean, color = type, group = type)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = CS15_mean - CS15_sd, ymax = CS15_mean + CS15_sd), width = 0.2) +
    ylim(c(-0.1, 1)) +
    custom_theme()


CS13metadata <- filter(metadata, stage == "CS13")
CS13sc_counts <- raw_sc_counts[, CS13metadata$cell]

# filter out dropout genes
CS13sc_counts <- CS13sc_counts[which(rowSums(CS13sc_counts) > 0), ]
# Clean single cell data using Seurat workflow
# throwing out cells with less than 200 unique features count
seurat_obj <- CreateSeuratObject(counts = CreateAssayObject(CS13sc_counts,
    assay = "RNA",
    min.features = 200,
    min.cells = 20
))

VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

seurat_obj <- subset(seurat_obj, nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000)
seurat_obj <- NormalizeData(seurat_obj)



Idents(seurat_obj) <- "seurat_annotations"
anteriorMB <- FindMarkers(seurat_obj, ident.1 = "CD16 Mono", ident.2 = "CD14 Mono")



CS13sc_norm <- seurat_obj@assays$RNA$data %>% as("dgCMatrix")

CS13metadata <- CS13metadata %>% dplyr::filter(cell %in% colnames(CS13sc_norm))
CS13metadata$brain_region %>% table()

common_genes <- intersect(VariableFeatures(seurat_obj), rownames(norm))
anteriorMB_cor <- WGCNA::cor(CS13sc_norm[common_genes, filter(CS13metadata, brain_region == "anterior")$cell], norm[common_genes, ])
anteriorMB_cor[which(is.na(anteriorMB_cor))] <- 0
meta$anteriorMB_mean <- anteriorMB_cor %>% colMeans()
meta$anteriorMB_sd <- anteriorMB_cor %>% colSds()

centerMB_corr <- WGCNA::cor(CS13sc_norm[common_genes, filter(CS13metadata, brain_region == "center")$cell], norm[common_genes, ])
centerMB_corr[which(is.na(centerMB_corr))] <- 0
meta$centerMB_mean <- centerMB_corr %>% colMeans()
meta$centerMB_sd <- centerMB_corr %>% colSds()

hindbrain_corr <- WGCNA::cor(CS13sc_norm[common_genes, filter(CS13metadata, brain_region == "hindbrain")$cell], norm[common_genes, ])
hindbrain_corr[which(is.na(hindbrain_corr))] <- 0
meta$hindbrain_mean <- hindbrain_corr %>% colMeans()
meta$hindbrain_sd <- hindbrain_corr %>% colSds()

prosencephalon_corr <- WGCNA::cor(CS13sc_norm[common_genes, filter(CS13metadata, brain_region == "prosencephalon")$cell], norm[common_genes, ])
prosencephalon_corr[which(is.na(prosencephalon_corr))] <- 0
meta$prosencephalon_mean <- prosencephalon_corr %>% colMeans()
meta$prosencephalon_sd <- prosencephalon_corr %>% colSds()

ggplot(meta, aes(x = day, y = hindbrain_mean, color = type, group = type)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = hindbrain_mean - hindbrain_sd, ymax = hindbrain_mean + hindbrain_sd), width = 0.2) +
    ylim(c(-0.1, 1)) +
    custom_theme()

ggplot(meta, aes(x = day, y = centerMB_mean, color = type, group = type)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = centerMB_mean - centerMB_sd, ymax = centerMB_mean + centerMB_sd), width = 0.2) +
    ylim(c(-0.1, 1)) +
    custom_theme()

ggplot(meta, aes(x = day, y = anteriorMB_mean, color = type, group = type)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = anteriorMB_mean - anteriorMB_sd, ymax = anteriorMB_mean + anteriorMB_sd), width = 0.2) +
    ylim(c(-0.1, 1)) +
    custom_theme()

ggplot(meta, aes(x = day, y = prosencephalon_mean, color = type, group = type)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = prosencephalon_mean - prosencephalon_sd, ymax = prosencephalon_mean + prosencephalon_sd), width = 0.2) +
    ylim(c(-0.1, 1)) +
    custom_theme()

CS13metadata$brain_region %>% table()
