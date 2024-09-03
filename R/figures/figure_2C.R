# Loading packages and functions
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(DESeq2)
library(Seurat)

# Setting working directory
rstudioapi::getSourceEditorContext()$path %>%
    str_split("/") %>%
    unlist() %>%
    head(-3) %>%
    str_c(collapse = "/") %>%
    str_c("/") %>%
    setwd()

# Loading custom functions
source("R/custom_fct.R")
# Loading data (path to change later)
rawcounts <- readcounts("data/rawcounts.csv", sep = ",", header = TRUE)
rawmeta <- read.table("data/meta.csv", sep = ",", header = TRUE)

# LON71_D12_2 does not have any reads in the count file
# though, the fastQC report shows that the sample is good
meta <- filter(rawmeta, sample != "LON71_D12_2", diff == "diff13", line %in% c("LON71", "WTC"))
View(meta)
# filtering out lowly expressed genes
counts <- rawcounts[, meta$sample][which(rowSums(rawcounts[, meta$sample]) >= 25), ]

# making DESeq object with lineage,days and type as covariates
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ line + day + type
)

# Normalization by variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)
rownames(vsd) <- gene_converter(rownames(vsd), "ENSEMBL", "SYMBOL")
rownames(vsd)
# loading single cell data from Zeng et al from week3,week4 and week5
sc_counts <- readMM("scdata/week345_wholebody.mtx") %>% t()
sc_counts <- t(sc_counts)

# keeping only cells related to Central Nervous System (SNC)
sc_meta <- read.table("scdata/indices_week345_wholebody.csv", header = TRUE, sep = ",")
genes <- read.table("scdata/genes_week345_wholebody.tsv", header = FALSE)$V1


rownames(sc_counts) <- genes

formated_indices <- c(1:nrow(meta)) %>% sapply(function(i) {
    c(str_split(meta[i, ]$index, "-")[[1]][1:2], meta[i, ]$week_stage) %>%
        paste(collapse = "-") %>%
        return()
})
colnames(sc_counts) <- formated_indices

authors_meta <- read.csv("scdata/week345_wholebody_metadata.csv")
authors_meta_f <- filter(authors_meta, celltype_region_num %in% c(1:8))


sc_counts_f_int <- sc_counts[, authors_meta_f$barcode]
gene_filter <- apply(sc_counts_f_int, 1, function(row) sum(row > 0) >= 20)
sc_counts_f <- sc_counts_f_int[names(which(gene_filter == TRUE)), ]

rm(sc_counts_f_int)
rm(sc_counts)

table(authors_meta_f$week_stage, authors_meta_f$celltype_region)

seurat_obj <- CreateSeuratObject(counts = sc_counts_f, meta.data = authors_meta_f)
seurat_obj <- NormalizeData(seurat_obj)
Idents(seurat_obj) <- "week_stage"

DE_weeks <- FindAllMarkers(seurat_obj, logfc.threshold = 1) %>% filter(p_val_adj < 0.05)
DE_weeks %>% View()

sc_norm <- seurat_obj@assays$RNA$data

common_genes <- intersect(rownames(vsd), rownames(sc_norm))

scDEGs <- read.csv("scdata/DEGs_week345.csv")
scDEGs <- filter(scDEGs, celltype_region %in% c(
    "Spinal Cord Motor Neuron 1",
    "Spinal Cord Motor Neuron 2",
    "Spinal Cord Motor Neuron 3",
    "Forebrain progenitor",
    "Midbrain progenitor",
    "Hindbrain progenitor",
    "Spinal cord progenitor",
    "Inhibitory Neuron"
) & Gene %in% common_genes)
scDEGs %>% View()
scDEGs$celltype_region %>% table()

genes1 <- filter(scDEGs, celltype_region %in% c("Forebrain progenitor"))$Gene
genes2 <- filter(scDEGs, celltype_region %in% c("Midbrain progenitor"))$Gene

test1 <- sc_norm[genes2, filter(authors_meta_f, week_stage == "W3-1" & celltype_region == "MB")$barcode]
testrank1 <- percent_rank(rowSums(test1))

test2 <- assay(vsd)[genes2, filter(meta, type == "ventral")$sample]
testrank2 <- percent_rank(rowSums(test2))

cor(testrank1, testrank2, method = "pearson")
