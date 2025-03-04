# Loading packages and functions
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(DESeq2)
library(Seurat)
library(reshape2)

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

# loading single cell data from Zeng et al from week3,week4 and week5
sc_counts <- readMM("scdata/week345_wholebody.mtx") %>% t()

# keeping only cells related to Central Nervous System (SNC)
cell_ids <- read.table("scdata/indices_week345_wholebody.csv", header = TRUE, sep = ",")
genes <- read.table("scdata/genes_week345_wholebody.tsv", header = FALSE)$V1
rownames(sc_counts) <- genes

formated_ids <- c(1:nrow(cell_ids)) %>% sapply(function(i) {
    c(str_split(cell_ids[i, ]$index, "-")[[1]][1:2], cell_ids[i, ]$week_stage) %>%
        paste(collapse = "-") %>%
        return()
})
colnames(sc_counts) <- formated_ids

authors_meta <- read.csv("scdata/week345_wholebody_metadata.csv")
authors_meta$celltype_region %>% unique()
authors_meta_f <- filter(
    authors_meta,
    (celltype_region == "FB" & week_stage == "W3-1") |
        (celltype_region == "FB" & week_stage == "W4-1") |
        (celltype_region == "Brain Neuron" & week_stage == "W5-1")
)
table(authors_meta_f$week_stage, authors_meta_f$celltype_region)

# sc_counts_f_int <- sc_counts[, authors_meta_f$barcode]
# gene_filter <- apply(sc_counts_f_int, 1, function(row) sum(row > 0) >= 20)
# sc_counts_f <- sc_counts_f_int[names(which(gene_filter == TRUE)), ]

rm(sc_counts_f_int)

SHH_cluster <- read.csv("results/tables/Figure_4/SHH_cluster.csv", header = TRUE)
SHH_cluster %>% head()
SHH_cluster %>% nrow()

setdiff(SHH_cluster$gene, rownames(sc_counts))
genes <- intersect(rownames(sc_counts), SHH_cluster$gene)

gene_long <- melt(as.matrix(sc_counts[genes, authors_meta_f$barcode]))

rm(sc_counts)

colnames(gene_long) <- c("gene", "barcode", "expression")
gene_long %>% head()
# Merge with week_stage information
merged_data <- merge(gene_long, dplyr::select(authors_meta_f, c("barcode", "week_stage")), by = "barcode")
merged_data
# Calculate percentage of cells with expression > 0 per week_stage and gene
week_gene_df <- merged_data %>%
    group_by(week_stage, gene) %>%
    summarise(percent_expressed = mean(expression > 0) * 100, .groups = "drop") %>%
    tidyr::spread(key = week_stage, value = percent_expressed) %>%
    as.data.frame()

week_gene_df %>% nrow()
write.csv(week_gene_df, "results/tables/Figure_4/sc_percent.csv", row.names = FALSE)
