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

# loading single cell data from Zeng et al from week3,week4 and week5
sc_counts <- readMM("scdata/week345_wholebody.mtx") %>% t()
authors_meta <- read.csv("scdata/week345_wholebody_metadata.csv")

cell_ids <- read.table("scdata/indices_week345_wholebody.csv", header = TRUE, sep = ",")
genes <- read.table("scdata/genes_week345_wholebody.tsv", header = FALSE)$V1
rownames(sc_counts) <- genes

# Formating cell ids for convenience
formated_ids <- c(1:nrow(cell_ids)) %>% sapply(function(i) {
    c(str_split(cell_ids[i, ]$index, "-")[[1]][1:2], cell_ids[i, ]$week_stage) %>%
        paste(collapse = "-") %>%
        return()
})
colnames(sc_counts) <- formated_ids

# keeping only cells tagged as Forebrain progenitor (FB) and Brain Neuron
authors_meta_f <- filter(authors_meta, celltype_region_num %in% c(1:8))

# keeping cells in the filtered metadata
sc_counts_f_int <- sc_counts[, authors_meta_f$barcode]
# keeping genes with at least 20 cells expressing it
gene_filter <- apply(sc_counts_f_int, 1, function(row) sum(row > 0) >= 20)
sc_counts_f <- sc_counts_f_int[names(which(gene_filter == TRUE)), ]

rm(sc_counts_f_int)
rm(sc_counts)

# SPON1 and SMIM32 not expressed at all
genes1 <- c("FOXA2", "PTCH1", "SIX3", "SHH", "GSC", "LRP2", "FREM1", "CHRD")
genes2 <- c("NKX2-1", "FGF10", "SLIT2", "DDC", "NOG")
genes3 <- c("NTNG1", "PITX2", "TMEFF2", "CLSTN2", "NEDD9", "SIM1", "KCND3", "LRRK2")
genes <- c(genes1, genes2, genes3) %>% rev()

# Keeping only cells related to FB and Brain Neuron
fb_bn_meta <- filter(authors_meta_f, celltype_region %in% c("FB", "Brain Neuron"))

# Saving cell counts for each celltypes and weeks
table(fb_bn_meta$week_stage, fb_bn_meta$celltype_region) %>% write.csv("results/tables/Figure2F_cellcounts.csv")

# Normalize counts to CPM, filter for genes and build dataframe
sc_norm <- sc_counts_f[, fb_bn_meta$barcode] / colSums(sc_counts_f[, fb_bn_meta$barcode]) * 1e6
sc_df <- sc_norm[genes, ] %>%
    t() %>%
    as.data.frame()
sc_df$week <- fb_bn_meta$week_stage
sc_df$celltype <- fb_bn_meta$celltype_region

# turn dataframe into long format
sc_df_long <- sc_df %>%
    pivot_longer(cols = -c(week, celltype), names_to = "gene", values_to = "expression") %>%
    as.data.frame()

#  Plotting
sc_df_long$gene <- factor(sc_df_long$gene, levels = genes)
sc_boxplot <- ggplot(filter(sc_df_long, expression != 0), aes(x = expression, y = gene, fill = gene)) +
    geom_boxplot() +
    custom_theme() +
    xlim(0, 1500) +
    xlab("Expression (CPM)") +
    theme(
        legend.position = "none",
        axis.text.y = element_text(size = 12, face = "bold.italic"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        strip.text = element_text(size = 20, face = "bold"),
    ) +
    facet_grid(~week, scales = "fixed")
png_save(
    plot = sc_boxplot,
    path = "results/images/Figure_2/Figure2F.png",
    width = 2600,
    height = 2900,
    gg = TRUE
)
