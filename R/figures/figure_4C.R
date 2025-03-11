# Loading packages and functions
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(DESeq2)
library(Seurat)
library(colorblindr)
library(ggpubr)

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
cellids <- read.table("scdata/indices_week345_wholebody.csv", header = TRUE, sep = ",")
genes <- read.table("scdata/genes_week345_wholebody.tsv", header = FALSE)$V1

# set rownames
rownames(sc_counts) <- genes

# format and set colnames
formated_indices <- c(1:nrow(cellids)) %>% sapply(function(i) {
    c(str_split(cellids[i, ]$index, "-")[[1]][1:2], cellids[i, ]$week_stage) %>%
        paste(collapse = "-") %>%
        return()
})
colnames(sc_counts) <- formated_indices

sc_meta <- read.csv("scdata/week345_wholebody_metadata.csv") %>%
    filter(celltype_region_num %in% c(1:8)) %>%
    dplyr::select(c("week_stage", "celltype_region_num", "celltype_region", "barcode"))

# Creating Seurat object for week4 data, normalizing, scaling, running PCA and UMAP
sc_meta_week4 <- filter(sc_meta, week_stage == "W4-1")
seurat_week4 <- CreateSeuratObject(counts = sc_counts[, sc_meta_week4$barcode], meta.data = sc_meta_week4, min.cells = 20)
seurat_week4 <- NormalizeData(seurat_week4)

seurat_week4 <- FindVariableFeatures(seurat_week4, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_week4)
seurat_week4 <- ScaleData(seurat_week4, features = all.genes)
seurat_week4 <- RunPCA(seurat_week4, features = VariableFeatures(object = seurat_week4))
seurat_week4 <- RunUMAP(seurat_week4, dims = 1:10)

umap_week4 <- seurat_week4[["umap"]]@cell.embeddings %>% as.data.frame()
sc_meta_week4 <- cbind(umap_week4[sc_meta_week4$barcode, ], sc_meta_week4)

################
################
# FIGURE S5B : Zeng et al. single cell data week 4 CNS cells UMAP colored by celltype

# Plotting UMAP for week4 data colored by celltype
umap_plot <- ggplot(data = sc_meta_week4, aes(x = umap_1, y = umap_2, color = celltype_region)) +
    geom_point(size = 2.5) +
    guides(color = guide_legend(
        override.aes = list(size = 5)
    )) +
    labs(color = "Cell identity") + # Modify the legend title # Legend point size
    theme(
        legend.title = element_text(size = 25), # Legend title size
        legend.text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)
    ) +
    scale_color_OkabeIto() +
    custom_theme()
ggsave(paste0("results/images/Figure_S5/Zeng_week4_celltypes.png"), plot = umap_plot, width = 12, height = 10, dpi = 300)

################
################
# FIGURE 4E : Co-expression of NKX2-1 and EMX2 genes with mice tested genes in week 4 CNS cells

# Expression by group of genes
# Define marker gene sets
marker_genes <- c("NKX2-1", "EMX2", "PAX6", "SFTA3", "CAPN6", "EPHB1", "AFF2", "SFRP1", "SALL1", "SHROOM3", "NUAK2", "CNTNAP2", "ZIC5")
NKX2_1_genes <- c("NKX2-1", "SFTA3", "CAPN6", "EPHB1", "AFF2", "SFRP1", "SALL1") #
EMX2_genes <- c("EMX2", "SHROOM3", "NUAK2", "CNTNAP2", "ZIC5")

# Filter Forebrain cells
w4_cellmatch <- sc_meta_week4

# Build dataframe to check expression of marker genes
w4df <- marker_genes %>%
    lapply(function(gene) {
        if (gene %in% rownames(seurat_week4@assays$RNA$counts)) {
            return(seurat_week4@assays$RNA$counts[gene, ] >= 1)
        } else {
            print(paste0(gene, " not in week4"))
            return(rep(FALSE, ncol(seurat_week4@assays$RNA$data)))
        }
    }) %>%
    do.call(cbind, .)
colnames(w4df) <- marker_genes
w4_cellmatch <- cbind(w4_cellmatch, w4df[w4_cellmatch$barcode, ])

# Identify cells expressing NKX2-1 and at least one other gene from NKX2_1_genes
w4_cellmatch$NKX2_1_all <- w4_cellmatch$"NKX2-1" & rowSums(w4_cellmatch[, NKX2_1_genes]) > 1

# Identify cells expressing EMX2 and at least one other gene from EMX2_genes
w4_cellmatch$EMX2_all <- w4_cellmatch$"EMX2" & rowSums(w4_cellmatch[, EMX2_genes]) > 1

# Define UMAP for NKX2-1 genes
title1 <- "Week4: Expression of NKX2-1 Gene Set"
w4_cellmatch$expression_NKX2_1 <- factor(ifelse(w4_cellmatch$NKX2_1_all, "NKX2-1 genes expressed", "None"), levels = c("None", "NKX2-1 genes expressed"))
colors1 <- c("None" = "grey", "NKX2-1 genes expressed" = "red")
plot1 <- ggplot(w4_cellmatch, aes(x = umap_1, y = umap_2)) +
    geom_point(aes(color = expression_NKX2_1), size = 4, alpha = 0.5) +
    geom_point(data = subset(w4_cellmatch, NKX2_1_all), aes(x = umap_1, y = umap_2), color = "red", size = 4) +
    ggtitle(title1) +
    scale_color_manual(values = colors1) +
    custom_theme() +
    theme(
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        panel.background = element_blank(), # Remove background
        plot.background = element_blank(), # Transparent background
        panel.grid = element_blank(), # Remove grid
        axis.line = element_blank(), # Remove axis lines
        axis.text = element_blank(), # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        axis.title = element_blank() # Remove axis labels
    )
ggsave("results/images/Figure_4/NKX2-1_expression.png", plot = plot1, width = 16, height = 12, dpi = 300)

# Define UMAP for EMX2 genes
title2 <- "Week4: Expression of EMX2 Gene Set"
w4_cellmatch$expression_EMX2 <- factor(ifelse(w4_cellmatch$EMX2_all, "EMX2 genes expressed", "None"), levels = c("None", "EMX2 genes expressed"))
colors2 <- c("None" = "grey", "EMX2 genes expressed" = "red")
plot2 <- ggplot(w4_cellmatch, aes(x = umap_1, y = umap_2)) +
    geom_point(aes(color = expression_EMX2), size = 4, alpha = 0.5) +
    geom_point(data = subset(w4_cellmatch, EMX2_all), aes(x = umap_1, y = umap_2), color = "red", size = 4) +
    ggtitle(title2) +
    scale_color_manual(values = colors2) +
    custom_theme() +
    theme(
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        panel.background = element_blank(), # Remove background
        plot.background = element_blank(), # Transparent background
        panel.grid = element_blank(), # Remove grid
        axis.line = element_blank(), # Remove axis lines
        axis.text = element_blank(), # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        axis.title = element_blank() # Remove axis labels
    )

ggsave("results/images/Figure_4/EMX2_expression.png", plot = plot2, width = 16, height = 12, dpi = 300)


################
################
# FIGURE S5C : Expresion in Forebrain progenitors at week 4 of NKX2-1 and EMX2 genes alongside
# mices tested genes

# Filter FB cells
FB_cells <- filter(sc_meta_week4, celltype_region == "FB", umap_1 <= -2, umap_2 <= -1)

sc_meta_week4$celltype_region %>% unique()

# Extract counts for FB cells
FB_counts <- sc_counts[, FB_cells$barcode]

# Check EMX2 expression
FB_cells$EMX2_expr <- FB_counts["EMX2", ] >= 1

# List of EMX2-related genes (excluding EMX2)
emx2_related_genes <- c("SHROOM3", "NUAK2", "CNTNAP2", "ZIC5")

# Define colors

# Loop through each EMX2-related gene
for (gene in emx2_related_genes) {
    color_palette <- c("#ffa20e", "#0fb9fc", "#198733", "gray")
    names(color_palette) <- c("Both", "Only EMX2", paste0("Only ", gene), "None")
    # Check expression of the gene
    FB_cells$gene_expr <- FB_counts[gene, ] >= 1

    # Classify expression
    FB_cells$expression_category <- case_when(
        FB_cells$gene_expr & !FB_cells$EMX2_expr ~ paste0("Only ", gene),
        !FB_cells$gene_expr & FB_cells$EMX2_expr ~ "Only EMX2",
        FB_cells$gene_expr & FB_cells$EMX2_expr ~ "Both",
        TRUE ~ "None"
    )

    # Create UMAP plot
    umap_plot <- ggplot(FB_cells, aes(x = umap_1, y = umap_2, color = expression_category)) +
        geom_point(size = 3) +
        scale_color_manual(values = color_palette) +
        guides(color = guide_legend(override.aes = list(size = 5))) +
        labs(color = "Expression", title = paste("EMX2 and", gene, "Expression")) +
        theme(
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            axis.text.x = element_text(size = 15),
            axis.text.y = element_text(size = 15)
        ) +
        custom_theme()
    # Save plot
    ggsave(paste0("results/images/Figure_S5/Zeng_week4_FB_EMX2__", gene, "_expression.png"),
        plot = umap_plot, width = 12, height = 10, dpi = 300
    )
}

# Check NKX2-1 expression
FB_cells$NKX2_1_expr <- FB_counts["NKX2-1", ] >= 1

# List of NKX2-1-related genes (excluding NKX2-1)
NKX2_1_related_genes <- c("SFTA3", "CAPN6", "EPHB1", "AFF2", "SFRP1", "SALL1")

# Loop through each NKX2-1-related gene
for (gene in NKX2_1_related_genes) {
    color_palette <- c("#ffa20e", "#0fb9fc", "#198733", "gray")
    names(color_palette) <- c("Both", "Only NKX2-1", paste0("Only ", gene), "None")
    # Check expression of the gene
    FB_cells$gene_expr <- FB_counts[gene, ] >= 1

    # Classify expression
    FB_cells$expression_category <- case_when(
        FB_cells$gene_expr & !FB_cells$NKX2_1_expr ~ paste0("Only ", gene),
        !FB_cells$gene_expr & FB_cells$NKX2_1_expr ~ "Only NKX2-1",
        FB_cells$gene_expr & FB_cells$NKX2_1_expr ~ "Both",
        TRUE ~ "None"
    )

    # Create UMAP plot
    umap_plot <- ggplot(FB_cells, aes(x = umap_1, y = umap_2, color = expression_category)) +
        geom_point(size = 3) +
        scale_color_manual(values = color_palette) +
        guides(color = guide_legend(override.aes = list(size = 5))) +
        labs(color = "Expression", title = paste("NKX2-1 and", gene, "Expression")) +
        theme(
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            axis.text.x = element_text(size = 15),
            axis.text.y = element_text(size = 15)
        ) +
        custom_theme()
    # Save plot
    ggsave(paste0("results/images/Figure_S5/Zeng_week4_FB_NKX2-1__", gene, "_expression.png"),
        plot = umap_plot, width = 12, height = 10, dpi = 300
    )
}
