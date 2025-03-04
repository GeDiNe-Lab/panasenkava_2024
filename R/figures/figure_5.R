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

# shorten some celltypes names
sc_meta$celltype_region[which(sc_meta$celltype_region == "Spinal Cord Motor Neuron 1")] <- "SCMN 1"
sc_meta$celltype_region[which(sc_meta$celltype_region == "Spinal Cord Motor Neuron 2")] <- "SCMN 2"
sc_meta$celltype_region[which(sc_meta$celltype_region == "Spinal Cord Motor Neuron 3")] <- "SCMN 3"

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
ggsave(paste0("results/images/Figure_5/Zeng_week4_celltypes.png"), plot = umap_plot, width = 12, height = 10, dpi = 300)

# Expression by group of genes
# Define marker gene sets
marker_genes <- c("NKX2-1", "EMX2", "PAX6", "SFTA3", "CAPN6", "EPHB1", "AFF2", "SFRP1", "SALL1", "SHROOM3", "NUAK2", "CNTNAP2", "ZIC5")
NKX2_1_genes <- c("NKX2-1", "SFTA3", "CAPN6", "EPHB1", "AFF2", "SFRP1", "SALL1") #
PAX6_genes <- c("PAX6", "SHROOM3", "NUAK2", "CNTNAP2", "ZIC5")
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

# Identify cells expressing PAX6 and at least one other gene from PAX6_genes
w4_cellmatch$PAX6_all <- w4_cellmatch$"PAX6" & rowSums(w4_cellmatch[, PAX6_genes]) > 1

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
plot1
ggsave("results/images/Figure_5/NKX2-1_expression.png", plot = plot1, width = 16, height = 12, dpi = 300)

# Define UMAP for PAX6 genes
title2 <- "Week4: Expression of PAX6 Gene Set"
w4_cellmatch$expression_PAX6 <- factor(ifelse(w4_cellmatch$PAX6_all, "PAX6 genes expressed", "None"), levels = c("None", "PAX6 genes expressed"))
colors2 <- c("None" = "grey", "PAX6 genes expressed" = "red")
plot2 <- ggplot(w4_cellmatch, aes(x = umap_1, y = umap_2)) +
    geom_point(aes(color = expression_PAX6), size = 4, alpha = 0.5) +
    geom_point(data = subset(w4_cellmatch, PAX6_all), aes(x = umap_1, y = umap_2), color = "red", size = 4) +
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

ggsave("results/images/Figure_5/PAX6_expression.png", plot = plot2, width = 16, height = 12, dpi = 300)

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

ggsave("results/images/Figure_5/EMX2_expression.png", plot = plot2, width = 16, height = 12, dpi = 300)


WGCNA <- read.table("results/tables/Figure_4/SHH_cluster.csv", sep = ",", header = TRUE)
WGCNA_genes <- intersect(WGCNA$gene, rownames(seurat_week4@assays$RNA$counts))
WGCNA_genes
WGCNA_pos_genes <- filter(WGCNA, cor > 0 & gene %in% WGCNA_genes)$gene
WGCNA_neg_genes <- c(filter(WGCNA, cor < 0 & gene %in% WGCNA_genes)$gene, "EMX2")

# Filter Forebrain cells
w4_cellmatch <- sc_meta_week4

# Build dataframe to check expression of marker genes
w4df <- c(WGCNA_genes, "EMX2", "EPHB1") %>%
    lapply(function(gene) {
        if (gene %in% rownames(seurat_week4@assays$RNA$counts)) {
            return(seurat_week4@assays$RNA$counts[gene, ] >= 1)
        } else {
            print(paste0(gene, " not in week4"))
            return(rep(FALSE, ncol(seurat_week4@assays$RNA$data)))
        }
    }) %>%
    do.call(cbind, .)
colnames(w4df) <- c(WGCNA_genes, "EMX2", "EPHB1")
w4_cellmatch <- cbind(w4_cellmatch, w4df[w4_cellmatch$barcode, ])


# Get percentage of expression
compute_gene_percentages <- function(gene_list, marker) {
    nkx2_1_cells <- w4_cellmatch[w4_cellmatch[, marker], ]
    percentages <- sapply(gene_list, function(gene) {
        if (gene %in% colnames(w4_cellmatch)) {
            mean(nkx2_1_cells[[gene]]) * 100
        } else {
            NA
        }
    })
    return(data.frame(Gene = gene_list, Percentage = percentages))
}

# Example usage
percentages_df_NKX <- compute_gene_percentages(c(WGCNA_pos_genes, "EPHB1"), "NKX2-1")
percentages_df_EMX2 <- compute_gene_percentages(WGCNA_neg_genes, "EMX2")


ggplot() +
    geom_jitter(data = percentages_df_EMX2, aes(x = "Negative genes/EMX2", y = Percentage), color = "grey") +
    geom_jitter(data = percentages_df_NKX, aes(x = "Positive genes/NKX2-1", y = Percentage), color = "grey") +
    geom_boxplot(data = percentages_df_EMX2, aes(x = "Negative genes/EMX2", y = Percentage), fill = NA) +
    geom_boxplot(data = percentages_df_NKX, aes(x = "Positive genes/NKX2-1", y = Percentage), fill = NA) +
    custom_theme()


write.csv(percentages_df_NKX, "results/tables/Figure_5/NKX2-1_percentages.csv", row.names = FALSE)
write.csv(percentages_df_EMX2, "results/tables/Figure_5/EMX2_percentages.csv", row.names = FALSE)
