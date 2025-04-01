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


#########################
#########################
#  Computing percent of WGCNA modules genes expressed in FB cells

WGCNA_genes <- read.csv("results/tables/Figure_4/SHH_cluster.csv", header = TRUE)

sc_meta_week4_FB <- filter(sc_meta, week_stage == "W4-1", celltype_region == "FB")
seurat_week4_FB <- CreateSeuratObject(counts = sc_counts[, sc_meta_week4_FB$barcode], meta.data = sc_meta_week4_FB)

FB_genes <- which(rowSums(seurat_week4_FB@assays$RNA$counts > 0) >= nrow(sc_meta_week4_FB) * 0.05) %>% names()
FB_WGCNA_genes <- intersect(FB_genes, WGCNA_genes$gene)

seurat_week4@meta.data$FB_DE <- ifelse(seurat_week4@meta.data$celltype_region == "FB", "FB", "Other")
Idents(seurat_week4) <- "FB_DE"
DE <- FindMarkers(seurat_week4, ident.1 = "FB", ident.2 = "Other")

DE_WGCNA_genes <- intersect(rownames(filter(DE, p_val_adj < 0.05, avg_log2FC >= 0)), FB_WGCNA_genes)

WGCNA_genes$sc_forebrain <- ifelse(WGCNA_genes$gene %in% FB_WGCNA_genes, "expressed", NA)
WGCNA_genes$sc_forebrain_DE <- ifelse(WGCNA_genes$gene %in% DE_WGCNA_genes, "DE", NA)

# percent of genes expressed and DE in Forebrain
## Among all genes
### expressed : 56.88 %
nrow(filter(WGCNA_genes, !is.na(sc_forebrain))) / nrow(WGCNA_genes) * 100
### DE 12.05 %
nrow(filter(WGCNA_genes, !is.na(sc_forebrain_DE))) / nrow(WGCNA_genes) * 100
## Among positively correlated genes with SHH
### expressed 48.08 %
nrow(filter(WGCNA_genes, !is.na(sc_forebrain), cor > 0)) / nrow(filter(WGCNA_genes, cor > 0)) * 100
### DE 5.75 %
nrow(filter(WGCNA_genes, !is.na(sc_forebrain_DE), cor > 0)) / nrow(filter(WGCNA_genes, cor > 0)) * 100
## Among negatively correlated genes with SHH
### expressed 65.10 %
nrow(filter(WGCNA_genes, !is.na(sc_forebrain), cor < 0)) / nrow(filter(WGCNA_genes, cor < 0)) * 100
### DE 17.93 %
nrow(filter(WGCNA_genes, !is.na(sc_forebrain_DE), cor < 0)) / nrow(filter(WGCNA_genes, cor < 0)) * 100

write.csv(WGCNA_genes, "results/tables/Figure_4/SHH_cluster.csv", row.names = FALSE)

#########################
#########################
#  Number of cells for each region in week 4
sc_week4_cellcount <- sc_meta_week4$celltype_region %>%
    table() %>%
    as.data.frame()
sc_week4_cellcount <- rbind(sc_week4_cellcount, data.frame("." = c("Total"), "Freq" = c(nrow(sc_meta_week4))))
colnames(sc_week4_cellcount) <- c("Region", "Number of cells")

write.csv(sc_week4_cellcount, "results/tables/Figure_4/week4_cellcount.csv", row.names = FALSE)

################
################
# FIGURE 4E : Expression in Zeng CNS sc data at week 4 of genes selected from Figure 4D

ventral_genes <- c("CAPN6", "EPHB1", "QKI", "SLIT1", "RGMA", "PDZRN3", "FRZB", "SFRP1", "SALL1") #
dorsal_genes <- c("SHROOM3", "NUAK2", "CNTNAP2", "ZIC5")

ctmat_v <- seurat_week4@assays$RNA$counts[ventral_genes, ]
sc_meta_week4$ventral_counts <- ctmat_v %>%
    apply(2, function(x) {
        return(x >= 5)
    }) %>%
    colSums()

ctmat_d <- seurat_week4@assays$RNA$counts[dorsal_genes, ]
sc_meta_week4$dorsal_counts <- ctmat_d %>%
    apply(2, function(x) {
        return(x >= 3)
    }) %>%
    colSums()

df <- sc_meta_week4

df <- df %>% arrange(ventral_counts)
ggplot(df, aes(x = umap_1, y = umap_2, fill = ventral_counts)) +
    geom_point(shape = 21, size = 2.5, color = "#1d1d1d", stroke = 0.1) +
    scale_fill_gradientn(
        colors = c("white", "#e9ffca", "#80AD3C", "#162500"), # 4-color gradient
        name = "Number of genes\n expressed (total 9)",
        guide = guide_colorbar(
            frame.colour = "black", # Adds a black border
            frame.linewidth = 0.5 # Thickness of the border
        )
    ) +
    custom_theme() +
    theme(
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)
    )
ggsave("results/images/Figure_4/ventral_UMAP.png", width = 10, height = 8, dpi = 300)

df <- df %>% arrange(dorsal_counts)
ggplot(df, aes(x = umap_1, y = umap_2, fill = dorsal_counts)) +
    geom_point(shape = 21, size = 2.5, color = "#1d1d1d", stroke = 0.1) +
    scale_fill_gradientn(
        colors = c("white", "#ceceff", "#7676c2", "#191971"), # 4-color gradient
        name = "Number of genes\n expressed (total 4)",
        guide = guide_colorbar(
            frame.colour = "black", # Adds a black border
            frame.linewidth = 0.5 # Thickness of the border
        )
    ) +
    custom_theme() +
    theme(
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)
    )
ggsave("results/images/Figure_4/dorsal_UMAP.png", width = 10, height = 8, dpi = 300)
