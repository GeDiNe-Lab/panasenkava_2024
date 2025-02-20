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
umap_plot
ggsave(paste0("results/images/Figure_5/Zeng_week4_celltypes.png"), plot = umap_plot, width = 12, height = 10, dpi = 300)

# genes selection
marker_genes <- c("SFTA3", "CAPN6", "EPHB1", "AFF2", "SFRP1", "SALL1", "SHROOM3", "NUAK2", "CNTNAP2", "ZIC5")

# Plotting UMAP for week4 data colored by marker genes keeping only Forebrain cells
w4_cellmatch <- filter(sc_meta_week4, celltype_region == "FB" & umap_1 < -2 & umap_2 < -1)
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

w4df <- union(marker_genes, c("NKX2-1", "PAX6")) %>%
    lapply(function(gene) {
        if (gene %in% rownames(seurat_week4@assays$RNA$counts)) {
            return(seurat_week4@assays$RNA$counts[gene, ] > 1)
        } else {
            print(paste0(gene, " not in week4"))
            return(rep(FALSE, ncol(seurat_week4@assays$RNA$data)))
        }
    }) %>%
    do.call(cbind, .)
colnames(w4df) <- union(marker_genes, c("NKX2-1", "PAX6"))

w4_cellmatch <- cbind(w4_cellmatch, w4df[w4_cellmatch$barcode, ])

pairwise_df <- expand.grid(marker_genes, c("NKX2-1", "PAX6"))
pairwise_df$Var1 <- as.vector(pairwise_df$Var1)
pairwise_df$Var2 <- as.vector(pairwise_df$Var2)


# execute this !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
for (i in c(1:nrow(pairwise_df))) {
    gene1 <- pairwise_df$Var1[i]
    gene2 <- pairwise_df$Var2[i]
    w4_cellmatch$gene1 <- w4_cellmatch[, pairwise_df$Var1[i]]
    w4_cellmatch$gene2 <- w4_cellmatch[, pairwise_df$Var2[i]]
    w4_cellmatch$duo <- sapply(c(1:nrow(w4_cellmatch)), function(i) {
        if (w4_cellmatch$gene1[i] == TRUE & w4_cellmatch$gene2[i] == TRUE) {
            return("1-both")
        } else if (w4_cellmatch$gene1[i] == TRUE) {
            return(paste0("2-", gene1))
        } else if (w4_cellmatch$gene2[i] == TRUE) {
            return(paste0("3-", gene2))
        } else {
            return("4-none")
        }
    })
    w4_cellmatch$expression <- w4_cellmatch$duo
    w4_cellmatch <- w4_cellmatch[order(w4_cellmatch$expression, decreasing = TRUE), ]

    # w4_cellmatch[, c("umap1", "umap2", "expression")] %>% head()
    colors <- c(palette_OkabeIto[1], "grey", palette_OkabeIto[2], palette_OkabeIto[4])
    names(colors) <- c("1-both", "4-none", paste0("2-", pairwise_df$Var1[i]), paste0("3-", pairwise_df$Var2[i]))
    plot <- ggplot(data = w4_cellmatch[, c("umap_1", "umap_2", "expression", "gene1", "gene2")], aes(x = umap_1, y = umap_2, color = expression)) +
        geom_point(size = 4) +
        ggtitle(paste0("Week4: gene1 = ", gene1, " and gene2 = ", gene2)) +
        scale_color_manual(values = colors) +
        custom_theme() +
        theme(
            legend.text = element_text(size = 15), # Change legend text size
            legend.title = element_text(size = 20), # Change legend title size
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20)
        )
    ggsave(paste0("results/images/Figure_5/", gene1, "_", gene2, ".png"), plot = plot, width = 8, height = 7, dpi = 300)
}
