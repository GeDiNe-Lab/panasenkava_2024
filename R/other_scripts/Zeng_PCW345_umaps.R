library(Matrix)
library(tidyverse)
library(Seurat)
source("/home/jules/Documents/phd/custom_fct.R")

meta <- read.table("/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/meta_final.csv", header = TRUE, sep = ",")
genes <- read.table("/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/genes_final.csv", header = TRUE, sep = ",")$x
counts <- readMM("/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/counts_w345.mtx")
colnames(counts) <- meta$index
rownames(counts) <- genes

meta <- meta %>% dplyr::filter(index %in% colnames(counts))

SNC_genes <- read.table("/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/SNC_DEGs.csv", sep = ";", header = TRUE)
SNC_genes <- filter(SNC_genes, pvals_adj < 0.05 & logfoldchanges > log(2))

counts_f <- counts[rownames(counts) %in% unique(SNC_genes$Gene), ]

set.seed(42)
kmean_clustering <- kmeans(centers = 2, counts_f %>% t(), iter.max = 10, nstart = 10)

meta$SN_cluster <- kmean_clustering$cluster[meta$index]

wilcox.test(
    rowSums(counts_f[, filter(meta, SN_cluster == 1)$index]),
    rowSums(counts_f[, filter(meta, SN_cluster == 2)$index]),
    alternative = "less"
)

meta_SN <- meta %>% dplyr::filter(SN_cluster == 2)
seurat_obj <- CreateSeuratObject(counts = counts[, meta_SN$index])
rm(counts)
rm(counts_f)
rm(meta)
gc()

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_obj)
seurat_obj@assays$RNA$counts <- NULL
gc()
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
# seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
# seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

umap_res <- seurat_obj[["umap"]]@cell.embeddings %>% as.data.frame()
meta_SN <- cbind(meta_SN, umap_res[meta_SN$index, ])

marker_genes <- read.table("/home/jules/Documents/phd/Data/Article_veranika/single_cell/Zeng_et_al/List4_UMAP.csv")$V1 %>% unique()
marker_genes <- c(marker_genes, c("EGR2"))
setdiff(marker_genes, rownames(seurat_obj@assays$RNA$data))
marker_genes <- intersect(marker_genes, rownames(seurat_obj@assays$RNA$data))
marker_genes
seurat_obj@assays$RNA$data <- seurat_obj@assays$RNA$data[marker_genes, meta_SN$index]
all.genes <- rownames(seurat_obj@assays$RNA$data)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

gc()
for (gene in marker_genes) {
    meta_SN$expression <- seurat_obj@assays$RNA$scale.data[gene, meta_SN$index]
    umap_plot <- ggplot(data = meta_SN, aes(x = umap_1, y = umap_2, color = expression)) +
        geom_point(size = 1) +
        ggtitle(gene) +
        scale_color_gradient(low = "lightblue", high = "darkred") +
        custom_theme()
    ggsave(paste0("/home/jules/Documents/phd/Data/results/genes_umap_Zeng_3/", gene, ".png"), plot = umap_plot, width = 14, height = 10, dpi = 300)
    gc()
}
week_umap_plot <- ggplot(data = meta_SN, aes(x = umap_1, y = umap_2, color = week_stage)) +
    geom_point(size = 1) +
    custom_theme()
ggsave("/home/jules/Documents/phd/Data/results/PCW345_Zeng_weeks.png", plot = week_umap_plot, width = 14, height = 10, dpi = 300)
