# Loading packages and functions
library(Matrix)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(org.Hs.eg.db)
library(ComplexHeatmap)

library(tibble)

rstudioapi::getSourceEditorContext()$path %>%
    str_split("/") %>%
    unlist() %>%
    head(-3) %>%
    str_c(collapse = "/") %>%
    str_c("/") %>%
    setwd()

source("R/custom_fct.R")
# Loading data (path to change later)
rawcounts <- readcounts("/home/jules/Documents/phd/Data/Article_veranika/bulk/counts.csv", sep = ",", header = TRUE)
meta <- read.table("/home/jules/Documents/phd/Data/Article_veranika/bulk/metadata.csv", sep = ",", header = TRUE)
# Keeping only necessary samples
meta <- filter(meta, type %in% c("dorsal", "ventral"), CRISPR %in% c("control"))

counts <- rawcounts[, meta$sample]
counts <- counts[which(rowSums(counts) > 50), ]
norm <- varianceStabilizingTransformation(counts)
norm <- limma::removeBatchEffect(norm, meta$line)

# Compute dorso_ventral DEGs
DEGs_DV <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ line + type
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "dorsal", "ventral")) %>%
    as.data.frame() %>%
    na.omit()

DEGs_DV_f <- filter(DEGs_DV, padj < 0.01, abs(log2FoldChange) > 1)
DEGs_DV_f$gene <- gene_converter(rownames(DEGs_DV_f), "ENSEMBL", "SYMBOL")
DEGs_DV_f <- filter(DEGs_DV_f, !is.na(gene))


norm_DEGs <- norm[rownames(DEGs_DV_f), ]
scaled_mat <- t(apply(norm_DEGs, 1, scale))
colnames(scaled_mat) <- colnames(norm_DEGs)

clustering <- hclust(dist(scaled_mat))
clusters_1 <- cutree(clustering, k = 2)
clusters_2 <- cutree(clustering, k = 5)

clusters_ha <- rowAnnotation(
    cluster_1 = as.character(clusters_1[clustering$order]),
    cluster_2 = as.character(clusters_2[clustering$order]),
    col = list(
        cluster_1 = c(
            "1" = "red",
            "2" = "blue"
        ),
        cluster_2 = c(
            "1" = "red",
            "2" = "blue",
            "3" = "green",
            "4" = "purple",
            "5" = "orange"
        )
    )
)

for (cluster in unique(clusters_1)) {
    print(cluster)
    GO_enrichment <- clusterProfiler::enrichGO(names(clusters_1[which(clusters_1 == cluster)]),
        OrgDb = "org.Hs.eg.db",
        keyType = "ENSEMBL",
        ont = "BP"
    )
    write.csv(GO_enrichment, paste0("results/tables/Figure_1/GO_enrichment_cluster_", cluster, ".csv"))
    goplot <- clusterProfiler::dotplot(GO_enrichment,
        title = paste0("GO enrichment on cluster", cluster, " (biological processes only)"),
        showCategory = 15
    )
    ggsave(paste0("results/images/Figure_1/F1_DE_GO_clust", cluster, ".png"), goplot, width = 8, height = 10)
}

sample_ha <- columnAnnotation(
    line = meta$line,
    type = meta$type,
    col = list(
        line = c("LON" = "#00f7ff", "WTC" = "#ff0000"),
        type = c("dorsal" = "#A1A1DE", "ventral" = "#80AD3C")
    )
)
png(filename = "results/images/Figure_1/F1_3_DE_HM.png", width = 1600, height = 1600, res = 250)
Heatmap(
    scaled_mat[clustering$order, ],
    name = "Normalized expression",
    row_title_gp = gpar(fontsize = 16, fontface = "bold"),
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    left_annotation = clusters_ha,
    bottom_annotation = sample_ha,
    row_names_side = "left",
    show_column_names = TRUE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    width = ncol(scaled_mat) * unit(4, "mm"),
    # height = nrow(mat) * unit(5, "mm"),
    col = colorRampPalette(c(
        "black",
        "purple",
        "orange",
        "yellow"
    ))(1000),
)
dev.off()

DEGs_DV_f$cluster_1 <- clusters_1[rownames(DEGs_DV_f)]
DEGs_DV_f$cluster_2 <- clusters_2[rownames(DEGs_DV_f)]

write.csv(DEGs_DV_f, "results/tables/Figure_1/dorsal_VS_ventral_DEGs.csv")
