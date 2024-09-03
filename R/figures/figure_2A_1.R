# Loading packages and functions
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(DESeq2)

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
vsd_blind <- vst(dds, blind = TRUE)

# PCA with all genes
  
genes <- nrow(assay(vsd_blind))
pca.data <- plotPCA.DESeqTransform(vsd_blind, intgroup = c("type", "day", "line"), returnData = TRUE, ntop = genes)
percentVar <- round(100 * attr(pca.data, "percentVar"))
pca_var <- attr(pca.data, "pca_var")

png(filename = "results/images/Figure_2A/F2A_1_PCA_1_2_days_allgenes.png", width = 1600, height = 1200, res = 250)
ggplot(pca.data, aes(PC1, PC2, color = type, shape = day)) +
    geom_point(size = 2, stroke = 1) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_color_manual(values = c("#A1A1DE", "#80AD3C")) +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6)) +
    custom_theme() +
    ggtitle("First and second PCs of dorsal and ventral kinetics all genes")
dev.off()

png(filename = "results/images/Figure_2A/F2A_1_PCA_1_2_line_allgenes.png", width = 1600, height = 1200, res = 250)
ggplot(pca.data, aes(PC1, PC2, color = type, shape = line)) +
    geom_point(size = 2, stroke = 1) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_color_manual(values = c("#A1A1DE", "#80AD3C")) +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6)) +
    custom_theme() +
    ggtitle("First and second PCs of dorsal and ventral kinetics all genes")
dev.off()

png(filename = "results/images/Figure_2A/F2A_1_PCA_2_3_allgenes.png", width = 1600, height = 1200, res = 250)
ggplot(pca.data, aes(PC2, PC3, color = type, shape = line)) +
    geom_point(size = 2, stroke = 1) +
    xlab(paste0("PC2: ", percentVar[2], "% variance")) +
    ylab(paste0("PC3: ", percentVar[3], "% variance")) +
    scale_color_manual(values = c("#A1A1DE", "#80AD3C")) +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6)) +
    custom_theme() +
    ggtitle("Second and third PCs of dorsal and ventral kinetics all genes")
dev.off()

# Variance explained by each PCs
png(filename = "results/images/Figure_2A/F2A_1_percentvar_allgenes.png", width = 1600, height = 1200, res = 250)
ggplot(data.frame(perc = percentVar, PC = factor(colnames(pca.data[1:20]), levels = colnames(pca.data[1:20]))), aes(x = PC, y = perc)) +
    geom_bar(stat = "identity") +
    custom_theme(diag_text = TRUE) +
    ylim(0, 100) +
    ggtitle("Variation explained by each PCs with all genes")
dev.off()

# Building matrix with first 5 PC and covariates
PC_covariate_all_genes <- cbind(pca.data[, 1:5], meta %>%
    dplyr::select(c("line", "type", "day")) %>%
    apply(2, function(x) {
        return(as.numeric(factor(x)) - 1)
    }) %>%
    as.matrix())

# Computing PC-covariate correlation and ANOVA
PC_covariate_all_genes_cor <- cor(PC_covariate_all_genes[, 1:5], PC_covariate_all_genes[, 6:ncol(PC_covariate_all_genes)]) %>% abs()
PC_covariate_all_genes_ANOVA <- c(6:ncol(PC_covariate_all_genes)) %>% lapply(function(i) {
    apply(PC_covariate_all_genes[, 1:5], 2, function(x) {
        aov(x ~ PC_covariate_all_genes[, i])
    }) %>% sapply(function(x) {
        summary(x)[[1]]$`Pr(>F)`[1]
    })
})
PC_covariate_all_genes_ANOVA <- Reduce(cbind, PC_covariate_all_genes_ANOVA)
colnames(PC_covariate_all_genes_ANOVA) <- colnames(PC_covariate_all_genes)[6:ncol(PC_covariate_all_genes)]

# Saving ANOVA results
write.csv(PC_covariate_all_genes_ANOVA, "results/tables/Figure_2A/F2_PC_covariate_ANOVA_all_genes.csv")

rownames(PC_covariate_all_genes_cor) <- paste0(rownames(PC_covariate_all_genes_cor), " (", percentVar[1:5], "%)")

png(filename = "results/images/Figure_2A/F2A_PC_covariate_correlation_all_genes.png", width = 2000, height = 1800, res = 250)
Heatmap(
    PC_covariate_all_genes_cor,
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.2f", PC_covariate_all_genes_cor[i, j]), x, y, gp = gpar(fontsize = 10, fontface = "bold", col = "#646464"))
    },
    name = "Absolute pearson correlation",
    row_title_gp = gpar(fontsize = 20, fontface = "bold"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    column_names_side = "top",
    column_names_rot = 0,
    column_names_centered = TRUE,
    show_column_names = TRUE,
    show_heatmap_legend = TRUE,
    width = ncol(PC_covariate_all_genes_cor) * unit(1.5, "cm"),
    height = nrow(PC_covariate_all_genes_cor) * unit(1, "cm"),
    col = colorRampPalette(c(
        "lightblue",
        "darkblue"
    ))(1000),
)
dev.off()

# PCA with top 500 variable genes
pca.data500 <- plotPCA.DESeqTransform(vsd_blind, intgroup = c("type", "day", "line"), returnData = TRUE, ntop = 500)
percentVar500 <- round(100 * attr(pca.data500, "percentVar"))

png(filename = "results/images/Figure_2A/F2A_1_PCA_1_2_days_500genes.png", width = 1600, height = 1200, res = 250)
ggplot(pca.data500, aes(PC1, PC2, color = type, shape = day)) +
    geom_point(size = 2, stroke = 1) +
    xlab(paste0("PC1: ", percentVar500[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar500[2], "% variance")) +
    scale_color_manual(values = c("#A1A1DE", "#80AD3C")) +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6)) +
    custom_theme() +
    ggtitle("First and second PCs of dorsal and ventral kinetics top 500 variable genes")
dev.off()

png(filename = "results/images/Figure_2A/F2A_1_PCA_1_2_line_500genes.png", width = 1600, height = 1200, res = 250)
ggplot(pca.data500, aes(PC1, PC2, color = type, shape = line)) +
    geom_point(size = 2, stroke = 1) +
    xlab(paste0("PC1: ", percentVar500[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar500[2], "% variance")) +
    scale_color_manual(values = c("#A1A1DE", "#80AD3C")) +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6)) +
    custom_theme() +
    ggtitle("First and second PCs of dorsal and ventral kinetics top 500 variable genes")
dev.off()

png(filename = "results/images/Figure_2A/F2A_1_PCA_2_3_500genes.png", width = 1600, height = 1200, res = 250)
ggplot(pca.data500, aes(PC2, PC3, color = type, shape = line)) +
    geom_point(size = 2, stroke = 1) +
    xlab(paste0("PC2: ", percentVar500[2], "% variance")) +
    ylab(paste0("PC3: ", percentVar500[3], "% variance")) +
    scale_color_manual(values = c("#A1A1DE", "#80AD3C")) +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6)) +
    custom_theme() +
    ggtitle("Second and third PCs of dorsal and ventral kinetics top 500 variable genes")
dev.off()

# Variance explained by each PCs
png(filename = "results/images/Figure_2A/F2A_1_percentvar_500genes.png", width = 1600, height = 1200, res = 250)
ggplot(data.frame(perc = percentVar500, PC = factor(colnames(pca.data500[1:20]), levels = colnames(pca.data500[1:20]))), aes(x = PC, y = perc)) +
    geom_bar(stat = "identity") +
    custom_theme(diag_text = TRUE) +
    ylim(0, 100) +
    ggtitle("Variation explained by each PCs with top 500 variable genes")
dev.off()

# Building dataframe with first 5 PC and covariates
PC_covariate_500 <- cbind(pca.data500[, 1:5], meta %>%
    dplyr::select(c("line", "type", "day")) %>%
    apply(2, function(x) {
        return(as.numeric(factor(x)) - 1)
    }) %>%
    as.matrix())

# Computing PC-covariate correlation and ANOVA
PC_covariate_500_cor <- cor(PC_covariate_500[, 1:5], PC_covariate_500[, 6:ncol(PC_covariate_500)]) %>% abs()
PC_covariate_500_ANOVA <- c(6:ncol(PC_covariate_500)) %>% lapply(function(i) {
    apply(PC_covariate_500[, 1:5], 2, function(x) {
        aov(x ~ PC_covariate_500[, i])
    }) %>% sapply(function(x) {
        summary(x)[[1]]$`Pr(>F)`[1]
    })
})
PC_covariate_500_ANOVA <- Reduce(cbind, PC_covariate_500_ANOVA)
colnames(PC_covariate_500_ANOVA) <- colnames(PC_covariate_500)[6:ncol(PC_covariate_500)]

# Saving ANOVA results
write.csv(PC_covariate_500_ANOVA, "results/tables/Figure_2A/F2_PC_covariate_ANOVA_500.csv")

rownames(PC_covariate_500_cor) <- paste0(rownames(PC_covariate_500_cor), " (", percentVar500[1:5], "%)")

png(filename = "results/images/Figure_2A/F2A_PC_covariate_correlation_500.png", width = 2000, height = 1800, res = 250)
Heatmap(
    PC_covariate_500_cor,
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.2f", PC_covariate_500_cor[i, j]), x, y, gp = gpar(fontsize = 10, fontface = "bold", col = "#646464"))
    },
    name = "Absolute pearson correlation",
    row_title_gp = gpar(fontsize = 20, fontface = "bold"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    column_names_side = "top",
    column_names_rot = 0,
    column_names_centered = TRUE,
    show_column_names = TRUE,
    show_heatmap_legend = TRUE,
    width = ncol(PC_covariate_500_cor) * unit(1.5, "cm"),
    height = nrow(PC_covariate_500_cor) * unit(1, "cm"),
    col = colorRampPalette(c(
        "lightblue",
        "darkblue"
    ))(1000),
)
dev.off()

# PC3 is associated with lineage difference, we want genes higlhy correlated with PC1 and PC2 but not PC3
top_PC1 <- cor(pca.data$PC1, t(assay(vsd_blind))) %>% as.vector()
names(top_PC1) <- rownames(vsd_blind)
top_PC1 <- sort(abs(top_PC1), decreasing = TRUE)[1:1000] %>% names()

top_PC2 <- cor(pca.data$PC2, t(assay(vsd_blind))) %>% as.vector()
names(top_PC2) <- rownames(vsd_blind)
top_PC2 <- sort(abs(top_PC2), decreasing = TRUE)[1:1000] %>% names()

top_PC3 <- cor(pca.data$PC3, t(assay(vsd_blind))) %>% as.vector()
names(top_PC3) <- rownames(vsd_blind)
top_PC3 <- sort(abs(top_PC3), decreasing = TRUE)[1:1000] %>% names()

# get Heatmap genes
hm_genes <- setdiff(union(top_PC1, top_PC2), top_PC3)

# Normalize by variance stabilizing transformation with covariates
vsd <- vst(dds, blind = FALSE)

# making matrix for heatmap, clustering and GO enrichment
vsd_hm <- assay(vsd)[hm_genes, ]

scaled_mat <- t(apply(vsd_hm, 1, scale))
colnames(scaled_mat) <- colnames(vsd_hm)

# get a fix sample order for the kinetic heatmap
sample_order <- c(
    filter(meta, type == "ventral")$sample[order(filter(meta, type == "ventral")$day, decreasing = TRUE)],
    filter(meta, type == "dorsal")$sample[order(filter(meta, type == "dorsal")$day)]
)

# hierarchical clustering using euclidian distance and "complete" method
clustering <- hclust(dist(scaled_mat[, sample_order]))
clusters <- cutree(clustering, k = 4)

# Subclustering within each cluster
sub_clusters_list <- unique(clusters) %>% lapply(function(cluster) {
    sub_mat <- scaled_mat[names(clusters[which(clusters == cluster)]), sample_order]
    sub_clustering <- hclust(dist(sub_mat))
    return(cutree(sub_clustering, k = 4))
})
names(sub_clusters_list) <- paste0("cluster_", unique(clusters))
sub_clusters <- sub_clusters_list %>%
    unname() %>%
    unlist()
sub_clusters <- sub_clusters[names(clusters)]

# cluster and subcluster annotation
clusters_ha <- rowAnnotation(
    cluster = as.character(clusters[clustering$order]),
    sub_cluster = as.character(sub_clusters[clustering$order]),
    col = list(
        cluster = c(
            "1" = "red",
            "2" = "blue",
            "3" = "green",
            "4" = "purple"
        ),
        sub_cluster = c(
            "1" = "black",
            "2" = "pink",
            "3" = "yellow",
            "4" = "brown"
        )
    )
)

#
#
# # Looks like there is no enriched GO terms for the clusters
#
#
# GO enrichment
for (cluster in unique(clusters)) {
    GO_enrichment <- clusterProfiler::enrichGO(names(clusters[which(clusters == cluster)]),
        OrgDb = "org.Hs.eg.db",
        keyType = "ENSEMBL",
        ont = "BP"
    )
    print(GO_enrichment)
    if (is.null(GO_enrichment) || nrow(GO_enrichment) == 0) {
        print(paste0("No enriched terms found for cluster ", cluster))
        next # Skip to the next cluster
    }
    write.csv(GO_enrichment, paste0("results/tables/Figure_2A/GO_enrichment_cluster_", cluster, ".csv"))
    goplot <- clusterProfiler::dotplot(GO_enrichment,
        title = paste0("GO enrichment on cluster", cluster, " (biological processes only)"),
        showCategory = 15
    )
    ggsave(paste0("results/images/Figure_2A/F2A_DE_GO_clust", cluster, ".png"), goplot, width = 8, height = 10)
}

# heatmap
png(filename = "results/images/Figure_2A/F2A_DE_HM.png", width = 2400, height = 1600, res = 250)
Heatmap(
    scaled_mat[clustering$order, sample_order],
    name = "Normalized expression",
    column_names_gp = gpar(fontsize = 6),
    cluster_rows = FALSE,
    left_annotation = clusters_ha,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    row_names_side = "left",
    show_column_names = TRUE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    width = ncol(scaled_mat) * unit(2, "mm"),
    # height = nrow(mat) * unit(5, "mm"),
    col = colorRampPalette(c(
        "black",
        "purple",
        "orange",
        "yellow"
    ))(1000),
)
dev.off()

# Making heatmap for LON71 lineage only
sample_order_LON <- c(
    filter(meta, type == "ventral" & line == "LON71")$sample[order(filter(meta, type == "ventral" & line == "LON71")$day, decreasing = TRUE)],
    filter(meta, type == "dorsal" & line == "LON71")$sample[order(filter(meta, type == "dorsal" & line == "LON71")$day)]
)

clustering_LON <- hclust(dist(scaled_mat[, sample_order_LON]))
clusters_LON <- cutree(clustering_LON, k = 4)

sub_clusters_LON_list <- unique(clusters_LON) %>% lapply(function(cluster) {
    sub_mat <- scaled_mat[names(clusters_LON[which(clusters_LON == cluster)]), sample_order_LON]
    sub_clustering_LON <- hclust(dist(sub_mat))
    return(cutree(sub_clustering_LON, k = 4))
})
sub_clusters_LON <- sub_clusters_LON_list %>%
    unname() %>%
    unlist()
sub_clusters_LON <- sub_clusters_LON[names(clusters_LON)]

clusters_LON_ha <- rowAnnotation(
    cluster = as.character(clusters_LON[clustering_LON$order]),
    sub_cluster = as.character(sub_clusters_LON[clustering_LON$order]),
    col = list(
        cluster = c(
            "1" = "red",
            "2" = "blue",
            "3" = "green",
            "4" = "purple"
        ),
        sub_cluster = c(
            "1" = "black",
            "2" = "pink",
            "3" = "yellow",
            "4" = "brown"
        )
    )
)

# GO enrichment for LON71 lineage only
for (cluster in unique(clusters_LON)) {
    GO_enrichment <- clusterProfiler::enrichGO(names(clusters_LON[which(clusters_LON == cluster)]),
        OrgDb = "org.Hs.eg.db",
        keyType = "ENSEMBL",
        ont = "BP"
    )
    print(GO_enrichment)
    if (is.null(GO_enrichment) || nrow(GO_enrichment) == 0) {
        print(paste0("No enriched terms found for cluster ", cluster))
        next # Skip to the next cluster
    }
    write.csv(GO_enrichment, paste0("results/tables/Figure_2A/GO_enrichment_cluster_", cluster, "_LON.csv"))
    goplot <- clusterProfiler::dotplot(GO_enrichment,
        title = paste0("GO enrichment on cluster", cluster, " (biological processes only)"),
        showCategory = 15
    )
    ggsave(paste0("results/images/Figure_2A/F2A_DE_GO_clust", cluster, "_LON.png"), goplot, width = 8, height = 10)
}

png(filename = "results/images/Figure_2A/F2A_DE_HM_LON.png", width = 2400, height = 1600, res = 250)
Heatmap(
    scaled_mat[clustering_LON$order, sample_order_LON],
    name = "Normalized expression",
    column_names_gp = gpar(fontsize = 6),
    cluster_rows = FALSE,
    left_annotation = clusters_LON_ha,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    row_names_side = "left",
    show_column_names = TRUE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    width = ncol(scaled_mat[, sample_order_LON]) * unit(2, "mm"),
    # height = nrow(mat) * unit(5, "mm"),
    col = colorRampPalette(c(
        "black",
        "purple",
        "orange",
        "yellow"
    ))(1000),
)
dev.off()

#  Making heatmap for WTC lineage only
sample_order_WTC <- c(
    filter(meta, type == "ventral" & line == "WTC")$sample[order(filter(meta, type == "ventral" & line == "WTC")$day, decreasing = TRUE)],
    filter(meta, type == "dorsal" & line == "WTC")$sample[order(filter(meta, type == "dorsal" & line == "WTC")$day)]
)

clustering_WTC <- hclust(dist(scaled_mat[, sample_order_WTC]))
clusters_WTC <- cutree(clustering_WTC, k = 4)

sub_clusters_WTC_list <- unique(clusters_WTC) %>% lapply(function(cluster) {
    sub_mat <- scaled_mat[names(clusters_WTC[which(clusters_WTC == cluster)]), sample_order_WTC]
    sub_clustering_WTC <- hclust(dist(sub_mat))
    return(cutree(sub_clustering_WTC, k = 4))
})
names(sub_clusters_WTC_list) <- paste0("cluster_", unique(clusters_WTC))
sub_clusters_WTC <- sub_clusters_WTC_list %>%
    unname() %>%
    unlist()
sub_clusters_WTC <- sub_clusters_WTC[names(clusters_WTC)]

clusters_WTC_ha <- rowAnnotation(
    cluster = as.character(clusters_WTC[clustering_WTC$order]),
    sub_cluster = as.character(sub_clusters_WTC[clustering_WTC$order]),
    col = list(
        cluster = c(
            "1" = "red",
            "2" = "blue",
            "3" = "green",
            "4" = "purple"
        ),
        sub_cluster = c(
            "1" = "black",
            "2" = "pink",
            "3" = "yellow",
            "4" = "brown"
        )
    )
)

# GO enrichment for WTC lineage only
for (cluster in unique(clusters_WTC)) {
    GO_enrichment <- clusterProfiler::enrichGO(names(clusters_WTC[which(clusters_WTC == cluster)]),
        OrgDb = "org.Hs.eg.db",
        keyType = "ENSEMBL",
        ont = "BP"
    )
    print(GO_enrichment)
    if (is.null(GO_enrichment) || nrow(GO_enrichment) == 0) {
        print(paste0("No enriched terms found for cluster ", cluster))
        next # Skip to the next cluster
    }
    write.csv(GO_enrichment, paste0("results/tables/Figure_2A/GO_enrichment_cluster_", cluster, "_WTC.csv"))
    goplot <- clusterProfiler::dotplot(GO_enrichment,
        title = paste0("GO enrichment on cluster", cluster, " (biological processes only)"),
        showCategory = 15
    )
    ggsave(paste0("results/images/Figure_2A/F2A_DE_GO_clust", cluster, "_WTC.png"), goplot, width = 8, height = 10)
}

png(filename = "results/images/Figure_2A/F2A_DE_HM_WTC.png", width = 2400, height = 1600, res = 250)
Heatmap(
    scaled_mat[clustering_WTC$order, sample_order_WTC],
    name = "Normalized expression",
    column_names_gp = gpar(fontsize = 6),
    cluster_rows = FALSE,
    left_annotation = clusters_WTC_ha,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    row_names_side = "left",
    show_column_names = TRUE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    width = ncol(scaled_mat[, sample_order_WTC]) * unit(2, "mm"),
    # height = nrow(mat) * unit(5, "mm"),
    col = colorRampPalette(c(
        "black",
        "purple",
        "orange",
        "yellow"
    ))(1000),
)
dev.off()

# making table of genes used in heatmap with cluster and subcluster annotation
genes_cluster <- data.frame(
    genes = hm_genes %>% gene_converter("ENSEMBL", "SYMBOL"),
    cluster = clusters[hm_genes],
    subcluster = sub_clusters[hm_genes],
    cluster_LON = clusters_LON[hm_genes],
    subclusterLON = sub_clusters_LON[hm_genes],
    cluster_WTC = clusters_WTC[hm_genes],
    subcluster_WTC = sub_clusters_WTC[hm_genes]
)
rownames(genes_cluster) <- hm_genes

write.csv(genes_cluster, "results/tables/Figure_2A/genes_cluster.csv")
