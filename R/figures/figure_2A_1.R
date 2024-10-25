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

# PCA with 3000 genes
genes <- 3000
pca.data <- plotPCA.DESeqTransform(vsd_blind, intgroup = c("type", "day", "line"), returnData = TRUE, ntop = genes)
percentVar <- round(100 * attr(pca.data, "percentVar"))
pca_var <- attr(pca.data, "factoextra")

png(filename = "results/images/Figure_2A/F2A_1_PCA_1_2_days_3000genes.png", width = 1600, height = 1200, res = 250)
ggplot(pca.data, aes(PC1, PC2, color = line, shape = day)) +
    geom_point(size = 2, stroke = 1) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_color_manual(values = c("#868686", "#000000")) +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6)) +
    custom_theme() +
    ggtitle("First and second PCs of dorsal and ventral kinetics all genes")
dev.off()

# Variance explained by each PCs
png(filename = "results/images/Figure_2A/F2A_1_percentvar_3000genes.png", width = 1600, height = 1200, res = 250)
ggplot(data.frame(perc = percentVar, PC = factor(colnames(pca.data[1:20]), levels = colnames(pca.data[1:20]))), aes(x = PC, y = perc)) +
    geom_bar(stat = "identity") +
    custom_theme(diag_text = TRUE) +
    ylim(0, 100) +
    ggtitle("Variation explained by each PCs with all genes")
dev.off()

# Building matrix with first 5 PC and covariates
PC_covariate_3000genes <- cbind(pca.data[, 1:5], meta %>%
    dplyr::select(c("line", "type", "day")) %>%
    apply(2, function(x) {
        return(as.numeric(factor(x)) - 1)
    }) %>%
    as.matrix())

# Computing PC-covariate correlation and ANOVA
PC_covariate_3000genes_cor <- cor(PC_covariate_3000genes[, 1:5], PC_covariate_3000genes[, 6:ncol(PC_covariate_3000genes)]) %>% abs()
PC_covariate_3000genes_ANOVA <- c(6:ncol(PC_covariate_3000genes)) %>% lapply(function(i) {
    apply(PC_covariate_3000genes[, 1:5], 2, function(x) {
        aov(x ~ PC_covariate_3000genes[, i])
    }) %>% sapply(function(x) {
        summary(x)[[1]]$`Pr(>F)`[1]
    })
})
PC_covariate_3000genes_ANOVA <- Reduce(cbind, PC_covariate_3000genes_ANOVA)
colnames(PC_covariate_3000genes_ANOVA) <- colnames(PC_covariate_3000genes)[6:ncol(PC_covariate_3000genes)]

# Saving ANOVA results
write.csv(PC_covariate_3000genes_ANOVA, "results/tables/Figure_2A/F2_PC_covariate_ANOVA_3000genes.csv")

rownames(PC_covariate_3000genes_cor) <- paste0(rownames(PC_covariate_3000genes_cor), " (", percentVar[1:5], "%)")

png(filename = "results/images/Figure_2A/F2A_PC_covariate_correlation_3000genes.png", width = 2000, height = 1800, res = 250)
Heatmap(
    PC_covariate_3000genes_cor,
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.2f", PC_covariate_3000genes_cor[i, j]), x, y, gp = gpar(fontsize = 10, fontface = "bold", col = "#646464"))
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
    width = ncol(PC_covariate_3000genes_cor) * unit(1.5, "cm"),
    height = nrow(PC_covariate_3000genes_cor) * unit(1, "cm"),
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
hm_genes
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
    col = list(
        cluster = c(
            "1" = "#b16060",
            "2" = "#4d6da5",
            "3" = "#78588c",
            "4" = "#5e9a5e"
        )
    )
)

sub_clusters_ha <- rowAnnotation(
    sub_cluster = as.character(sub_clusters[clustering$order]),
    col = list(
        sub_cluster = c(
            "1" = "black",
            "2" = "pink",
            "3" = "yellow",
            "4" = "brown"
        )
    )
)

# GO enrichment
for (cluster in unique(clusters)) {
    print(cluster)
    GO_enrichment <- clusterProfiler::enrichGO(names(clusters[which(clusters == cluster)]),
        OrgDb = "org.Hs.eg.db",
        keyType = "ENSEMBL",
        ont = "BP"
    )
    GO_results <- GO_enrichment@result
    GO_results$GeneRatio <- sapply(GO_enrichment@result$GeneRatio, function(x) {
        eval(parse(text = x))
    }) %>% unname()
    GO_results_f <- GO_results[order(GO_results$GeneRatio, decreasing = TRUE)[1:15], ]
    GO_results_f$Description <- str_wrap(GO_results_f$Description, width = 40)
    GO_results_f$Description <- factor(GO_results_f$Description, levels = rev(GO_results_f$Description))
    goplot <- ggplot(GO_results_f, aes(x = GeneRatio, y = Description, fill = p.adjust)) +
        geom_bar(stat = "identity") +
        custom_theme() +
        scale_fill_gradient(low = "#e06663", high = "#327eba") +
        ggtitle(paste0("GO enrichment on cluster", cluster, " (biological processes only)"))
    write.csv(GO_enrichment, paste0("results/tables/Figure_2A/GO_enrichment_cluster_", cluster, ".csv"))
    ggsave(paste0("results/images/Figure_2A/F2A_DE_GO_clust", cluster, ".png"), goplot, width = 15, height = 10)
}

# heatmap
png(filename = "results/images/Figure_2A/F2A_DE_HM.png", width = 2400, height = 1600, res = 250)
Heatmap(
    scaled_mat[clustering$order, sample_order],
    name = "Normalized expression",
    column_names_gp = gpar(fontsize = 6),
    cluster_rows = FALSE,
    left_annotation = sub_clusters_ha,
    right_annotation = clusters_ha,
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
    col = list(
        cluster = c(
            "1" = "#b16060",
            "2" = "#4d6da5",
            "3" = "#78588c",
            "4" = "#5e9a5e"
        )
    )
)

sub_clusters_LON_ha <- rowAnnotation(
    sub_cluster = as.character(sub_clusters_LON[clustering_LON$order]),
    col = list(
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
    print(cluster)
    GO_enrichment <- clusterProfiler::enrichGO(names(clusters[which(clusters == cluster)]),
        OrgDb = "org.Hs.eg.db",
        keyType = "ENSEMBL",
        ont = "BP"
    )
    GO_results <- GO_enrichment@result
    GO_results$GeneRatio <- sapply(GO_enrichment@result$GeneRatio, function(x) {
        eval(parse(text = x))
    }) %>% unname()
    GO_results_f <- GO_results[order(GO_results$GeneRatio, decreasing = TRUE)[1:15], ]
    GO_results_f$Description <- factor(GO_results_f$Description, levels = rev(GO_results_f$Description))
    goplot <- ggplot(GO_results_f, aes(x = GeneRatio, y = Description, fill = p.adjust)) +
        geom_bar(stat = "identity") +
        custom_theme() +
        scale_fill_gradient(low = "#e06663", high = "#327eba") +
        ggtitle(paste0("GO enrichment on cluster", cluster, " (biological processes only)"))
    write.csv(GO_enrichment, paste0("results/tables/Figure_2A/GO_enrichment_cluster_", cluster, "_LON.csv"))
    ggsave(paste0("results/images/Figure_2A/F2A_DE_GO_clust", cluster, "_LON.png"), goplot, width = 20, height = 10)
}

png(filename = "results/images/Figure_2A/F2A_DE_HM_LON.png", width = 2400, height = 1600, res = 250)
Heatmap(
    scaled_mat[clustering_LON$order, sample_order_LON],
    name = "Normalized expression",
    column_names_gp = gpar(fontsize = 6),
    cluster_rows = FALSE,
    left_annotation = sub_clusters_LON_ha,
    right_annotation = clusters_LON_ha,
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
    col = list(
        cluster = c(
            "1" = "#b16060",
            "2" = "#4d6da5",
            "3" = "#78588c",
            "4" = "#5e9a5e"
        )
    )
)

sub_clusters_WTC_ha <- rowAnnotation(
    sub_cluster = as.character(sub_clusters_WTC[clustering_WTC$order]),
    col = list(
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
    GO_results <- GO_enrichment@result
    GO_results$GeneRatio <- sapply(GO_enrichment@result$GeneRatio, function(x) {
        eval(parse(text = x))
    }) %>% unname()
    GO_results_f <- GO_results[order(GO_results$GeneRatio, decreasing = TRUE)[1:5], ]
    GO_results_f$Description <- str_wrap(GO_results_f$Description, width = 60)
    GO_results_f$Description <- factor(GO_results_f$Description, levels = rev(GO_results_f$Description))
    goplot <- ggplot(GO_results_f, aes(x = GeneRatio, y = reorder(Description, GeneRatio), fill = p.adjust)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = Description),
            hjust = 1.1, # Move text inside the bar, adjust as needed
            color = "black", # Make the text white for better visibility
            size = 10
        ) + # Adjust size to fit the text inside the bar
        custom_theme() +
        scale_fill_gradient(low = "#e06663", high = "#327eba") +
        ggtitle(paste0("GO enrichment on cluster", cluster, " (biological processes only)")) +
        theme(
            axis.text.y = element_blank(), # Remove y-axis text
            axis.ticks.y = element_blank(),
            legend.text = element_text(size = 12), # Adjusts the legend text size
            legend.title = element_text(size = 14), # Adjusts the legend title size
            legend.key.size = unit(1.5, "lines")
        )
    write.csv(GO_enrichment, paste0("results/tables/Figure_2A/GO_enrichment_cluster_diapo", cluster, "_WTC.csv"))
    ggsave(paste0("results/images/Figure_2A/F2A_DE_GO_clust_diapo", cluster, "_WTC.png"), goplot, width = 20, height = 10)
}

png(filename = "results/images/Figure_2A/F2A_DE_HM_WTC.png", width = 2400, height = 1600, res = 250)
Heatmap(
    scaled_mat[clustering_WTC$order, sample_order_WTC],
    name = "Normalized expression",
    column_names_gp = gpar(fontsize = 6),
    cluster_rows = FALSE,
    left_annotation = sub_clusters_WTC_ha,
    right_annotation = clusters_WTC_ha,
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






rawcounts <- readcounts("data/rawcounts.csv", sep = ",", header = TRUE)
rawmeta <- read.table("data/meta.csv", sep = ",", header = TRUE)

# LON71_D12_2 does not have any reads in the count file
# though, the fastQC report shows that the sample is good
lp_meta <- filter(rawmeta, (sample != "LON71_D12_2" & diff == "diff13" & line %in% c("LON71", "WTC")) | (sample == "WTC6cipc"))
lp_meta[which(lp_meta$sample == "WTC6cipc"), "day"] <- "day00"
View(lp_meta)
# filtering out lowly expressed genes
lp_counts <- rawcounts[, lp_meta$sample][which(rowSums(rawcounts[, lp_meta$sample]) >= 25), ]
lp_meta$day
lp_meta$type
lp_meta$line
# making DESeq object with lineage,days and type as covariates
lp_vsd <- DESeqDataSetFromMatrix(
    countData = lp_counts,
    colData = lp_meta,
    design = ~ line + day
) %>% vst(blind = FALSE)


# genes_ventral <- c("NKX2-1", "HTR1D", "PTCH1", "FGF19", "FOXA2", "FREM1", "LINC00261", "CAPN6")
genes_ventral1 <- c("FOXA2", "FREM1", "SHH", "NKX2-1", "PTCH1", "LINC00261", "PLCL1", "CAPN6", "LRRK2", "DDC", "SMIM32", "GADL1", "DRC1", "CLSTN2", "SLIT2")
genes_ventral2 <- c("LRRK2", "DDC", "SMIM32", "GADL1", "DRC1", "CLSTN2", "SLIT2")
genes_dorsal1 <- c("CNTN2", "PAX6", "PAX3", "CNTNAP2", "GLI3", "EMX2", "NELL2", "GDF7", "SYT4", "GRIP2")
genes_dorsal2 <- c("EMX2", "NELL2", "GDF7", "SYT4", "GRIP2")

vAN_meta <- filter(lp_meta, type != "dorsal")
vAN_vsd <- assay(lp_vsd)[, vAN_meta$sample]
scaled_vAN_vsd <- vAN_vsd - min(vAN_vsd)

vAN_df_1 <- t(scaled_vAN_vsd[genes_ventral1 %>% gene_converter("SYMBOL", "ENSEMBL"), ])
colnames(vAN_df_1) <- colnames(vAN_df_1) %>% gene_converter("ENSEMBL", "SYMBOL")
vAN_df_1 <- cbind(dplyr::select(vAN_meta, "day"), vAN_df_1)

head(vAN_df_1)

df_long_1 <- vAN_df_1 %>%
    pivot_longer(
        cols = -day, # All columns except 'day'
        names_to = "gene", # New column for gene names
        values_to = "expression"
    ) # New column for expression values

# Calculate mean and standard error for each gene on each day
df_summary_1 <- df_long_1 %>%
    group_by(day, gene) %>%
    summarise(
        expression_mean = mean(expression, na.rm = TRUE),
        expression_se = sd(expression, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
    )
library(paletteer)
ggplot(df_summary_1, aes(x = day, y = expression_mean, group = gene, color = gene)) +
    geom_point(size = 3) +
    geom_line(size = 2) +
    geom_errorbar(aes(ymin = expression_mean - expression_se, ymax = expression_mean + expression_se), width = 0.1) +
    scale_color_paletteer_d("ggsci::category20_d3") + # 20 distinct colors
    ylim(0, 7.5) +
    labs(
        x = "Days",
        y = "Mean scaled normalized expression"
    ) +
    custom_theme()
ggsave("/home/jules/Documents/phd/projects/panasenkava_2024/results/images/Figure_2A/F2A_lineplots_genes_ventral.png", width = 12, height = 10)


dAN_meta <- filter(lp_meta, type != "ventral")
dAN_vsd <- assay(lp_vsd)[, dAN_meta$sample]
scaled_dAN_vsd <- dAN_vsd - min(dAN_vsd)

dAN_df_1 <- t(scaled_dAN_vsd[genes_dorsal1 %>% gene_converter("SYMBOL", "ENSEMBL"), ])
colnames(dAN_df_1) <- colnames(dAN_df_1) %>% gene_converter("ENSEMBL", "SYMBOL")
dAN_df_1 <- cbind(dplyr::select(dAN_meta, "day"), dAN_df_1)

head(dAN_df_1)

df_long_1 <- dAN_df_1 %>%
    pivot_longer(
        cols = -day, # All columns except 'day'
        names_to = "gene", # New column for gene names
        values_to = "expression"
    ) # New column for expression values

# Calculate mean and standard error for each gene on each day
df_summary_1 <- df_long_1 %>%
    group_by(day, gene) %>%
    summarise(
        expression_mean = mean(expression, na.rm = TRUE),
        expression_se = sd(expression, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
    )

ggplot(df_summary_1, aes(x = day, y = expression_mean, group = gene, color = gene)) +
    geom_point(size = 3) +
    geom_line(size = 2) +
    geom_errorbar(aes(ymin = expression_mean - expression_se, ymax = expression_mean + expression_se), width = 0.1) +
    scale_color_paletteer_d("ggsci::category20_d3") + # 20 distinct colors
    ylim(0, 7.5) +
    labs(
        x = "Days",
        y = "Mean scaled normalized expression"
    ) +
    custom_theme()
ggsave("/home/jules/Documents/phd/projects/panasenkava_2024/results/images/Figure_2A/F2A_lineplots_genes_dorsal.png", width = 12, height = 10)
