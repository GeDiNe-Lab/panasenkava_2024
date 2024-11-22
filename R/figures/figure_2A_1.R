# Loading packages and functions
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(ggrepel)
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
meta <- filter(rawmeta, sample != "LON71_D12_2", diff == "diff13", line %in% c("LON71", "WTC"), ((manip == "veranika" & day != "day12") | (manip == "lauryane" & day == "day12")))
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
    geom_point(size = 3, stroke = 1.5) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_color_manual(values = c("#868686", "#000000")) +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6)) +
    custom_theme()
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

# Replace values
clusters_LON <- ifelse(clusters_LON == 1, "Early differentiation",
    ifelse(clusters_LON == 2, "Anterior forebrain",
        ifelse(clusters_LON == 3, "Ventral forebrain", "Dorsal forebrain")
    )
)
clusters_LON
# cluster and subcluster annotation
clusters_LON_ha <- rowAnnotation(
    cluster = as.character(clusters_LON[clustering_LON$order]),
    col = list(
        cluster = c(
            "Early differentiation" = "#b16060",
            "Anterior forebrain" = "#4d6da5",
            "Ventral forebrain" = "#5e9a5e",
            "Dorsal forebrain" = "#78588c"
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
    GO_enrichment <- clusterProfiler::enrichGO(names(clusters_LON[which(clusters_LON == cluster)]),
        OrgDb = "org.Hs.eg.db",
        keyType = "ENSEMBL",
        ont = "BP"
    )
    GO_results <- GO_enrichment@result
    GO_results$GeneRatio <- sapply(GO_enrichment@result$GeneRatio, function(x) {
        eval(parse(text = x))
    }) %>% unname()
    GO_results_f <- GO_results[order(GO_results$GeneRatio, decreasing = TRUE)[1:5], ]

    GO_results_f$Description <- str_wrap(GO_results_f$Description, width = 42) %>% str_to_upper()
    GO_results_f$Description <- factor(GO_results_f$Description, levels = rev(GO_results_f$Description))

    goplot <- ggplot(GO_results_f, aes(x = GeneRatio, y = reorder(Description, GeneRatio), fill = p.adjust)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = Description),
            hjust = 1.01, # Move text inside the bar, adjust as needed
            color = "black", # Make the text white for better visibility
            size = 13
        ) + # Adjust size to fit the text inside the bar
        custom_theme() +
        scale_fill_gradient(name = "p-value", low = "#e06663", high = "#327eba") +
        theme(
            axis.title.x = element_text(size = 30), # Adjusts the x-axis title size
            axis.text.x = element_text(size = 20),
            axis.text.y = element_blank(), # Remove y-axis text
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            legend.text = element_text(size = 20), # Adjusts the legend text size
            legend.title = element_text(size = 30), # Adjusts the legend title size
            legend.key.size = unit(2, "lines")
        )
    goplot
    write.csv(GO_enrichment, paste0("results/tables/Figure_2A/GO_enrichment_cluster_", cluster, "_LON.csv"))
    ggsave(paste0("results/images/Figure_2A/F2A_DE_GO_clust", cluster, "_LON.png"), goplot, width = 19, height = 10)
}

png(filename = "results/images/Figure_2A/F2A_DE_HM_LON.png", width = 2400, height = 1600, res = 20)
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
clusters_WTC <- cutree(clustering_WTC, h = 9)

# Find the heights for merging clusters
merge_heights <- sort(clustering_WTC$height, decreasing = TRUE)

# Interval for h to get two clusters
h_start <- merge_heights[4] # Second largest height where split happens
h_end <- merge_heights[3] # Largest height before all merge into one cluster

# Output interval
cat("Interval for h to get k=2: (", h_start, ",", h_end, "]\n")



# Replace values
clusters_WTC <- ifelse(clusters_WTC == 1, "Early differentiation",
    ifelse(clusters_WTC == 2, "Anterior forebrain",
        ifelse(clusters_WTC == 3, "Ventral forebrain", "Dorsal forebrain")
    )
)

# cluster and subcluster annotation
clusters_WTC_ha <- rowAnnotation(
    cluster = as.character(clusters_WTC[clustering_WTC$order]),
    col = list(
        cluster = c(
            "Early differentiation" = "#b16060",
            "Anterior forebrain" = "#4d6da5",
            "Ventral forebrain" = "#5e9a5e",
            "Dorsal forebrain" = "#78588c"
        )
    )
)

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

    GO_results_f$Description <- str_wrap(GO_results_f$Description, width = 42) %>% str_to_upper()
    GO_results_f$Description <- factor(GO_results_f$Description, levels = rev(GO_results_f$Description))

    goplot <- ggplot(GO_results_f, aes(x = GeneRatio, y = reorder(Description, GeneRatio), fill = p.adjust)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = Description),
            hjust = 1.01, # Move text inside the bar, adjust as needed
            color = "black", # Make the text white for better visibility
            size = 13
        ) + # Adjust size to fit the text inside the bar
        custom_theme() +
        scale_fill_gradient(name = "p-value", low = "#e06663", high = "#327eba") +
        theme(
            axis.title.x = element_text(size = 30), # Adjusts the x-axis title size
            axis.text.x = element_text(size = 20),
            axis.text.y = element_blank(), # Remove y-axis text
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            legend.text = element_text(size = 20), # Adjusts the legend text size
            legend.title = element_text(size = 30), # Adjusts the legend title size
            legend.key.size = unit(2, "lines")
        )
    write.csv(GO_enrichment, paste0("results/tables/Figure_2A/GO_enrichment_cluster_diapo", cluster, "_WTC.csv"))
    ggsave(paste0("results/images/Figure_2A/F2A_DE_GO_clust", cluster, "_WTC.png"), goplot, width = 19, height = 10)
}


png(filename = "results/images/Figure_2A/F2A_DE_HM_WTC_test.png", width = 2400, height = 1600, res = 260)
Heatmap(
    scaled_mat[clustering_WTC$order, sample_order_WTC],
    name = "Normalized expression",
    column_names_gp = gpar(fontsize = 6),
    cluster_rows = FALSE,
    left_annotation = sub_clusters_WTC_ha,
    right_annotation = clusters_WTC_ha,
    cluster_columns = TRUE,
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


library(paletteer)
library(DEGreport)
rawcounts <- readcounts("data/rawcounts.csv", sep = ",", header = TRUE)
rawmeta <- read.table("data/meta.csv", sep = ",", header = TRUE)

dbd_ventral <- read.csv("/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2A/Volcano_DEG_dbd_ventral.csv", header = TRUE)
genes <- dbd_ventral$gene

lp_meta <- filter(rawmeta, (sample != "LON71_D12_2" & diff == "diff13" & line == "WTC" & type == "ventral" & ((manip == "veranika" & day != "day12") | (manip == "lauryane" & day == "day12"))))
# filtering out lowly expressed genes
lp_counts <- rawcounts[, lp_meta$sample][which(rowSums(rawcounts[, lp_meta$sample]) >= 25), ]
lp_meta %>% View()
# making DESeq object with lineage,days and type as covariates
lp_vsd <- DESeqDataSetFromMatrix(
    countData = lp_counts,
    colData = lp_meta,
    design = ~day
) %>% vst(blind = FALSE)
filtered <- assay(lp_vsd)[rownames(lp_vsd) %in% gene_converter(genes, "SYMBOL", "ENSEMBL"), ]
rownames(lp_meta) <- lp_meta$sample
lp_meta$day <- as.factor(lp_meta$day)

clusters <- degPatterns(
    filtered,
    meta = lp_meta,
    time = "day",
    reduce = TRUE,
    nClusters = 10,
)
df_cluster <- clusters$df
clusters$normalize %>% View()


df_cluster$symbol <- df_cluster$genes %>% gene_converter("ENSEMBL", "SYMBOL")
df_cluster$refseq <- df_cluster$genes %>% gene_converter("ENSEMBL", "REFSEQ")
df_cluster %>% View()
filter(df_cluster, cluster == 9)$symbol

write.csv(df_cluster, file = "/home/jules/Documents/phd/projects/panasenkava_2024/results/images/Figure_2A/DEGpattern.csv", quote = FALSE, row.names = FALSE)
assay(lp_vsd)[filter(df_cluster, cluster == 9)$genes, c("WTC6c_V12_L1", "WTC6c_V12_L2")] %>% View()

degPlotCluster(table, time,
    color = NULL, min_genes = 10,
    process = FALSE, points = TRUE, boxes = TRUE, smooth = TRUE,
    lines = TRUE, facet = TRUE, cluster_column = "cluster"
)


gene_converter(c("CAPN6", "FREM1", "FGF19", "HTR1D", "SPON1", "GPC3", "LRP2", "SOX3"), "SYMBOL", "REFSEQ")











# LON71_D12_2 does not have any reads in the count file
# though, the fastQC report shows that the sample is good
lp_meta <- filter(rawmeta, (sample != "LON71_D12_2" & diff == "diff13" & line == "WTC" & ((manip == "veranika" & day != "day12") | (manip == "lauryane" & day == "day12"))) | (sample == "WTC6cipc"))
lp_meta[which(lp_meta$sample == "WTC6cipc"), "day"] <- "day00"
# filtering out lowly expressed genes
lp_counts <- rawcounts[, lp_meta$sample][which(rowSums(rawcounts[, lp_meta$sample]) >= 25), ]

# making DESeq object with lineage,days and type as covariates
lp_vsd <- DESeqDataSetFromMatrix(
    countData = lp_counts,
    colData = lp_meta,
    design = ~day
) %>% vst(blind = FALSE)

dbd_ventral <- read.csv("/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2A/Volcano_DEG_dbd_ventral.csv", header = TRUE)
dbd_ventral <- filter(dbd_ventral, gene %in% gene_converter(rownames(lp_vsd), "ENSEMBL", "SYMBOL"))

dbd_dorsal <- read.csv("/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2A/Volcano_DEG_dbd_dorsal.csv", header = TRUE)
dbd_dorsal <- filter(dbd_dorsal, gene %in% gene_converter(rownames(lp_vsd), "ENSEMBL", "SYMBOL"))

selection_v_known <- c("FGF10", "LRP2", "FGF9", "FOXA2", "NKX2-2", "RAX", "NKX2-1", "SIX3", "SHH", "PTCH1", "GADL1", "DRC1", "CLSTN2")
selection_v_new <- c("CAPN6", "PLCL1", "FRZB", "FREM1", "LINC01833", "LINC00261", "SPON1", "DDC", "SLIT2", "LRRK2", "SMIM32")

selection_d <- c("CNTN2", "CNTNAP2", "EMX2", "GDF7", "GLI3", "GRIP2", "NELL2", "PAX3", "PAX6", "SYT4")

loess_marker_genes <- c("NKX2-1", "PTCH1", "FOXA2", "FGF19", "FREM1")


mean_df <- assay(lp_vsd)[!rownames(assay(lp_vsd)) %in% c(loess_marker_genes %>% gene_converter("SYMBOL", "ENSEMBL")), ]

DESeq2::DEpa

DE <- "UP_then_DOWN"

# by_day
# genes_ventral <- filter(
#     dbd_ventral,
#     day_04_02 == "UP",
#     # day_06_04 == "UP",
#     # day_08_06 == "UP",
#     # day_10_08 == "UP",
#     # day_12_10 != "UP"
# )$gene
# genes_ventral
# # early
genes_ventral <- filter(
    dbd_ventral,
    (day_04_02 == "UP" | day_06_04 == "UP"),
    (day_04_02 != "DOWN" & day_06_04 != "DOWN"),
    day_08_06 != "UP",
    day_10_08 != "UP",
    day_12_10 != "UP"
)$gene
# genes_ventral
# tardif
# genes_ventral <- filter(
#     dbd_ventral,
#     day_04_02 != "UP",
#     day_06_04 != "UP",
#     (day_08_06 == "UP" | day_10_08 == "UP" | day_12_10 == "UP"),
#     (day_08_06 != "DOWN" & day_10_08 != "DOWN" & day_12_10 != "DOWN")
# )$gene

genes_ventral <- loess_marker_genes

lp_meta <- filter(lp_meta, type != "dorsal")
# scaled_lp_vsd <- assay(lp_vsd)[, lp_meta$sample] - min(assay(lp_vsd)[, lp_meta$sample])
lp_df_1 <- t(assay(lp_vsd)[genes_ventral %>% gene_converter("SYMBOL", "ENSEMBL"), lp_meta$sample])
colnames(lp_df_1) <- colnames(lp_df_1) %>% gene_converter("ENSEMBL", "SYMBOL")
lp_df_1 <- cbind(dplyr::select(lp_meta, "day"), lp_df_1)
lp_df_1

# Reshape dataframe to long format and compute mean by day for each gene
df_long <- lp_df_1 %>%
    pivot_longer(
        cols = -day, # All columns except "day"
        names_to = "gene",
        values_to = "value"
    ) %>%
    group_by(day, gene) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") # Compute mean


df_long$selection <- ifelse(df_long$gene %in% selection_v_known, "known", "no")
df_long$selection <- ifelse(df_long$gene %in% selection_v_new, "new", df_long$selection)

df_long$day_numeric <- match(df_long$day, c("day00", "day02", "day04", "day06", "day08", "day10", "day12"))

clustering <- rep(0, length(colnames(dplyr::select(lp_df_1, !c("day")))))
# clustering <- hclust(as.dist(1 - cor(dplyr::select(lp_df_1, !c("day"))))) %>%
#     dynamicTreeCut::cutreeDynamic(cutHeight = 0.9, minClusterSize = 45)
clustering %>% table()
names(clustering) <- colnames(dplyr::select(lp_df_1, !c("day")))
df_long$clustering <- clustering[df_long$gene]

df_fit <- data.frame(day = unique(df_long$day))
df_fit <- cbind(
    df_fit, Reduce(cbind, lapply(unique(clustering), function(cluster) {
        return(unique(predict(loess(mean_value ~ day_numeric, data = filter(df_long, clustering == cluster), span = 0.25))))
    }))
)
colnames(df_fit) <- c("day", paste0("fit", (as.character(unique(clustering)))))


loess_fit <- df_fit



label_positions <- df_long %>%
    filter(selection != "no") %>%
    group_by(gene) %>%
    slice_tail(n = 1)

for (cluster in unique(clustering)) {
    fit <- colnames(df_fit)[2 + cluster]
    df_fit$fit <- df_fit[, fit]
    lineplot <- ggplot() +
        # Plot "Not selected" genes in light grey
        geom_line(
            data = filter(df_long, clustering == cluster, selection == "no"),
            aes(x = day, y = mean_value, group = gene),
            color = "lightgrey", size = 0.5
        ) +
        # Plot "Selected" genes with unique colors using the color palette
        geom_line(
            data = filter(df_long, clustering == cluster, selection != "no"),
            aes(x = day, y = mean_value, group = gene, color = selection),
            size = 1
        ) +
        # Add labels for "Selected" genes at the last point
        geom_text_repel(
            data = filter(label_positions, clustering == cluster),
            aes(x = day, y = mean_value, label = gene, color = selection),
            size = 5, nudge_x = 1, direction = "y",
            hjust = 0, segment.color = "grey50", segment.size = 0.5, force = 0.1
        ) +
        # Overlay the overall mean in black with a thicker line
        geom_line(
            data = df_fit, aes(x = day, y = fit, group = 1),
            color = "black", size = 1.5 # Thicker line for the mean
        ) +
        # Add a label for the mean line
        annotate("text",
            x = max(df_fit$day), y = max(df_fit$fit),
            label = "Mean Expression", color = "black", size = 5, hjust = 1, vjust = 0
        ) +
        # Customize labels and theme
        labs(
            y = "Scaled mean expression by day",
        ) +
        # Apply color palette for selected genes
        scale_color_paletteer_d("ggsci::category20_d3") + # Use the "category20_d3" palette
        custom_theme(diag_text = TRUE) +
        theme(
            legend.position = "none", # Remove legend
            axis.title.x = element_blank(),
            axis.text.x = element_text(size = 16),
        )

    ggsave(plot = lineplot, paste0("/home/jules/Documents/phd/projects/panasenkava_2024/results/images/Figure_2A/F2A_lineplots_WTC_genes_ventral_", DE, cluster, ".png"), width = 12, height = 10)
}







rawcounts <- readcounts("data/rawcounts.csv", sep = ",", header = TRUE)
rawmeta <- read.table("data/meta.csv", sep = ",", header = TRUE)

# LON71_D12_2 does not have any reads in the count file
# though, the fastQC report shows that the sample is good
lp_meta <- filter(rawmeta, (sample != "LON71_D12_2" & diff == "diff13" & line == "WTC" & ((manip == "veranika" & day != "day12") | (manip == "lauryane" & day == "day12"))) | (sample == "WTC6cipc"))
lp_meta[which(lp_meta$sample == "WTC6cipc"), "day"] <- "day00"
# filtering out lowly expressed genes
lp_counts <- rawcounts[, lp_meta$sample][which(rowSums(rawcounts[, lp_meta$sample]) >= 25), ]

# making DESeq object with lineage,days and type as covariates
lp_vsd <- DESeqDataSetFromMatrix(
    countData = lp_counts,
    colData = lp_meta,
    design = ~day
) %>% vst(blind = FALSE)

dbd_ventral <- read.csv("/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2A/Volcano_DEG_dbd_ventral.csv", header = TRUE)
dbd_ventral <- filter(dbd_ventral, gene %in% gene_converter(rownames(lp_vsd), "ENSEMBL", "SYMBOL"))

dbd_dorsal <- read.csv("/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2A/Volcano_DEG_dbd_dorsal.csv", header = TRUE)
dbd_dorsal <- filter(dbd_dorsal, gene %in% gene_converter(rownames(lp_vsd), "ENSEMBL", "SYMBOL"))

selection_v_known <- c("FGF10", "LRP2", "FGF9", "FOXA2", "NKX2-2", "RAX", "NKX2-1", "SIX3", "SHH", "PTCH1", "GADL1", "DRC1", "CLSTN2")
selection_v_new <- c("CAPN6", "PLCL1", "FRZB", "FREM1", "LINC01833", "LINC00261", "SPON1", "DDC", "SLIT2", "LRRK2", "SMIM32")
# selection_v_new <- c("")

selection_d <- c("CNTN2", "CNTNAP2", "EMX2", "GDF7", "GLI3", "GRIP2", "NELL2", "PAX3", "PAX6", "SYT4")

excluded <- "dorsal"
dbd <- dbd_ventral

DE_list <- list(
    filter(
        dbd,
        (day_04_02 == "UP" | day_06_04 == "UP"),
        (day_04_02 != "DOWN" & day_06_04 != "DOWN"),
        day_08_06 != "UP",
        day_10_08 != "UP",
        day_12_10 != "UP"
    )$gene,
    filter(dbd, day_04_02 == "UP")$gene,
    filter(dbd, day_04_02 != "UP", day_06_04 == "UP")$gene,
    filter(dbd, day_04_02 != "UP", day_06_04 != "UP", day_08_06 == "UP")$gene,
    filter(dbd, day_04_02 != "UP", day_06_04 != "UP", day_08_06 != "UP", day_10_08 == "UP")$gene,
    filter(dbd, day_04_02 != "UP", day_06_04 != "UP", day_08_06 != "UP", day_10_08 != "UP", day_12_10 == "UP")$gene
)
names(DE_list) <- c("UP_then_DOWN", "from_day_2", "from_day_4", "from_day_6", "from_day_8", "from_day_10")

for (DE in names(DE_list)) {
    filtered_genes <- DE_list[DE] %>%
        unlist() %>%
        unname()
    filtered_genes
    lp_meta <- filter(lp_meta, type != excluded)
    scaled_lp_vsd <- assay(lp_vsd)[, lp_meta$sample] - min(assay(lp_vsd)[, lp_meta$sample])
    lp_df_1 <- t(scaled_lp_vsd[filtered_genes %>% gene_converter("SYMBOL", "ENSEMBL"), ])
    colnames(lp_df_1) <- colnames(lp_df_1) %>% gene_converter("ENSEMBL", "SYMBOL")
    lp_df_1 <- cbind(dplyr::select(lp_meta, "day"), lp_df_1)


    # Reshape dataframe to long format and compute mean by day for each gene
    df_long <- lp_df_1 %>%
        pivot_longer(
            cols = -day, # All columns except "day"
            names_to = "gene",
            values_to = "value"
        ) %>%
        group_by(day, gene) %>%
        summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") # Compute mean



    df_long$selection <- ifelse(df_long$gene %in% selection_v_known, "known", "no")
    df_long$selection <- ifelse(df_long$gene %in% selection_v_new, "new", df_long$selection)
    print(df_long$selection)
    df_long$day_numeric <- match(df_long$day, c("day00", "day02", "day04", "day06", "day08", "day10", "day12"))

    clustering <- rep(0, length(colnames(dplyr::select(lp_df_1, !c("day")))))
    names(clustering) <- colnames(dplyr::select(lp_df_1, !c("day")))
    df_long$clustering <- clustering[df_long$gene]

    df_fit <- data.frame(day = unique(df_long$day))
    df_fit <- cbind(
        df_fit, Reduce(cbind, lapply(unique(clustering), function(cluster) {
            return(unique(predict(loess(mean_value ~ day_numeric, data = filter(df_long, clustering == cluster), span = 0.25))))
        }))
    )
    colnames(df_fit) <- c("day", paste0("fit", (as.character(unique(clustering)))))

    label_positions <- df_long %>%
        filter(selection != "no") %>%
        group_by(gene) %>%
        slice_tail(n = 1)

    for (cluster in unique(clustering)) {
        fit <- colnames(df_fit)[2 + cluster]
        df_fit$fit <- df_fit[, fit]
        lineplot <- ggplot() +
            # Plot "Not selected" genes in light grey
            geom_line(
                data = filter(df_long, clustering == cluster, selection == "no"),
                aes(x = day, y = mean_value, group = gene),
                color = "lightgrey", size = 0.5
            ) +
            # Plot "Selected" genes with unique colors using the color palette
            geom_line(
                data = filter(df_long, clustering == cluster, selection != "no"),
                aes(x = day, y = mean_value, group = gene, color = selection),
                size = 1
            ) +
            scale_color_manual(
                values = c("new" = "#ff7f0e", "known" = "#1f77b4") # Define colors for "new" and "known"
            ) +
            # Add labels for "Selected" genes at the last point
            geom_text_repel(
                data = filter(label_positions, clustering == cluster),
                aes(x = day, y = mean_value, label = gene, color = selection),
                size = 9, nudge_x = 2, direction = "y",
                hjust = 0, segment.color = "grey50", segment.size = 1, force = 2, force.pull = 0.2
            ) +
            # Overlay the overall mean in black with a thicker line
            geom_line(
                data = df_fit, aes(x = day, y = fit, group = 1),
                color = "black", size = 1.5 # Thicker line for the mean
            ) +
            # Customize labels and theme
            labs(
                y = "Scaled mean expression by day",
            ) +
            # Apply color palette for selected genes
            custom_theme(diag_text = TRUE) +
            theme(
                legend.position = "none", # Remove legend
                axis.title.x = element_blank(),
                axis.text.x = element_text(size = 16),
                plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt") # Increase the right margin
            )

        ggsave(plot = lineplot, paste0("/home/jules/Documents/phd/projects/panasenkava_2024/results/images/Figure_2A/specific_timepoints_lineplots/F2A_lineplots_WTC_genes_ventral_", DE, ".png"), width = 12, height = 10)
    }
}
