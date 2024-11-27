# Loading packages and functions
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(ggrepel)
library(DESeq2)
library(paletteer)
library(DEGreport)
library(ggsignif)

# Setting working directory
rstudioapi::getSourceEditorContext()$path %>%
    str_split("/") %>%
    unlist() %>%
    head(-3) %>%
    str_c(collapse = "/") %>%
    str_c("/") %>%
    setwd()
source("R/custom_fct.R")

# Loading custom functions
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




#  Lineplots :
dbd_ventral <- read.csv("/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2A/Volcano_DEG_dbd_ventral.csv", header = TRUE)
genes <- dbd_ventral$gene

lp_meta <- filter(rawmeta, (sample != "LON71_D12_2" & diff == "diff13" & line == "WTC" & type == "ventral" & ((manip == "veranika" & day != "day12") | (manip == "lauryane" & day == "day12"))))
# filtering out lowly expressed genes
lp_counts <- rawcounts[, lp_meta$sample][which(rowSums(rawcounts[, lp_meta$sample]) >= 25), ]
# lp_counts <- rawcounts[hm_genes, lp_meta$sample]

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
clusters$df$symbol <- clusters$df$genes %>% gene_converter("ENSEMBL", "SYMBOL")
write.csv(clusters$df, file = "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2A/DEGpattern_WTC.csv", quote = FALSE, row.names = FALSE)


ggsave("results/images/Figure_2A/F2A_DEpattern_WTC_full.png", clusters$plot, width = 15, height = 10)

c1 <- filter(clusters$normalize, cluster == 1)
c7 <- filter(clusters$normalize, cluster == 7)

sign_comp <- list(
    c("day02", "day04"),
    c("day04", "day06"),
    c("day06", "day08"),
    c("day08", "day10"),
    c("day10", "day12")
)
source("R/custom_fct.R")

c1_plot <- MyDegPlotCluster(table = c1, time = "day", sign_comp = sign_comp)
ggsave("results/images/Figure_2A/F2A_DEpattern_LON_c1.png", c1_plot, width = 15, height = 10)

c7_plot <- MyDegPlotCluster(table = c7, time = "day", sign_comp = sign_comp)
ggsave("results/images/Figure_2A/F2A_DEpattern_LON_c7.png", c7_plot, width = 15, height = 10)

clusters$df$symbol <- clusters$df$genes %>% gene_converter("ENSEMBL", "SYMBOL")
write.csv(clusters$df, file = "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2A/DEGpattern_LON.csv", quote = FALSE, row.names = FALSE)

early <- c("PTCH1", "FREM1", "FOXA2")
mid <- c("SHH", "SIX3", "LRP2", "FGF9", "FGF10", "NKX2-1", "NKX2-2")
late <- c("GADL1", "DRC1", "CLSTN2")

symbol_vsd <- assay(lp_vsd)
rownames(symbol_vsd) <- gene_converter(rownames(symbol_vsd), "ENSEMBL", "SYMBOL")
filtered_vsd <- symbol_vsd[Reduce(union, list(early, mid, late)), ]
scaled_filtered_vsd <- filtered_vsd - min(filtered_vsd)

early_vsd <- scaled_filtered_vsd[early, ]
early_df <- cbind(dplyr::select(lp_meta, c("day")), t(early_vsd)[lp_meta$sample, ])
early_df_long <- reshape2::melt(early_df, id.vars = "day", variable.name = "gene", value.name = "expression")
early_df_long <- summarise(group_by(early_df_long, gene, day), mean_expression = mean(expression), sd_expression = sd(expression))

early_plot <- ggplot(early_df_long, aes(x = day, y = mean_expression, color = gene, group = gene, label = gene)) +
    geom_point() +
    geom_line(size = 1) +
    geom_errorbar(
        aes(
            ymin = mean_expression - sd_expression,
            ymax = mean_expression + sd_expression
        ),
        width = 0, # Width of the horizontal bar on the error bar
        size = 1 # Thickness of the error bars
    ) +
    geom_text_repel(
        data = filter(early_df_long, day == "day12"),
        aes(x = as.numeric(day) + 0.75),
        nudge_x = 0,
        nudge_y = 0.1,
        direction = "y",
        size = 10,
        segment.color = NA
    ) +
    ylim(0, 3.5) +
    ylab("Scaled normalized expression") +
    custom_theme() +
    theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 20),
        plot.margin = margin(t = 10, r = 30, b = 10, l = 10)
    )
ggsave("results/images/Figure_2A/F2A_early_LON.png", early_plot, width = 15, height = 10)

mid_vsd <- scaled_filtered_vsd[mid, ]
mid_df <- cbind(dplyr::select(lp_meta, c("day")), t(mid_vsd)[lp_meta$sample, ])
mid_df_long <- reshape2::melt(mid_df, id.vars = "day", variable.name = "gene", value.name = "expression")
mid_df_long <- summarise(group_by(mid_df_long, gene, day), mean_expression = mean(expression), sd_expression = sd(expression))

mid_plot <- ggplot(mid_df_long, aes(x = day, y = mean_expression, color = gene, group = gene, label = gene)) +
    geom_point() +
    geom_line(size = 1) +
    geom_errorbar(
        aes(
            ymin = mean_expression - sd_expression,
            ymax = mean_expression + sd_expression
        ),
        width = 0, # Width of the horizontal bar on the error bar
        size = 1 # Thickness of the error bars
    ) +
    geom_text_repel(
        data = filter(mid_df_long, day == "day12"),
        aes(x = as.numeric(day) + 0.75),
        nudge_x = 0,
        nudge_y = 0.1,
        direction = "y",
        size = 10,
        segment.color = NA
    ) +
    ylim(0, 3.5) +
    ylab("Scaled normalized expression") +
    custom_theme() +
    theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 20),
        plot.margin = margin(t = 10, r = 30, b = 10, l = 10) # Add margin on the right
    )
ggsave("results/images/Figure_2A/F2A_mid_LON.png", mid_plot, width = 15, height = 10)

late_vsd <- scaled_filtered_vsd[late, ]
late_df <- cbind(dplyr::select(lp_meta, c("day")), t(late_vsd)[lp_meta$sample, ])
late_df_long <- reshape2::melt(late_df, id.vars = "day", variable.name = "gene", value.name = "expression")
late_df_long <- summarise(group_by(late_df_long, gene, day), mean_expression = mean(expression), sd_expression = sd(expression))

late_plot <- ggplot(late_df_long, aes(x = day, y = mean_expression, color = gene, group = gene, label = gene)) +
    geom_point() +
    geom_line(size = 1) +
    geom_errorbar(
        aes(
            ymin = mean_expression - sd_expression,
            ymax = mean_expression + sd_expression
        ),
        width = 0, # Width of the horizontal bar on the error bar
        size = 1 # Thickness of the error bars
    ) +
    geom_text_repel(
        data = filter(late_df_long, day == "day12"),
        aes(x = as.numeric(day) + 0.75),
        nudge_x = 0,
        nudge_y = 0.1,
        direction = "y",
        size = 10,
        segment.color = NA
    ) +
    ylim(0, 3.5) +
    ylab("Scaled normalized expression") +
    custom_theme() +
    theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 20),
        plot.margin = margin(t = 10, r = 30, b = 10, l = 10)
    )
ggsave("results/images/Figure_2A/F2A_late_LON.png", late_plot, width = 15, height = 10)
