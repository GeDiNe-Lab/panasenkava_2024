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
# Loading data
rawcounts <- readcounts("data/rawcounts.csv", sep = ",", header = TRUE)
rawmeta <- read.table("data/meta.csv", sep = ",", header = TRUE)

# Keeping only necessary samples
meta <- filter(rawmeta, type %in% c("dorsal", "ventral"), diff %in% c("diff9", "diff12"), CRISPR %in% c("no", "control"), sample != "L9C1_2")

# Get defined marker genes
markers_symbol <- read.table("data/GeneListFig1_04_09_24.csv", sep = ",", header = TRUE)
markers_symbol <- markers_symbol$gene
# Getting marker's ENSEMBL ID
markers <- markers_symbol %>% gene_converter("SYMBOL", "ENSEMBL")
# Make sure to keep only markers in the matrix
markers <- intersect(markers, rownames(rawcounts))

# Get rows corresponding to markers
retained_row <- rawcounts[rownames(rawcounts) %in% markers, meta$sample]

# filtering out lowly expressed genes and keeping only selected samples
counts <- rawcounts[, meta$sample][rowSums(rawcounts[, meta$sample]) >= 25, ]

# putting back potentially filtered out markers (posterior markers for example as they should not be expressed)
if (length(which(!rownames(retained_row) %in% rownames(counts))) == 1) {
    counts <- rbind(counts, retained_row[which(!rownames(retained_row) %in% rownames(counts)), ])
    rownames(counts)[nrow(counts)] <- rownames(retained_row)[which(!rownames(retained_row) %in% rownames(counts))]
} else if (length(which(!rownames(retained_row) %in% rownames(counts))) > 1) {
    counts <- rbind(counts, retained_row[which(!rownames(retained_row) %in% rownames(counts)), ])
}

# making DESeq object with lineage and type as covariates for design
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ line + type
)
# Normalization by variance stabilizing transformation without covariates
vsd_blind <- vst(dds, blind = TRUE)

# PCA plot of top 3000 most variable genes
pca.data <- plotPCA.DESeqTransform(vsd_blind, intgroup = c("line", "type"), returnData = TRUE, ntop = nrow(vsd_blind))
percentVar <- round(100 * attr(pca.data, "percentVar"))
fe_info <- attr(pca.data, "factoextra")

png(filename = "results/images/Figure_1/F1_2_PCA.png", width = 1600, height = 1200, res = 250)
ggplot(pca.data, aes(PC1, PC2, color = type, shape = line)) +
    geom_point(size = 5, stroke = 1) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_color_manual(values = c("#A1A1DE", "#80AD3C")) +
    custom_theme() +
    ggtitle("PCA of dorsal and ventral samples at day12")
dev.off()

# Building dataframe with first 5 PC and covariates
PC_covariate <- cbind(pca.data[, 1:5], meta %>%
    dplyr::select(c("line", "type")) %>%
    apply(2, function(x) {
        return(as.numeric(factor(x)) - 1)
    }) %>%
    as.matrix())

# Computing PC-covariate correlation and ANOVA
PC_covariate_cor <- cor(PC_covariate[, 1:5], PC_covariate[, 6:ncol(PC_covariate)]) %>% abs()
PC_covariate_ANOVA <- c(6:ncol(PC_covariate)) %>% lapply(function(i) {
    apply(PC_covariate[, 1:5], 2, function(x) {
        aov(x ~ PC_covariate[, i])
    }) %>% sapply(function(x) {
        summary(x)[[1]]$`Pr(>F)`[1]
    })
})
PC_covariate_ANOVA <- Reduce(cbind, PC_covariate_ANOVA)
colnames(PC_covariate_ANOVA) <- colnames(PC_covariate)[6:ncol(PC_covariate)]

#  Saving ANOVA results
write.csv(PC_covariate_ANOVA, "results/tables/Figure_1/F1_PC_covariate_ANOVA.csv")

# Adding variance explained by each PC
rownames(PC_covariate_cor) <- paste0(rownames(PC_covariate_cor), " (", percentVar[1:5], "%)")

png(filename = "results/images/Figure_1/F1_2B_PC_covariate_correlation.png", width = 2000, height = 1800, res = 250)
Heatmap(
    PC_covariate_cor,
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.2f", PC_covariate_cor[i, j]), x, y, gp = gpar(fontsize = 10, fontface = "bold", col = "#646464"))
    },
    name = "absolute Pearson correlation",
    row_title_gp = gpar(fontsize = 20, fontface = "bold"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    column_names_side = "top",
    column_names_rot = 0,
    column_names_centered = TRUE,
    show_column_names = TRUE,
    show_heatmap_legend = TRUE,
    width = ncol(PC_covariate_cor) * unit(1.5, "cm"),
    height = nrow(PC_covariate_cor) * unit(1, "cm"),
    col = colorRampPalette(c(
        "lightblue",
        "darkblue"
    ))(1000),
)
dev.off()

# PCs variation percentages :
png(filename = "results/images/Figure_1/F1_2a_percentVar.png", width = 1600, height = 1200, res = 250)
ggplot(data.frame(perc = percentVar, PC = factor(colnames(pca.data[1:20]), levels = colnames(pca.data[1:20]))), aes(x = PC, y = perc)) +
    geom_bar(stat = "identity") +
    custom_theme(diag_text = TRUE) +
    ylim(0, 100) +
    ggtitle("Variation explained by each PC")
dev.off()

# Normalization by variance stabilizing transformation taking covariates in account
vsd <- vst(dds, blind = FALSE)

# subset rows for markers
vsd_symbol <- assay(vsd[markers, ])
# convert rownames to gene symbols
rownames(vsd_symbol) <- rownames(vsd_symbol) %>% gene_converter("ENSEMBL", "SYMBOL")

# Define groupings and labels for blocks
row_split <- factor(
    c(
        rep("iPSCs", 2),
        rep("forebrain", 7),
        rep("dorsal\nforebrain", 7),
        rep("ventral\nforebrain", 9),
        rep("Midbrain &\nHindbrain", 4)
    ),
    levels = c("iPSCs", "forebrain", "dorsal\nforebrain", "ventral\nforebrain", "Midbrain &\nHindbrain")
)

# Colors for each group
group_colors <- c(
    "iPSCs" = "#b16060",
    "forebrain" = "#4d6da5",
    "ventral\nforebrain" = "#80AD3C",
    "dorsal\nforebrain" = "#A1A1DE",
    "Midbrain &\nHindbrain" = "#d09322"
)

# Block annotation for color blocks (without outline)
marker_ha <- rowAnnotation(
    markers = anno_block(
        gp = gpar(fill = group_colors[levels(row_split)], col = NA), # Remove outline with `col = NA`
        which = "row",
        width = unit(2, "mm") # Thinner blocks
    ),
    labels = anno_block(
        gp = gpar(fill = "white", col = "white"), # Remove outline with `col = NA`
        labels = levels(row_split), # Add group labels
        labels_gp = gpar(fontsize = 10, fontface = "bold"), # Customize label appearance
        labels_rot = -90, # Horizontal labels
        labels_just = "center"
    )
)

# heatmap sample annotation
sample_ha <- columnAnnotation(
    line = meta[order(meta$type), ]$line,
    type = meta[order(meta$type), ]$type,
    col = list(
        line = c("LON71" = "#c1c1c1", "WTC" = "#7d7d7d"),
        type = c("dorsal" = "#A1A1DE", "ventral" = "#80AD3C")
    ),
    show_legend = FALSE,
    annotation_name_gp = gpar(fontsize = 0)
)


png(filename = "results/images/Figure_1/F1_1_marker_HM.png", width = 2000, height = 1800, res = 250)
Heatmap(
    vsd_symbol[, meta[order(meta$type), ]$sample],
    name = "Normalized expression",
    row_title = NULL,
    row_split = row_split,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_side = "left",
    show_column_names = FALSE,
    bottom_annotation = sample_ha,
    right_annotation = marker_ha,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    width = ncol(vsd_symbol) * unit(4, "mm"),
    height = nrow(vsd_symbol) * unit(4, "mm"),
    col = colorRampPalette(c(
        "black",
        "purple",
        "orange",
        "yellow"
    ))(1000),
)
dev.off()

# Compute ventral vs dorsal DEGs
DEGs_vAN_vs_dAN <- dds %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "ventral", "dorsal")) %>%
    as.data.frame() %>%
    na.omit()

DEGs_vAN_vs_dAN$gene <- rownames(DEGs_vAN_vs_dAN) %>% gene_converter("ENSEMBL", "SYMBOL")
# DEGs_vAN_vs_dAN_f <- DEGs_vAN_vs_dAN %>% filter(!is.na(gene))
DEGs_vAN_vs_dAN_f <- filter(DEGs_vAN_vs_dAN, padj < 0.01, abs(log2FoldChange) >= 1)
# filter out genes with inverted log2FoldChange in LON and WTC
DEGs_vAN_vs_dAN_f <- DEGs_vAN_vs_dAN_f[!rownames(DEGs_vAN_vs_dAN_f) %in% rownames(invert), ]
View(DEGs_vAN_vs_dAN_f)

# Heatmap of vAN vs dAN DEGs
vsd_DEGs <- assay(vsd[rownames(DEGs_vAN_vs_dAN_f), ])
scaled_mat <- t(apply(vsd_DEGs, 1, scale))
colnames(scaled_mat) <- colnames(vsd_DEGs)

ordered_samples <- c("L9D_1", "L9D_2", "L9D_3", "L9D_4", "L9D_5", "W6C12D_1", "W6C12D_2", "W6C12D_3", "W6C12D_4", "W6C12D_5", "W6C12D_6", "L9V_1", "L9V_2", "L9V_3", "L9V_4", "L9V_5", "W6C12V_1", "W6C12V_2", "W6C12V_3", "W6C12V_4", "W6C12V_5", "W6C12V_6") %>% as.factor()

# hierarchical clustering using euclidian distance and "complete" method
clustering <- hclust(dist(scaled_mat))
clusters <- cutree(clustering, k = 2)

# cluster and subcluster annotation
clusters_ha <- rowAnnotation(
    cluster = as.character(clusters[clustering$order]),
    col = list(
        cluster = c(
            "1" = "#b16060",
            "2" = "#4d6da5"
        )
    )
)

# sample annotation
rownames(meta) <- meta$sample
sample_ha <- columnAnnotation(
    line = meta[ordered_samples, ]$line,
    type = meta[ordered_samples, ]$type,
    col = list(
        line = c("LON71" = "#c1c1c1", "WTC" = "#7d7d7d"),
        type = c("dorsal" = "#A1A1DE", "ventral" = "#80AD3C")
    )
)



# performing GO enrichment on clusters
for (cluster in unique(clusters)) {
    cluster <- c(1, 2)
    GO_enrichment <- clusterProfiler::enrichGO(names(clusters[which(clusters %in% cluster)]),
        OrgDb = "org.Hs.eg.db",
        keyType = "ENSEMBL",
        ont = "BP"
    )
    GO_results <- GO_enrichment@result
    GO_results$GeneRatio <- sapply(GO_enrichment@result$GeneRatio, function(x) {
        eval(parse(text = x))
    }) %>% unname()
    GO_results$rank <- rank(-GO_results$GeneRatio, ties.method = "first")
    View(GO_results)
    GO_results_f <- GO_results[order(GO_results$GeneRatio, decreasing = TRUE)[1:20], ]

    GO_results_f$Description <- str_wrap(GO_results_f$Description, width = 40) %>% str_to_upper()
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
    # write.csv(GO_enrichment, paste0("results/tables/Figure_1/GO_enrichment_cluster_", cluster, ".csv"))
    ggsave(paste0("results/images/Figure_1/test", cluster, ".png"), goplot, width = 17, height = 10)
}

png(filename = "results/images/Figure_1/F1_3_DE_HM_noinvert_noclust.png", width = 1600, height = 1600, res = 250)
Heatmap(
    scaled_mat[clustering$order, ordered_samples],
    name = "Normalized expression",
    row_title_gp = gpar(fontsize = 16, fontface = "bold"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
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

# adding cluster and subcluster to DEGs table
DEGs_vAN_vs_dAN$clusters <- clusters[rownames(DEGs_vAN_vs_dAN)]

write.csv(DEGs_vAN_vs_dAN, "results/tables/Figure_1/DEGs_vAN_vs_dAN.csv")
