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

######################################
######################################
#  Data loading, cleaning and normalization
# Loading data
rawcounts <- readcounts("data/rawcounts.csv", sep = ",", header = TRUE)
rawmeta <- read.table("data/meta.csv", sep = ",", header = TRUE)

# Keeping only necessary samples
meta <- filter(rawmeta, type %in% c("dorsal", "ventral"), diff %in% c("diff9", "diff12"), CRISPR %in% c("no", "control"), sample != "L9C1_2")

# set gene markers list for Figure1D
markers <- c(
    "ENSG00000111704",
    "ENSG00000146477",
    "ENSG00000181449",
    "ENSG00000004848",
    "ENSG00000116106",
    "ENSG00000138685",
    "ENSG00000082701",
    "ENSG00000185551",
    "ENSG00000165588",
    "ENSG00000144857",
    "ENSG00000170370",
    "ENSG00000106571",
    "ENSG00000163132",
    "ENSG00000115507",
    "ENSG00000135903",
    "ENSG00000007372",
    "ENSG00000070193",
    "ENSG00000129514",
    "ENSG00000125798",
    "ENSG00000177283",
    "ENSG00000136352",
    "ENSG00000164690",
    "ENSG00000180730",
    "ENSG00000138083",
    "ENSG00000184302",
    "ENSG00000163064",
    "ENSG00000164778",
    "ENSG00000105991",
    "ENSG00000120094"
)

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

# Creating DESEQ2 dataset object
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ line + type
)
# Normalization by variance stabilizing transformation with and without covariates
vsd_blind <- vst(dds, blind = TRUE)
vsd <- vst(dds, blind = FALSE)

######################################
######################################
# FIGURE 1C : PCA, variance percentage and covariates ANOVA

# PCA plot of top 3000 most variable genes (DESeq2 default = 500)
pca.data <- plotPCA.DESeqTransform(vsd_blind, intgroup = c("line", "type"), returnData = TRUE, ntop = 3000)
percentVar <- round(100 * attr(pca.data, "percentVar"))

pca_plot <- ggplot(pca.data, aes(PC1, PC2, color = type, shape = line)) +
    geom_point(size = 5, stroke = 1) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_color_manual(values = c("#A1A1DE", "#80AD3C")) +
    custom_theme()
ggsave("results/images/Figure_1/Figure1C.png", plot = pca_plot, width = 1600, height = 1200, units = "px", dpi = 250)

# Get PCA/covariates ANOVA results
PC_covariate_ANOVA <- pca_anova(
    pca_data = pca.data,
    metadata = meta,
    covariates = c("line", "type")
)
#  Saving ANOVA results
write.csv(PC_covariate_ANOVA, "results/tables/Figure_1/Figure1_ANOVA.csv")

# PCs variation percentages :
percentVar_plot <- ggplot(data.frame(perc = percentVar, PC = factor(colnames(pca.data[1:15]), levels = colnames(pca.data[1:15]))), aes(x = PC, y = perc)) +
    geom_bar(stat = "identity") +
    custom_theme(diag_text = TRUE) +
    ylim(0, 100) +
    ylab("% of variance explained") +
    xlab("Principal Components")
ggsave("results/images/Figure_1/Figure1supp_percentVar.png", plot = percentVar_plot, width = 10, height = 5)

######################################
######################################
# FIGURE 1 D : Heatmap of chosen markers

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

marker_ht_plot <- Heatmap(
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

png_save(
    plot = marker_ht_plot,
    path = "results/images/Figure_1/Figure1D.png",
    width = 2000,
    height = 1800
)

######################################
######################################
# FIGURE 1 E : Heatmap of DEGs between ventral and dorsal

# Compute ventral vs dorsal DEGs
DEGs_vAN_vs_dAN <- dds %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "ventral", "dorsal")) %>%
    as.data.frame() %>%
    na.omit()

# Add column with gene symbols
DEGs_vAN_vs_dAN$gene <- rownames(DEGs_vAN_vs_dAN) %>% gene_converter("ENSEMBL", "SYMBOL")
DEGs_vAN_vs_dAN_f <- filter(DEGs_vAN_vs_dAN, padj < 0.01, abs(log2FoldChange) >= 1)

# Heatmap of vAN vs dAN DEGs
vsd_DEGs <- assay(vsd[rownames(DEGs_vAN_vs_dAN_f), ])
scaled_mat <- t(apply(vsd_DEGs, 1, scale))
colnames(scaled_mat) <- colnames(vsd_DEGs)

# Put sample in a specific order
ordered_samples <- c("L9D_1", "L9D_2", "L9D_3", "L9D_4", "L9D_5", "W6C12D_1", "W6C12D_2", "W6C12D_3", "W6C12D_4", "W6C12D_5", "W6C12D_6", "L9V_1", "L9V_2", "L9V_3", "L9V_4", "L9V_5", "W6C12V_1", "W6C12V_2", "W6C12V_3", "W6C12V_4", "W6C12V_5", "W6C12V_6") %>% as.factor()

# hierarchical clustering using euclidian distance and "complete" method
clustering <- hclust(dist(scaled_mat))
clusters <- cutree(clustering, k = 2)
# renaming clusters
clusters <- ifelse(clusters == 1, "vAN", "dAN")

# Getting figure_genelist_1 :
fig_genelist_1 <- clusters %>% as.data.frame()
colnames(fig_genelist_1) <- "Figure1E"
write.csv(fig_genelist_1, "results/tables/Figure_1/fig_genelist_1.csv")

row_split <- factor(clusters[clustering$order], levels = c("vAN", "dAN"))

# Colors for each group
group_colors <- c(
    "vAN" = "#b16060",
    "dAN" = "#4d6da5"
)

# Block annotation for color blocks (without outline)
DEGs_ha <- rowAnnotation(
    DEGs = anno_block(
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

# sample annotation
rownames(meta) <- meta$sample
sample_ha <- columnAnnotation(
    line = meta[ordered_samples, ]$line,
    type = meta[ordered_samples, ]$type,
    col = list(
        line = c("LON71" = "#c1c1c1", "WTC" = "#7d7d7d"),
        type = c("dorsal" = "#A1A1DE", "ventral" = "#80AD3C")
    ),
    show_legend = FALSE
)

DEG_ht_plot <- Heatmap(
    scaled_mat[clustering$order, ordered_samples],
    name = "Normalized expression",
    row_title = NULL,
    cluster_rows = FALSE,
    row_split = row_split,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    right_annotation = DEGs_ha,
    bottom_annotation = sample_ha,
    row_names_side = "left",
    show_column_names = TRUE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    width = ncol(scaled_mat) * unit(4, "mm"),
    col = colorRampPalette(c(
        "black",
        "purple",
        "orange",
        "yellow"
    ))(1000),
)

png_save(
    plot = DEG_ht_plot,
    path = "results/images/Figure_1/Figure1E.png",
    width = 1600,
    height = 1600
)

# adding cluster to DEGs table
DEGs_vAN_vs_dAN$clusters <- clusters[rownames(DEGs_vAN_vs_dAN)]
# Saving DE results
write.csv(DEGs_vAN_vs_dAN, "results/tables/Figure_1/DE_vAN_vs_dAN.csv")

######################################
######################################
# FIGURE 1 F : GO term enrichment on DEGs clusters
source("R/custom_fct.R")

# GO term plot for vAN cluster
GO_terms_vAN <- plot_go_term(
    names(clusters[which(clusters %in% "vAN")]),
    path = "results/images/Figure_1/Figure1F_vAN",
    range = c(1, 2, 3, 5, 6, 8, 13, 14),
    imgw = 18
)
write.csv(GO_terms_vAN, "results/tables/Figure_1/GO_terms_vAN.csv")
# GO term plot for dAN cluster
GO_terms_dAN <- plot_go_term(
    names(clusters[which(clusters %in% "dAN")]),
    path = "results/images/Figure_1/Figure1F_dAN",
    range = c(1, 2, 4, 6, 10, 12, 13, 15),
    imgw = 18
)
write.csv(GO_terms_dAN, "results/tables/Figure_1/GO_terms_dAN.csv")
