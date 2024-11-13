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

# filtering out lowly expressed genes and keeping only seleted samples
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
pca.data <- plotPCA.DESeqTransform(vsd_blind, intgroup = c("line", "type"), returnData = TRUE, ntop = 3000)
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

# heatmap gene annotations
marker_ha <- rowAnnotation(
    markers = c(
        rep("pluripotency", 2),
        rep("anterior neuroectoderm", 7),
        rep("dorsal", 7),
        rep("ventral", 9),
        rep("posterior neuroectoderm", 4)
    ),
    col = list(
        markers = c(
            "pluripotency" = "#b16060",
            "anterior neuroectoderm" = "#4d6da5",
            "ventral" = "#5e9a5e",
            "dorsal" = "#78588c",
            "posterior neuroectoderm" = "#d09322"
        )
    ),
    show_annotation_name = FALSE
)
# heatmap sample annotation
sample_ha <- columnAnnotation(
    line = meta[order(meta$type), ]$line,
    type = meta[order(meta$type), ]$type,
    col = list(
        line = c("LON71" = "#c1c1c1", "WTC" = "#7d7d7d"),
        type = c("dorsal" = "#A1A1DE", "ventral" = "#80AD3C")
    )
)
png(filename = "results/images/Figure_1/F1_1_marker_HM.png", width = 2000, height = 1800, res = 250)
Heatmap(
    vsd_symbol[, meta[order(meta$type), ]$sample],
    name = "Normalized expression",
    row_title_gp = gpar(fontsize = 16, fontface = "bold"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_side = "left",
    show_column_names = TRUE,
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



DEGs_LON <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, line == "LON71")$sample],
    colData = filter(meta, line == "LON71"),
    design = ~type
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "ventral", "dorsal")) %>%
    as.data.frame() %>%
    na.omit()
# Add "LON" prefix to column names
colnames(DEGs_LON) <- paste0("LON", colnames(DEGs_LON))

DEGs_WTC <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, line == "WTC")$sample],
    colData = filter(meta, line == "WTC"),
    design = ~type
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "ventral", "dorsal")) %>%
    as.data.frame() %>%
    na.omit()
colnames(DEGs_WTC) <- paste0("WTC", colnames(DEGs_WTC))

# Get genes DE in a different way in the two lines
common <- intersect(rownames(DEGs_WTC), rownames(DEGs_LON))
DEGs_LON_WTC <- cbind(DEGs_LON[common, ], DEGs_WTC[common, ])
invert <- filter(DEGs_LON_WTC, LONlog2FoldChange * WTClog2FoldChange < 0)

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


# Heatmap of vAN vs dAN DEGs
vsd_DEGs <- assay(vsd[rownames(DEGs_vAN_vs_dAN_f), ])
scaled_mat <- t(apply(vsd_DEGs, 1, scale))
colnames(scaled_mat) <- colnames(vsd_DEGs)

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
clusters %>% table()

# sample annotation
sample_ha <- columnAnnotation(
    line = meta$line,
    type = meta$type,
    col = list(
        line = c("LON71" = "#c1c1c1", "WTC" = "#7d7d7d"),
        type = c("dorsal" = "#A1A1DE", "ventral" = "#80AD3C")
    )
)

# performing GO enrichment on clusters
for (cluster in unique(clusters)) {
    cluster <- 2
    GO_enrichment <- clusterProfiler::enrichGO(names(clusters[which(clusters == cluster)]),
        OrgDb = "org.Hs.eg.db",
        keyType = "ENSEMBL",
        ont = "BP"
    )
    GO_results <- GO_enrichment@result
    GO_results$GeneRatio <- sapply(GO_enrichment@result$GeneRatio, function(x) {
        eval(parse(text = x))
    }) %>% unname()
    GO_results_f <- GO_results[order(GO_results$GeneRatio, decreasing = TRUE)[c(1, 2, 4, 5, 8, 9, 11, 12, 13, 14)], ]

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
    goplot
    write.csv(GO_enrichment, paste0("results/tables/Figure_1/GO_enrichment_cluster_", cluster, ".csv"))
    ggsave(paste0("results/images/Figure_1/F1_DE_GO_clust", cluster, ".png"), goplot, width = 17, height = 10)
}

png(filename = "results/images/Figure_1/F1_3_DE_HM_noinvert.png", width = 1600, height = 1600, res = 250)
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

# adding cluster and subcluster to DEGs table
DEGs_vAN_vs_dAN$clusters <- clusters[rownames(DEGs_vAN_vs_dAN)]
# DEGs_vAN_vs_dAN$sub_clusters <- sub_clusters[rownames(DEGs_vAN_vs_dAN)]

write.csv(DEGs_vAN_vs_dAN, "results/tables/Figure_1/DEGs_vAN_vs_dAN.csv")


DEG_vAN_vs_dAN_LON <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, line == "LON71")$sample],
    colData = filter(meta, line == "LON71"),
    design = ~type
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "ventral", "dorsal")) %>%
    as.data.frame() %>%
    na.omit()
DEG_vAN_vs_dAN_LON$gene <- rownames(DEG_vAN_vs_dAN_LON) %>% gene_converter("ENSEMBL", "SYMBOL")
DEG_vAN_vs_dAN_LON_f <- filter(DEG_vAN_vs_dAN_LON, padj < 0.01, abs(log2FoldChange) >= 1, !is.na(gene))

DEG_vAN_vs_dAN_WTC <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, line == "WTC")$sample],
    colData = filter(meta, line == "WTC"),
    design = ~type
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "ventral", "dorsal")) %>%
    as.data.frame() %>%
    na.omit()
DEG_vAN_vs_dAN_WTC$gene <- rownames(DEG_vAN_vs_dAN_WTC) %>% gene_converter("ENSEMBL", "SYMBOL")
DEG_vAN_vs_dAN_WTC_f <- filter(DEG_vAN_vs_dAN_WTC, padj < 0.01, abs(log2FoldChange) >= 1, !is.na(gene))

common_DE <- intersect(rownames(DEG_vAN_vs_dAN_LON_f), rownames(DEG_vAN_vs_dAN_WTC_f))


vsd_DEGs <- assay(vsd[common_DE, ])
scaled_mat <- t(apply(vsd_DEGs, 1, scale))
colnames(scaled_mat) <- colnames(vsd_DEGs)

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
clusters %>% table()

# sample annotation
sample_ha <- columnAnnotation(
    line = meta$line,
    type = meta$type,
    col = list(
        line = c("LON71" = "#c1c1c1", "WTC" = "#7d7d7d"),
        type = c("dorsal" = "#A1A1DE", "ventral" = "#80AD3C")
    )
)

png(filename = "results/images/Figure_1/F1_3_DE_HM_variante.png", width = 1600, height = 1600, res = 250)
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
