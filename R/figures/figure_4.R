# Loading packages and functions
library(Matrix)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(WGCNA)
library(tibble)
library(paletteer)

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
rawcounts <- readcounts("data/rawcounts.csv", sep = ",", header = TRUE)
rawmeta <- read.table("data/meta.csv", sep = ",", header = TRUE)

# L9C1_2 is an outlier and is removed
meta <- filter(rawmeta, type %in% c("cyclo", "ventral") & sample != "L9C1_2" & sequencing == "batch1")
rownames(meta) <- meta$sample
counts <- rawcounts[, meta$sample][which(rowSums(rawcounts[, meta$sample]) >= 25), ]
meta$`Cyclopamine dose` <- as.factor(meta$cyclo_dose_quant)

# Normalizing data using DESeq2 variance stabilizing transformation
# Making DESeq object with lineage and type as covariates for design
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~cyclo_dose_qual
)

# Normalization by variance stabilizing transformation with and without covariates
vsd_blind <- vst(dds, blind = TRUE)
vsd <- vst(dds, blind = FALSE)

######################################
######################################
# FIGURE 4A : PCA, variance percentage and covariates ANOVA

# PCA plot of top 3000 most variable genes (DESeq2 default = 500)
pca.data <- plotPCA.DESeqTransform(vsd_blind, intgroup = c("type", "cyclo_dose_qual", "Cyclopamine dose"), returnData = TRUE, ntop = 3000)
percentVar <- round(100 * attr(pca.data, "percentVar"))

pca_plot <- ggplot(pca.data, aes(PC1, PC2, color = type, shape = Cyclopamine.dose)) +
    geom_point(size = 3, stroke = 2) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_color_manual(values = c("#ecb039", "#80AD3C")) +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6)) +
    custom_theme()
ggsave("results/images/Figure_4/Figure4supp_PCA.png", plot = pca_plot, width = 1600, height = 1200, units = "px", dpi = 250)

# Get PCA/covariates ANOVA results
PC_covariate_ANOVA <- pca_anova(
    pca_data = pca.data,
    metadata = meta,
    covariates = c("cyclo_dose_qual", "type")
)
#  Saving ANOVA results
write.csv(PC_covariate_ANOVA, "results/tables/Figure_4/Figure4_ANOVA.csv")

# PCs variation percentages :
percentVar_plot <- ggplot(data.frame(perc = percentVar, PC = factor(colnames(pca.data[1:15]), levels = colnames(pca.data[1:15]))), aes(x = PC, y = perc)) +
    geom_bar(stat = "identity") +
    custom_theme(diag_text = TRUE) +
    ylim(0, 100) +
    ylab("% of variance explained") +
    xlab("Principal Components")
ggsave("results/images/Figure_4/Figure4supp_percentVar.png", plot = percentVar_plot, width = 10, height = 5)

######################################
######################################
# FIGURE 4 B/C : WGCNA analysis, Heatmap and GO enrichment

# keeping only genes with higher variance (50% quantile)
vsd_var <- assay(vsd)[which(rowVars(assay(vsd)) >= quantile(apply(assay(vsd), 1, var), 0.5)), ]

# WGCNA analysis
# soft thresholding
sft <- pickSoftThreshold(t(vsd_var),
    powerVector = c(1:20),
    verbose = 5
)
# taking a soft trehsold of 12 :
# highest SFT.R.sq and mean.k. < 100
sft$fitIndices

#  perform WGCNA
adj <- adjacency(t(vsd_var), power = 12)
TOM <- TOMsimilarity(adj)
dissTOM <- 1 - TOM
gene.tree <- hclust(as.dist(dissTOM), method = "average")
dynamic.modules <- cutreeDynamic(
    dendro = gene.tree,
    distM = dissTOM,
    deepSplit = 2,
    pamRespectsDendro = FALSE,
    minClusterSize = 30
)
dynamic.colors <- labels2colors(dynamic.modules)

# merging similar modules
MEList <- moduleEigengenes(t(vsd_var),
    colors = dynamic.colors
)
MEs <- MEList$eigengenes
ME.diss <- 1 - cor(MEs)

merge <- mergeCloseModules(t(vsd_var),
    dynamic.colors,
    cutHeight = 0.25,
    iterate = TRUE,
    verbose = 3
)
mergedColors <- merge$colors


# get cluster into a dataframe
merged_clusters <- rownames(vsd_var)
merged_clusters_sym <- gene_converter(merged_clusters, "ENSEMBL", "SYMBOL")
cluster_df <- data.frame(ENSEMBLE = merged_clusters, gene = merged_clusters_sym, module = dynamic.colors, merged_module = mergedColors)

# get blue module (SHH modules) dataframe with weighted correlation values
blue_clust_df_f <- filter(cluster_df, merged_module == "blue", !is.na(gene))
blue_clust_df_f$cor_abs <- adj["ENSG00000164690", blue_clust_df_f$ENSEMBLE]
# get back the sign of the correlation
blue_clust_df_f$cor <- sapply(c(1:nrow(blue_clust_df_f)), function(x) {
    if (cor(assay(vsd)["ENSG00000164690", ], assay(vsd)[blue_clust_df_f$ENSEMBLE[x], ]) > 0) {
        return(blue_clust_df_f$cor_abs[x])
    } else {
        return(-blue_clust_df_f$cor_abs[x])
    }
})

# get get the top 500 genes most co-expressed with SHH in the WGCNA modules to show on the STRINGdb plot
write.csv(blue_clust_df_f[order(blue_clust_df_f$cor_abs, decreasing = TRUE), ] %>% head(500), "results/tables/Figure_3/SHH_cluster_500.csv")
write.csv(blue_clust_df_f[order(blue_clust_df_f$cor_abs, decreasing = TRUE), ], "results/tables/Figure_3/SHH_cluster.csv")

# Heatmap of the genes in the blue cluster
# Preparation for heatmap, clustering and GO enrichment
sample_order <- meta$sample[order(meta$cyclo_dose_quant)]
scaled_mat <- t(apply(assay(vsd)[filter(cluster_df, merged_module == "blue", !is.na(gene))$ENSEMBLE, sample_order], 1, scale))
colnames(scaled_mat) <- colnames(assay(vsd)[, sample_order])

# hierarchical clustering using euclidian distance and "complete" method
clustering <- hclust(dist(scaled_mat))
clusters <- cutree(clustering, k = 2)

# renaming clusters
clusters <- ifelse(clusters == 1, "Cluster 1", "Cluster 2")

# Making genelist for Figure4B
fig_genelist_5 <- clusters %>% as.data.frame()
colnames(fig_genelist_5) <- "Figure4B"
write.csv(fig_genelist_5, "results/tables/Figure_4/fig_genelist_5.csv")

row_split <- factor(
    clusters[clustering$order],
    levels = c("Cluster 2", "Cluster 1")
)

# Colors for each group
group_colors <- c(
    "Cluster 2" = "#4d6da5",
    "Cluster 1" = "#b16060"
)

# Block annotation for color blocks (without outline)
sc_percent <- read.csv("results/tables/Figure_4/sc_percent.csv", header = TRUE)

sc_percent$gene <- gene_converter(sc_percent$gene, "SYMBOL", "ENSEMBL")
missing <- setdiff(rownames(scaled_mat), sc_percent$gene)

sc_percent <- sc_percent <- rbind(sc_percent, data.frame(gene = missing, `W3-1` = 0, `W4-1` = 0, `W5-1` = 0))
rownames(sc_percent) <- sc_percent$gene
sc_percent <- sc_percent[names(clusters[clustering$order]), ]

clusters_ha <- rowAnnotation(
    clusters = anno_block(
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
    ),
    W3 = sc_percent$W3.1,
    W4 = sc_percent$W4.1,
    W5 = sc_percent$W5.1
)

WGCNA_ht_plot <- Heatmap(
    scaled_mat[clustering$order, sample_order],
    name = "Normalized expression",
    column_names_gp = gpar(fontsize = 6),
    row_split = row_split,
    row_title = NULL,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    right_annotation = clusters_ha,
    show_row_names = FALSE,
    row_names_side = "left",
    show_column_names = TRUE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
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

png_save(
    plot = WGCNA_ht_plot,
    path = "results/images/Figure_4/Figure4A.png",
    width = 1800,
    height = 1600
)

# GO term plot for Cluster 1
GO_terms_Clust1 <- plot_go_term(
    names(clusters[which(clusters %in% "Cluster 1")]),
    path = "results/images/Figure_4/Figure4B_Clust1",
    range = c(1:10),
    imgw = 22,
    imgh = 16,
    cut = 45
)
write.csv(GO_terms_Clust1, "results/tables/Figure_4/GO_terms_Clust1.csv")

GO_terms_Clust2 <- plot_go_term(
    names(clusters[which(clusters %in% "Cluster 2")]),
    path = "results/images/Figure_4/Figure4B_Clust2",
    range = c(1, 2, 3, 4, 5, 6, 7, 9, 10, 11),
    imgw = 22,
    imgh = 16,
    cut = 45
)
write.csv(GO_terms_Clust2, "results/tables/Figure_4/GO_terms_Clust2.csv")

####################################
####################################
# FIGURE SUPP : Cyclopamine differential expression

dds_WGCNA <- DESeqDataSetFromMatrix(
    countData = counts[filter(cluster_df, merged_module == "blue", !is.na(gene))$ENSEMBLE, ],
    colData = meta,
    design = ~cyclo_dose_qual
)

DE_WGCNA <- dds_WGCNA %>%
    DESeq() %>%
    results(alpha = 0.01, contrast = c("cyclo_dose_qual", "high", "none")) %>%
    as.data.frame()
DE_WGCNA$cluster <- clusters[rownames(DE_WGCNA)]


DE_WGCNA %>% View()
DE_WGCNA$gene <- gene_converter(rownames(DE_WGCNA), "ENSEMBL", "SYMBOL")
DE_WGCNA_f <- filter(DE_WGCNA, padj < 0.05 & abs(log2FoldChange) >= 2 & !is.na(gene))

# DESEQ2 object with only cyclo dose as covariate
dds_quant <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~cyclo_dose_quant
)
DE_cyclo <- dds_quant %>%
    DESeq() %>%
    results(alpha = 0.01, name = "cyclo_dose_quant") %>%
    as.data.frame()
DE_cyclo %>% View()
DE_cyclo$gene <- gene_converter(rownames(DE_cyclo), "ENSEMBL", "SYMBOL")
DE_cyclo_f <- filter(DE_cyclo, padj < 0.05 & abs(log2FoldChange) >= 2 & !is.na(gene))


# DEGs high vs no cyclo, high vs low cyclo, low vs no cyclo
DE_high_vs_no <- dds %>%
    DESeq() %>%
    results(alpha = 0.01, contrast = c("cyclo_dose_qual", "high", "none")) %>%
    as.data.frame()
DE_high_vs_no$gene <- gene_converter(rownames(DE_high_vs_no), "ENSEMBL", "SYMBOL")
DE_high_vs_no_f <- filter(DE_high_vs_no, padj < 0.05 & abs(log2FoldChange) >= 2 & !is.na(gene))

DE_high_vs_low <- dds %>%
    DESeq() %>%
    results(alpha = 0.01, contrast = c("cyclo_dose_qual", "high", "low")) %>%
    as.data.frame()
DE_high_vs_low$gene <- gene_converter(rownames(DE_high_vs_low), "ENSEMBL", "SYMBOL")
DE_high_vs_low_f <- filter(DE_high_vs_low, padj < 0.05 & abs(log2FoldChange) >= 2 & !is.na(gene))

DE_low_vs_no <- dds %>%
    DESeq() %>%
    results(alpha = 0.01, contrast = c("cyclo_dose_qual", "low", "none")) %>%
    as.data.frame()
DE_low_vs_no$gene <- gene_converter(rownames(DE_low_vs_no), "ENSEMBL", "SYMBOL")
DE_low_vs_no_f <- filter(DE_low_vs_no, padj < 0.05 & abs(log2FoldChange) >= 2 & !is.na(gene))

# Building heatmap genes annotation (DE, genes in WGCNA blue cluster,  heatmap cluster)
colnames(DE_high_vs_no) <- paste0("HvsN_", colnames(DE_high_vs_no))
colnames(DE_high_vs_low) <- paste0("HvsL_", colnames(DE_high_vs_low))
colnames(DE_low_vs_no) <- paste0("LvsN_", colnames(DE_low_vs_no))
colnames(DE_cyclo) <- paste0("cyclo_", colnames(DE_cyclo))

cyclo_genes_df <- Reduce(cbind, list(
    dplyr::select(DE_high_vs_no, c("HvsN_padj", "HvsN_log2FoldChange")),
    dplyr::select(DE_high_vs_low, c("HvsL_padj", "HvsL_log2FoldChange")),
    dplyr::select(DE_low_vs_no, c("LvsN_padj", "LvsN_log2FoldChange")),
    dplyr::select(DE_cyclo, c("cyclo_padj", "cyclo_log2FoldChange"))
))
cyclo_genes_df$HvsN_thr <- ifelse(abs(cyclo_genes_df$HvsN_log2FoldChange) >= 1 & cyclo_genes_df$HvsN_padj < 0.01, TRUE, FALSE)
cyclo_genes_df$HvsL_thr <- ifelse(abs(cyclo_genes_df$HvsL_log2FoldChange) >= 1 & cyclo_genes_df$HvsL_padj < 0.01, TRUE, FALSE)
cyclo_genes_df$LvsN_thr <- ifelse(abs(cyclo_genes_df$LvsN_log2FoldChange) >= 1 & cyclo_genes_df$LvsN_padj < 0.01, TRUE, FALSE)
cyclo_genes_df$cyclo_thr <- ifelse(abs(cyclo_genes_df$cyclo_log2FoldChange) >= 1 & cyclo_genes_df$cyclo_padj < 0.01, TRUE, FALSE)

cyclo_genes_df$WGCNA <- ifelse(rownames(cyclo_genes_df) %in% blue_clust_df_f$ENSEMBLE, "yes", "no")
cyclo_genes_df$cluster <- rownames(cyclo_genes_df) %>% sapply(function(gene) {
    if (gene %in% names(clusters)) {
        return(clusters[gene])
    } else {
        return("no")
    }
})
cyclo_genes_df$sub_cluster <- rownames(cyclo_genes_df) %>% sapply(function(gene) {
    if (gene %in% names(sub_clusters)) {
        return(sub_clusters[gene])
    } else {
        return("no")
    }
})
cyclo_genes_df %>% View()
cyclo_genes_df$genes <- gene_converter(rownames(cyclo_genes_df), "ENSEMBL", "SYMBOL")

# Saving heatmap genes annotation table
write.csv(cyclo_genes_df, "results/tables/Figure_4/cyclo_genes_df.csv")
save(cyclo_genes_df, file = "results/tables/Figure_4/cyclo_genes_df.RData")


####################################
####################################
# FIGURE SUPP : comparisons with kinetic
shh_cluster <- read.table("results/tables/Figure_4/SHH_cluster.csv", sep = ",", header = TRUE)
kinetic_genes_df <- read.table("results/tables/Figure_2/DEGpattern_WTC.csv", sep = ",", header = TRUE)

cyclo_ventral_genes <- filter(shh_cluster, cor > 0)$gene
kinetic_genes <- filter(kinetic_genes_df, cluster %in% c(1, 2))$symbol

# Create a dataframe with unique genes and their category
unique_genes <- unique(c(cyclo_ventral_genes, kinetic_genes))
gene_category <- sapply(unique_genes, function(gene) {
    if (gene %in% cyclo_ventral_genes & gene %in% kinetic_genes) {
        return("both")
    } else if (gene %in% cyclo_ventral_genes) {
        return("cyclo")
    } else {
        return("kinetic")
    }
})

gene_df <- data.frame(
    gene = unique_genes,
    category = gene_category
)

# Save the dataframe
write.csv(gene_df, "results/tables/Figure_4/intersect_genes.csv")

# Venn diagram of cyclo_ventral_genes and kinetic_genes
library(VennDiagram)
venn.plot <- venn.diagram(
    x = list(
        CycloVentral = cyclo_ventral_genes,
        Kinetic = kinetic_genes
    ),
    category.names = c("Cyclo Ventral Genes", "Kinetic Genes"),
    filename = NULL,
    output = TRUE,
    imagetype = "png",
    height = 800,
    width = 800,
    resolution = 300,
    compression = "lzw",
    lwd = 2,
    col = c("red", "blue"),
    fill = c("red", "blue"),
    alpha = 0.5,
    cex = 1.5,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 1.5,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-20, 20),
    cat.dist = c(0.05, 0.05),
    cat.fontfamily = "sans",
    cat.col = c("red", "blue")
)
venn.plot
# Save Venn diagram
png("results/images/Figure_4/Figure4supp_VennDiagram.png", width = 1800, height = 1800, res = 250)
grid.draw(venn.plot)
dev.off()



################
################
# FIGURE SUPP : DE pattern
pat_meta <- meta
rownames(pat_meta) <- meta$sample
pat_meta$cyclo_dose_quant <- as.factor(pat_meta$cyclo_dose_quant)

patterns <- degPatterns(
    assay(vsd)[names(clusters), ],
    meta = pat_meta,
    time = "cyclo_dose_quant",
    reduce = TRUE,
    nClusters = 10,
)

sign_comp <- list(
    c(0, 0.125),
    c(0.125, 0.25),
    c(0.25, 0.5),
    c(0.5, 1)
)

plot_list <- lapply(c(1, 2, 3, 5, 6, 7, 8, 9), function(i) {
    print(i)
    return(MyDegPlotCluster(table = filter(patterns$normalize, cluster == i), time = "cyclo_dose_quant", sign_comp = sign_comp, cluster_i = i))
})
combined_plot <- wrap_plots(plot_list, ncol = 3)
ggsave(paste0("results/images/Figure_4/Figure4suppPattern.png"), combined_plot, width = 35, height = 25)
