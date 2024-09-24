# Loading packages and functions
library(Matrix)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(WGCNA)
library(tibble)

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

meta$cyclo_dose_quant_factor <- as.factor(meta$cyclo_dose_quant)
# making DESeq object
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~cyclo_dose_qual
)

# Normalization without covariates
vsd_blind <- vst(dds, blind = TRUE)

#  PCA
pca.data <- plotPCA.DESeqTransform(vsd_blind, intgroup = c("type", "cyclo_dose_qual", "cyclo_dose_quant_factor"), returnData = TRUE, ntop = 3000)
percentVar <- round(100 * attr(pca.data, "percentVar"))

png(filename = "results/images/Figure_3/F3_PCA_1_2.png", width = 1600, height = 1200, res = 250)
ggplot(pca.data, aes(PC1, PC2, color = type, shape = cyclo_dose_quant_factor)) +
    geom_point(size = 2, stroke = 1) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_color_manual(values = c("#ecb039", "#80AD3C")) +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6)) +
    custom_theme()
dev.off()

# Building matrix with first 5 PC and covariates
PC_covariate <- cbind(pca.data[, 1:5], meta %>%
    dplyr::select(c("cyclo_dose_qual", "type")) %>%
    apply(2, function(x) {
        return(as.numeric(factor(x)) - 1)
    }) %>%
    as.matrix())

# Correlation between PC and covariates and ANOVA
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

# Saving ANOVA results
write.csv(PC_covariate_ANOVA, "results/tables/Figure_3/F3_PC_covariate_ANOVA.csv")

rownames(PC_covariate_cor) <- paste0(rownames(PC_covariate_cor), " (", percentVar[1:5], "%)")

png(filename = "results/images/Figure_3/F3_PC_covariate_correlation.png", width = 2000, height = 1800, res = 250)
Heatmap(
    PC_covariate_cor,
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.2f", PC_covariate_cor[i, j]), x, y, gp = gpar(fontsize = 10, fontface = "bold", col = "#646464"))
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
    width = ncol(PC_covariate_cor) * unit(2.5, "cm"),
    height = nrow(PC_covariate_cor) * unit(1.5, "cm"),
    col = colorRampPalette(c(
        "lightblue",
        "darkblue"
    ))(1000),
)
dev.off()

#  Variance explained by each PC
png(filename = "results/images/Figure_3/F3_PCA_percentVar.png", width = 1600, height = 1200, res = 250)
ggplot(data.frame(perc = percentVar, PC = factor(colnames(pca.data[1:20]), levels = colnames(pca.data[1:20]))), aes(x = PC, y = perc)) +
    geom_bar(stat = "identity") +
    custom_theme(diag_text = TRUE) +
    ylim(0, 100) +
    ggtitle("Variation explained by each PC")
dev.off()

# Normalization with covariates
vsd <- vst(dds, blind = FALSE)

# Genes tested in mices
ventral <- c(
    "AFF2", "ATP2C2", "AUTS2", "BAHCC1", "CAPN6", "CNTN6", "EDNRA", "FOXA1", "FOXA2", "FRZB", "GPM6B", "GRIK3", "HTR1D", "LDB2", "LINC00261", "MBIP", "MPPED1", "MYRF", "NAALAD2", "NACC2", "NKX2-1", "NKX2-1-AS1", "NKX2-2", "NTN1", "PDZRN3", "PLCL1", "PNMA2", "PNRC2", "PPM1L", "PTCH1", "QKI", "RGMA", "RORA", "RPS6KA6", "RXRA", "SERPINF1", "SERPINI1", "SFRP1", "SHH", "SLC38A2", "SLC38A4", "SLIT1", "SMIM32", "SPON1", "SPTSSB", "TMTC2", "TRIM9", "USP2"
)
dorsal <- c("PAX6", "ADD3", "ATP2B1", "CNN3", "COLGALT2", "EPHA4", "FZD3", "GLI2", "GLI3", "HOMER1", "NLGN1", "NUAK2", "OPTN", "PALLD", "PDP1", "PLK2", "PRDX6", "SLC3A2", "VCL", "ZIC2", "ZIC5", "ZNF385B", "ZFHX4")

# Known SHH related genes
known_genes <- c("GLI2", "GLI3", "ZIC2", "FOXA1", "FOXA2", "NKX2-1", "PAX6", "PTCH1", "SHH")

# keeping only genes with higher variance (50% quantile)
vsd_var <- assay(vsd)[which(rowVars(assay(vsd)) >= quantile(apply(assay(vsd), 1, var), 0.5)), ]

# WGCNA analysis
# soft thresholding
sft <- pickSoftThreshold(t(vsd_var),
    powerVector = c(1:20),
    verbose = 5
)
ggplot(sft$fitIndices, aes(Power, SFT.R.sq, label = Power)) +
    geom_point() +
    geom_text(nudge_y = 0.1) +
    geom_hline(yintercept = 0.8, color = "red") +
    labs(x = "Power", y = "Scale free topology model fit, signed R^2") +
    theme_classic()

ggplot(sft$fitIndices, aes(Power, mean.k., label = Power)) +
    geom_point() +
    geom_text(nudge_y = 0.1) +
    labs(x = "Power", y = "Mean connectivity") +
    theme_classic()

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
View(cluster_df)
# get blue module (shh modules) dataframe with weighted correlation values
blue_clust_df_f <- filter(cluster_df, merged_module == "blue", !is.na(gene))
table(cluster_df$merged_module)
blue_clust_df_f$cor_abs <- adj["ENSG00000164690", blue_clust_df_f$ENSEMBLE]
# get back the sign of the correlation
blue_clust_df_f$cor <- sapply(c(1:nrow(blue_clust_df_f)), function(x) {
    if (cor(assay(vsd)["ENSG00000164690", ], assay(vsd)[blue_clust_df_f$ENSEMBLE[x], ]) > 0) {
        return(blue_clust_df_f$cor_abs[x])
    } else {
        return(-blue_clust_df_f$cor_abs[x])
    }
})
View(blue_clust_df_f)

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

# Subclustering of each cluster
sub_clusters_list <- unique(clusters) %>% lapply(function(cluster) {
    sub_mat <- scaled_mat[names(clusters[which(clusters == cluster)]), sample_order]
    sub_clustering <- hclust(dist(sub_mat))
    return(cutree(sub_clustering, k = 5))
})
names(sub_clusters_list) <- paste0("cluster_", unique(clusters))
sub_clusters <- sub_clusters_list %>%
    unname() %>%
    unlist()
sub_clusters <- sub_clusters[names(clusters)]

# Heatmap gene annotation
clusters_ha <- rowAnnotation(
    cluster = as.character(clusters[clustering$order]),
    sub_cluster = as.character(sub_clusters[clustering$order]),
    col = list(
        cluster = c(
            "1" = "#b16060",
            "2" = "#4d6da5"
        ),
        sub_cluster = c(
            "1" = "black",
            "2" = "pink",
            "3" = "yellow",
            "4" = "brown",
            "5" = "grey"
        )
    )
)

# GO enrichment for each cluster
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
    write.csv(GO_enrichment, paste0("results/tables/Figure_3/GO_enrichment_cluster_", cluster, ".csv"))
    ggsave(paste0("results/images/Figure_3/GO_enrichment_cluster_", cluster, ".png"), goplot, width = 15, height = 10)
}

png(filename = "results/images/Figure_3/F3_cyclo_genes_HM.png", width = 2400, height = 1600, res = 250)
Heatmap(
    scaled_mat[clustering$order, sample_order],
    name = "Normalized expression",
    column_names_gp = gpar(fontsize = 6),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    left_annotation = clusters_ha,
    show_row_names = FALSE,
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
View(meta)

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

# Building heatmap genes annotation (DE, genes in WGCNA blue cluster,  heatmap cluster, heatmap subcluster)
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
write.csv(cyclo_genes_df, "results/tables/Figure_3/cyclo_genes_df.csv")
save(cyclo_genes_df, file = "results/tables/Figure_3/cyclo_genes_df.RData")

filter(cyclo_genes_df, HvsN_thr == "TRUE", !is.na(genes))$genes
STRING_edges <- read.table("/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_3/string_interactions_WGCNA.tsv", header = FALSE)
View(STRING_edges)
STRING_edges <- STRING_edges[, c(1, 2, 13)]
colnames(STRING_edges) <- c("gene1", "gene2", "combined_score")

# SYMBOL on STRINGdb does not agree with consensus symbol for this ENSEMBL id
STRING_edges$gene1[which(STRING_edges$gene1 == "MYLPF")] <- "MYL11"
STRING_edges$gene1[which(STRING_edges$gene1 == "MYLPF")] <- "MYL11"
STRING_edges$gene2[which(STRING_edges$gene2 == "MYLPF")] <- "MYL11"
STRING_edges$gene2[which(STRING_edges$gene2 == "MYLPF")] <- "MYL11"

SHH_corr_t <- cor(assay(vsd)["ENSG00000164690", ], t(assay(vsd)[blue_clust_df_f$ENSEMBLE, ]))
colnames(SHH_corr_t) <- gene_converter(colnames(SHH_corr_t), "ENSEMBL", "SYMBOL")
SHH_corr <- as.vector(SHH_corr_t)
names(SHH_corr) <- colnames(SHH_corr_t)
SHH_corr
STRING_edges$type <- ifelse(SHH_corr[STRING_edges$gene1] > 0, "positive", "negative")
STRING_edges$known <- ifelse(STRING_edges$gene1 %in% known_genes, "known", "unknown")

write.csv(STRING_edges, file = "results/tables/Figure_3/STRING_edges_format_WGCNA.csv", row.names = FALSE)



scaled_vsd <- assay(vsd) - min(assay(vsd))
subset_assay_pos <- scaled_vsd[blue_clust_df_f[order(blue_clust_df_f$cor, decreasing = TRUE), ]$ENSEMBLE[2:11], ]
rownames(subset_assay_pos) <- rownames(subset_assay_pos) %>% gene_converter("ENSEMBL", "SYMBOL")
subset_assay_neg <- scaled_vsd[blue_clust_df_f[order(blue_clust_df_f$cor, decreasing = FALSE), ]$ENSEMBLE[1:10], ]
rownames(subset_assay_neg) <- rownames(subset_assay_neg) %>% gene_converter("ENSEMBL", "SYMBOL")

# Building dataframe for lineplots
df_neg <- t(subset_assay_neg) %>% as.data.frame()
df_neg$cyclo_dose_quant <- meta[rownames(df_neg), ]$cyclo_dose_quant
df_neg_mean <- df_neg %>%
    group_by(cyclo_dose_quant) %>%
    summarise(across(everything(), mean)) %>%
    as.data.frame()
df_neg_sd <- df_neg %>%
    group_by(cyclo_dose_quant) %>%
    summarise(across(everything(), sd)) %>%
    as.data.frame()
df_neg <- df_neg_mean %>% reshape2::melt(id.vars = "cyclo_dose_quant")
df_neg$sd <- reshape2::melt(df_neg_sd, id.vars = "cyclo_dose_quant")$value
df_neg$cyclo_dose_quant <- factor(df_neg$cyclo_dose_quant, levels = c(0, 0.125, 0.25, 0.5, 1))
colnames(df_neg) <- c("cyclo_dose_quant", "gene", "expression_mean", "expression_sd")

df_pos <- t(subset_assay_pos) %>% as.data.frame()
df_pos$cyclo_dose_quant <- meta[rownames(df_pos), ]$cyclo_dose_quant
df_pos_mean <- df_pos %>%
    group_by(cyclo_dose_quant) %>%
    summarise(across(everything(), mean)) %>%
    as.data.frame()
df_pos_sd <- df_pos %>%
    group_by(cyclo_dose_quant) %>%
    summarise(across(everything(), sd)) %>%
    as.data.frame()
df_pos <- df_pos_mean %>% reshape2::melt(id.vars = "cyclo_dose_quant")
df_pos$sd <- reshape2::melt(df_pos_sd, id.vars = "cyclo_dose_quant")$value
df_pos$cyclo_dose_quant <- factor(df_pos$cyclo_dose_quant, levels = c(0, 0.125, 0.25, 0.5, 1))
colnames(df_pos) <- c("cyclo_dose_quant", "gene", "expression_mean", "expression_sd")

png(filename = "results/images/Figure_3/F3_lineplot_negative.png", width = 2400, height = 1600, res = 250)
ggplot(df_neg, aes(x = cyclo_dose_quant, y = expression_mean, group = gene, color = gene, shape = gene)) +
    geom_line() +
    geom_point(size = 2) +
    scale_shape_manual(values = 1:10) +
    ylim(0, max(cbind(subset_assay_neg, subset_assay_pos))) +
    geom_errorbar(aes(ymin = expression_mean - expression_sd, ymax = expression_mean + expression_sd), width = 0.2) +
    custom_theme(diag_text = TRUE)
dev.off()

png(filename = "results/images/Figure_3/F3_lineplot_positive.png", width = 2400, height = 1600, res = 250)
ggplot(df_pos, aes(x = cyclo_dose_quant, y = expression_mean, group = gene, color = gene, shape = gene)) +
    geom_line() +
    geom_point(size = 2) +
    scale_shape_manual(values = 1:10) +
    ylim(0, max(cbind(subset_assay_neg, subset_assay_pos))) +
    geom_errorbar(aes(ymin = expression_mean - expression_sd, ymax = expression_mean + expression_sd), width = 0.2) +
    custom_theme(diag_text = TRUE)
dev.off()
