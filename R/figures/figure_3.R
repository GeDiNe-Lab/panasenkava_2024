# Loading packages and functions
library(Matrix)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(VennDiagram)
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
View(meta)
# making DESeq object
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~cyclo_dose_qual
)

# Normalization without covariates
vsd_blind <- vst(dds, blind = TRUE)

#  PCA
pca.data <- plotPCA.DESeqTransform(vsd_blind, intgroup = c("type", "cyclo_dose_qual"), returnData = TRUE, ntop = nrow(assay(vsd_blind)))
percentVar <- round(100 * attr(pca.data, "percentVar"))

png(filename = "results/images/Figure_3/F3_PCA_1_2.png", width = 1600, height = 1200, res = 250)
ggplot(pca.data, aes(PC1, PC2, color = type, shape = cyclo_dose_qual)) +
    geom_point(size = 2, stroke = 1) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_color_manual(values = c("#ecb039", "#80AD3C")) +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6)) +
    custom_theme() +
    ggtitle("First and second PCs of dorsal and ventral kinetics all genes")
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

# getting genes correlated to cyclopamine dose
cyclo_genes <- cor(t(assay(vsd)[, meta$sample[order(meta$cyclo_dose_quant)]]), sort(meta$cyclo_dose_quant)) %>% as.data.frame()
colnames(cyclo_genes) <- c("cor")
cyclo_genes$genes <- rownames(cyclo_genes) %>% gene_converter("ENSEMBL", "SYMBOL")
cyclo_genes$abscor <- abs(cyclo_genes$cor)
cyclo_genes_f <- filter(cyclo_genes, abscor > 0.5)

# Preparation for heatmap, clustering and GO enrichment
sample_order <- meta$sample[order(meta$cyclo_dose_quant)]
scaled_mat <- t(apply(assay(vsd)[rownames(cyclo_genes_f), sample_order], 1, scale))
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
            "1" = "red",
            "2" = "blue"
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
    write.csv(GO_enrichment, paste0("results/tables/Figure_3/GO_enrichment_cluster_", cluster, ".csv"))
    goplot <- clusterProfiler::dotplot(GO_enrichment,
        title = paste0("GO enrichment on cluster", cluster, " (biological processes only)"),
        showCategory = 15
    )
    ggsave(paste0("results/images/Figure_3/GO_enrichment_cluster_", cluster, ".png"), goplot, width = 8, height = 10)
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

# DEGs high vs no cyclo, high vs low cyclo, low vs no cyclo
DE_high_vs_no <- dds %>%
    DESeq() %>%
    results(alpha = 0.01, contrast = c("cyclo_dose_qual", "high", "no_cyclo")) %>%
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
    results(alpha = 0.01, contrast = c("cyclo_dose_qual", "low", "no_cyclo")) %>%
    as.data.frame()
DE_low_vs_no$gene <- gene_converter(rownames(DE_low_vs_no), "ENSEMBL", "SYMBOL")
DE_low_vs_no_f <- filter(DE_low_vs_no, padj < 0.05 & abs(log2FoldChange) >= 2 & !is.na(gene))

# Building heatmap genes annotation (DE, correlation with cyclo, cluster, subcluster)
colnames(DE_high_vs_no) <- paste0("HvsN_", colnames(DE_high_vs_no))
colnames(DE_high_vs_low) <- paste0("HvsL_", colnames(DE_high_vs_low))
colnames(DE_low_vs_no) <- paste0("LvsN_", colnames(DE_low_vs_no))
cyclo_genes_df <- Reduce(cbind, list(
    dplyr::select(DE_high_vs_no, c("HvsN_padj", "HvsN_log2FoldChange")),
    dplyr::select(DE_high_vs_low, c("HvsL_padj", "HvsL_log2FoldChange")),
    dplyr::select(DE_low_vs_no, c("LvsN_padj", "LvsN_log2FoldChange"))
))
cyclo_genes_df$HvsN_thr <- ifelse(abs(cyclo_genes_df$HvsN_log2FoldChange) >= 2 & cyclo_genes_df$HvsN_padj < 0.01, "DE", "no")
cyclo_genes_df$HvsL_thr <- ifelse(abs(cyclo_genes_df$HvsL_log2FoldChange) >= 2 & cyclo_genes_df$HvsL_padj < 0.01, "DE", "no")
cyclo_genes_df$LvsN_thr <- ifelse(abs(cyclo_genes_df$LvsN_log2FoldChange) >= 2 & cyclo_genes_df$LvsN_padj < 0.01, "DE", "no")

cyclo_genes_df$cyclo_cor <- ifelse(rownames(cyclo_genes_df) %in% rownames(cyclo_genes_f), "yes", "no")
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
cyclo_genes_df$genes <- gene_converter(rownames(cyclo_genes_df), "ENSEMBL", "SYMBOL")

# Saving heatmap genes annotation table
write.csv(cyclo_genes_df, "results/")

# Genes tested in mices
ventral <- c(
    "AFF2", "ATP2C2", "AUTS2", "BAHCC1", "CAPN6", "CNTN6", "EDNRA", "FOXA1", "FOXA2", "FRZB", "GPM6B", "GRIK3", "HTR1D", "LDB2", "LINC00261", "MBIP", "MPPED1", "MYRF", "NAALAD2", "NACC2", "NKX2-1", "NKX2-1-AS1", "NKX2-2", "NTN1", "PDZRN3", "PLCL1", "PNMA2", "PNRC2", "PPM1L", "PTCH1", "QKI", "RGMA", "RORA", "RPS6KA6", "RXRA", "SERPINF1", "SERPINI1", "SFRP1", "SHH", "SLC38A2", "SLC38A4", "SLIT1", "SMIM32", "SPON1", "SPTSSB", "TMTC2", "TRIM9", "USP2"
)
dorsal <- c("PAX6", "ADD3", "ATP2B1", "CNN3", "COLGALT2", "EPHA4", "FZD3", "GLI2", "GLI3", "HOMER1", "NLGN1", "NUAK2", "OPTN", "PALLD", "PDP1", "PLK2", "PRDX6", "SLC3A2", "VCL", "ZIC2", "ZIC5", "ZNF385B", "ZFHX4")

# Known SHH related genes
known_genes <- c("GLI2", "GLI3", "ZIC2", "FOXA1", "FOXA2", "NKX2-1", "PAX6", "PTCH1")

# co-expression
# getting co-expression matrix (Pearson correlation)
# using WGCNA::cor for faster computation
corr <- WGCNA::cor(t(assay(vsd)))

# getting genes correlation with SHH, NKX2-1 and PAX6
corr_df <- as.data.frame(corr[, c("ENSG00000164690", "ENSG00000136352", "ENSG00000007372")])
colnames(corr_df) <- c("SHH", "NKX2.1", "PAX6")
corr_df$abs_max_cor <- apply(corr_df, 1, function(x) max(abs(x)))
corr_df$gene <- gene_converter(rownames(corr_df), "ENSEMBL", "SYMBOL")
corr_df <- filter(corr_df, !is.na(gene))

# filtering genes with |cor| >= 0.85
corr_df_f <- filter(corr_df, abs_max_cor >= 0.85)

# Saving correlation table
write.csv(corr_df_f, file = "results/tables/Figure_3/Pearson_0_85.csv")

# Getting top 10 genes positively and negatively correlated with SHH
scaled_vsd <- assay(vsd) - min(assay(vsd))
subset_assay_pos <- scaled_vsd[rownames(corr_df_f[order(corr_df_f$SHH, decreasing = TRUE), ][2:11, ]), ]
rownames(subset_assay_pos) <- rownames(subset_assay_pos) %>% gene_converter("ENSEMBL", "SYMBOL")
subset_assay_neg <- scaled_vsd[rownames(corr_df_f[order(corr_df_f$SHH, decreasing = FALSE), ][1:10, ]), ]
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

# Loading the table obtained through STRING export "as tabular text output"
STRING_edges <- read.table("/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_3/string_interactions.tsv", header = FALSE)
colnames(STRING_edges) <- c("gene1", "gene2", "score")
corr_df_f_symbol <- corr_df_f
rownames(corr_df_f_symbol) <- corr_df_f$gene

# Adding genes showing high correlation but which are not in STRING since they are non-coding transcripts
# Confidence score for these relationships is set to 0.4 which is the minimum here
SHH_high_cor <- data.frame(
    gene1 = setdiff(filter(corr_df_f, abs(SHH) >= 0.95)$gene, union(unique(STRING_edges$gene1), unique(STRING_edges$gene2))),
    gene2 = "SHH",
    score = 0.4
)
NKX2.1_high_cor <- data.frame(
    gene1 = setdiff(filter(corr_df_f, abs(NKX2.1) >= 0.95)$gene, union(unique(STRING_edges$gene1), unique(STRING_edges$gene2))),
    gene2 = "NKX2-1",
    score = 0.4
)

STRING_edges <- Reduce(rbind, list(STRING_edges, SHH_high_cor, NKX2.1_high_cor))

# Adding SHH relationship sign info
STRING_edges$SHH_cor <- sapply(STRING_edges$gene1, function(x) {
    return(corr_df_f_symbol[x, "SHH"] > 0)
})
write.csv(STRING_edges, file = "results/tables/Figure_3/STRING_edges_corr.csv", row.names = FALSE)

#  Plotting higly DE genes volcano plots
DE_cyclo <- list(
    high_vs_no = DE_high_vs_no_f,
    high_vs_low = DE_high_vs_low_f,
    low_vs_no = DE_low_vs_no_f
)

for (contrast in names(DE_cyclo)) {
    print(contrast)
    DE <- DE_cyclo[[contrast]]
    ggplot(DE, aes(x = log2FoldChange, y = -log10(padj), label = gene)) +
        geom_text(size = 2) +
        custom_theme() +
        geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
        labs(
            x = "log2FoldChange",
            y = "-log10(padj)",
            title = paste0("DE: ", contrast, " ", nrow(DE), " DE genes total"),
            subtitle = "|log2FC| >= 2 & FDR < 0.01"
        )
    ggsave(filename = paste0("results/images/Figure_3/Volcano_", contrast, ".png"), units = "px", width = 1800, height = 1400, dpi = 250)
}
