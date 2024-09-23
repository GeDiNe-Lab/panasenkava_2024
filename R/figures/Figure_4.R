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

# Loading data
rawcounts <- readcounts("data/rawcounts.csv", sep = ",", header = TRUE)
rawmeta <- read.table("data/meta.csv", sep = ",", header = TRUE)

#  selecting only necessery samples
meta <- filter(rawmeta, type == "ventral" & CRISPR %in% c("control", "hetero", "homo"))
View(meta)

filtercounts <- rawcounts[, filter(meta, type == "ventral")$sample][rowSums(rawcounts[, filter(meta, type == "ventral")$sample]) >= 25, ]
# putting back TBR1 since we want to look at it expression and it is filtered out by the counts filter
filtercounts <- rbind(filtercounts, rawcounts["ENSG00000136535", colnames(filtercounts)])
rownames(filtercounts)[length(rownames(filtercounts))] <- "ENSG00000136535"

# Making DESeq objects for ventral samples
dds_vAN <- DESeqDataSetFromMatrix(
    countData = filtercounts,
    colData = filter(meta, type == "ventral"),
    design = ~CRISPR
)
# Normalization with variance stabilizing transformation without covariates
vsd_vAN_blind <- vst(dds_vAN, blind = TRUE)

# PCA for ventral samples for top 3000 variable genes
pca.data_vAN <- plotPCA.DESeqTransform(vsd_vAN_blind, intgroup = c("sample", "type", "CRISPR"), returnData = TRUE, ntop = 3000, pcsToUse = 1:10)
percentVar_vAN <- round(100 * attr(pca.data_vAN, "percentVar"))
png(filename = "results/images/Figure_4/F4_PCA_vAN.png", width = 1600, height = 1200, res = 250)
ggplot(pca.data_vAN, aes(PC1, PC2, color = type, shape = CRISPR)) +
    geom_point(size = 3, stroke = 1) +
    xlab(paste0("PC1: ", percentVar_vAN[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar_vAN[2], "% variance")) +
    scale_color_manual(values = c("#80AD3C")) +
    custom_theme() +
    ggtitle("PCA of ventral samples in CRISPR line")
dev.off()

# Building matrix with first 5 PC and covariates
PC_covariate_vAN <- cbind(pca.data_vAN[, 1:5], filter(meta, type == "ventral") %>%
    dplyr::select(c("CRISPR")) %>%
    apply(2, function(x) {
        return(as.numeric(factor(x)) - 1)
    }) %>% as.data.frame())
PC_covariate_vAN$NKX21 <- assay(vsd_vAN_blind)["ENSG00000136352", ]
PC_covariate_vAN$SHH <- assay(vsd_vAN_blind)["ENSG00000164690", ]
PC_covariate_vAN <- as.matrix(PC_covariate_vAN)

# Computing correlation and ANOVA between first 5 PC and covariates
PC_covariate_vAN_cor <- cor(PC_covariate_vAN[, 1:5], PC_covariate_vAN[, 6:ncol(PC_covariate_vAN)]) %>% abs()
PC_covariate_vAN_ANOVA <- c(6:ncol(PC_covariate_vAN)) %>% lapply(function(i) {
    apply(PC_covariate_vAN[, 1:5], 2, function(x) {
        aov(x ~ PC_covariate_vAN[, i])
    }) %>% sapply(function(x) {
        summary(x)[[1]]$`Pr(>F)`[1]
    })
})
PC_covariate_vAN_ANOVA <- Reduce(cbind, PC_covariate_vAN_ANOVA)
colnames(PC_covariate_vAN_ANOVA) <- colnames(PC_covariate_vAN)[6:ncol(PC_covariate_vAN)]

# Saving ANOVA results
write.csv(PC_covariate_vAN_ANOVA, "results/tables/Figure_4/F4_PC_covariate_vAN_ANOVA.csv")

rownames(PC_covariate_vAN_cor) <- paste0(rownames(PC_covariate_vAN_cor), " (", percentVar_vAN[1:5], "%)")

png(filename = "results/images/Figure_4/F4_PC_covariate_vAN_correlation.png", width = 2000, height = 1800, res = 250)
Heatmap(
    PC_covariate_vAN_cor,
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.2f", PC_covariate_vAN_cor[i, j]), x, y, gp = gpar(fontsize = 10, fontface = "bold", col = "#646464"))
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
    width = ncol(PC_covariate_vAN_cor) * unit(2, "cm"),
    height = nrow(PC_covariate_vAN_cor) * unit(1, "cm"),
    col = colorRampPalette(c(
        "lightblue",
        "darkblue"
    ))(1000),
)
dev.off()

# Variance explained by each PC
png(filename = "results/images/Figure_4/F4_PCA_vAN_percentVar.png", width = 1600, height = 1200, res = 250)
ggplot(data.frame(perc = percentVar_vAN, PC = factor(colnames(pca.data_vAN[1:10]), levels = colnames(pca.data_vAN[1:10]))), aes(x = PC, y = perc)) +
    geom_bar(stat = "identity") +
    custom_theme(diag_text = TRUE) +
    ylim(0, 100) +
    ggtitle("Variation explained by each PC")
dev.off()


# Normalization with variance stabilizing transformation with covariates
vsd_vAN <- vst(dds_vAN, blind = FALSE)

# DEGs between control and heterozygous for ventral, dorsal, and ventral+dorsal samples
DE_vAN_control_vs_het <- dds_vAN %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("CRISPR", "control", "hetero")) %>%
    as.data.frame() %>%
    na.omit()
DE_vAN_control_vs_het$gene <- gene_converter(rownames(DE_vAN_control_vs_het), "ENSEMBL", "SYMBOL")
DE_vAN_control_vs_het_f <- filter(DE_vAN_control_vs_het, !is.na(gene))
DE_vAN_control_vs_het_f <- filter(DE_vAN_control_vs_het_f, padj < 0.01, abs(log2FoldChange) >= 2)

# DEGs between control and homozygous for ventral, dorsal, and ventral+dorsal samples
DE_vAN_control_vs_homo <- dds_vAN %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("CRISPR", "control", "homo")) %>%
    as.data.frame() %>%
    na.omit()
DE_vAN_control_vs_homo$gene <- gene_converter(rownames(DE_vAN_control_vs_homo), "ENSEMBL", "SYMBOL") %>% as.vector()
DE_vAN_control_vs_homo_f <- filter(DE_vAN_control_vs_homo, !is.na(gene))
DE_vAN_control_vs_homo_f <- filter(DE_vAN_control_vs_homo_f, padj < 0.01, abs(log2FoldChange) >= 2)

# DEGs between heterozygous and homozygous for ventral, dorsal, and ventral+dorsal samples
DE_vAN_het_vs_homo <- dds_vAN %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("CRISPR", "hetero", "homo")) %>%
    as.data.frame() %>%
    na.omit()
DE_vAN_het_vs_homo$gene <- gene_converter(rownames(DE_vAN_het_vs_homo), "ENSEMBL", "SYMBOL")
DE_vAN_het_vs_homo_f <- filter(DE_vAN_het_vs_homo, !is.na(gene))
DE_vAN_het_vs_homo_f <- filter(DE_vAN_het_vs_homo_f, padj < 0.01, abs(log2FoldChange) >= 2)

#  Saving DE tables
write.csv(DE_vAN_control_vs_het, "results/tables/Figure_4/DE_vAN_control_vs_het.csv")
write.csv(DE_vAN_control_vs_homo, "results/tables/Figure_4/DE_vAN_control_vs_homo.csv")
write.csv(DE_vAN_het_vs_homo, "results/tables/Figure_4/DE_vAN_het_vs_homo.csv")

# Making volcanp plots
DE_CRISPR <- list(
    vAN_control_vs_het = DE_vAN_control_vs_het_f,
    vAN_control_vs_homo = DE_vAN_control_vs_homo_f,
    vAN_het_vs_homo = DE_vAN_het_vs_homo_f
)

for (contrast in names(DE_CRISPR)) {
    print(contrast)
    DE <- DE_CRISPR[[contrast]]
    ggplot(DE, aes(x = log2FoldChange, y = -log10(padj), label = gene)) +
        geom_text(size = 2) +
        custom_theme() +
        geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
        geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
        labs(x = "log2FoldChange", y = "-log10(padj)", title = paste0("DE: ", contrast), subtitle = "|log2FC| >= 2 & FDR < 0.01")
    ggsave(filename = paste0("results/images/Figure_4/Volcano_", contrast, ".png"), units = "px", width = 1800, height = 1400, dpi = 250)
}

# Making marker genes barplots (and also barplots of gene tested in mices not in the WGCNA cluster)
norm <- assay(vsd_vAN)
rownames(norm) <- gene_converter(rownames(norm), "ENSEMBL", "SYMBOL")

# scaling data so the minimum normalized counts is set to 0
meta$SHH <- norm["SHH", ] - min(norm)
meta$NKX21 <- norm["NKX2-1", ] - min(norm)
meta$PAX6 <- norm["PAX6", ] - min(norm)
meta$TBR1 <- norm["TBR1", ] - min(norm)

# meta$SFTA3 <- norm["SFTA3", ] - min(norm)
meta$EPHB1 <- norm["EPHB1", ] - min(norm)
meta$TRIM9 <- norm["TRIM9", ] - min(norm)
meta$EPHA4 <- norm["EPHA4", ] - min(norm)
meta$ZFHX4 <- norm["ZFHX4", ] - min(norm)

meta$CRISPR_type <- paste(meta$CRISPR, meta$type, sep = "_")
grouped_df <- meta %>%
    group_by(CRISPR_type, type) %>%
    summarize(
        SHH_mean = mean(SHH),
        SHH_sd = sd(SHH),
        NKX21_mean = mean(NKX21),
        NKX21_sd = sd(NKX21),
        PAX6_mean = mean(PAX6),
        PAX6_sd = sd(PAX6),
        TBR1_mean = mean(TBR1),
        TBR1_sd = sd(TBR1),
        EPHB1_mean = mean(EPHB1),
        EPHB1_sd = sd(EPHB1),
        TRIM9_mean = mean(TRIM9),
        TRIM9_sd = sd(TRIM9),
        EPHA4_mean = mean(EPHA4),
        EPHA4_sd = sd(EPHA4),
        ZFHX4_mean = mean(ZFHX4),
        ZFHX4_sd = sd(ZFHX4)
    )
grouped_df <- grouped_df[order(grouped_df$type, decreasing = TRUE), ]
grouped_df$CRISPR_type <- factor(c("vAN +/+", "vAN +/-", "vAN -/-"), levels = c("vAN +/+", "vAN +/-", "vAN -/-"))

png(filename = "results/images/Figure_4/F4_SHH_barplot.png", width = 1600, height = 1400, res = 250)
ggplot(grouped_df, aes(x = CRISPR_type, y = SHH_mean, fill = type)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = SHH_mean - SHH_sd, ymax = SHH_mean + SHH_sd), width = 0.2) +
    ylim(-1, 7) +
    scale_fill_manual(values = c("#80AD3C")) +
    custom_theme(diag_text = TRUE, hide_legend = TRUE)
dev.off()

png(filename = "results/images/Figure_4/F4_NKX21_barplot.png", width = 1600, height = 1400, res = 250)
ggplot(grouped_df, aes(x = CRISPR_type, y = NKX21_mean, fill = type)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = NKX21_mean - NKX21_sd, ymax = NKX21_mean + NKX21_sd), width = 0.2) +
    ylim(-1, 7) +
    scale_fill_manual(values = c("#80AD3C")) +
    custom_theme(diag_text = TRUE, hide_legend = TRUE)
dev.off()

png(filename = "results/images/Figure_4/F4_PAX6_barplot.png", width = 1600, height = 1400, res = 250)
ggplot(grouped_df, aes(x = CRISPR_type, y = PAX6_mean, fill = type)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = PAX6_mean - PAX6_sd, ymax = PAX6_mean + PAX6_sd), width = 0.2) +
    ylim(-1, 7) +
    scale_fill_manual(values = c("#80AD3C")) +
    custom_theme(diag_text = TRUE, hide_legend = TRUE)
dev.off()

png(filename = "results/images/Figure_4/F4_TBR1_barplot.png", width = 1600, height = 1400, res = 250)
ggplot(grouped_df, aes(x = CRISPR_type, y = TBR1_mean, fill = type)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = TBR1_mean - TBR1_sd, ymax = TBR1_mean + TBR1_sd), width = 0.2) +
    ylim(-1, 7) +
    scale_fill_manual(values = c("#80AD3C")) +
    custom_theme(diag_text = TRUE, hide_legend = TRUE)
dev.off()

png(filename = "results/images/Figure_4/F4_EPHB1_barplot.png", width = 1600, height = 1400, res = 250)
ggplot(grouped_df, aes(x = CRISPR_type, y = EPHB1_mean, fill = type)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = EPHB1_mean - EPHB1_sd, ymax = EPHB1_mean + EPHB1_sd), width = 0.2) +
    ylim(-1, 7) +
    scale_fill_manual(values = c("#80AD3C")) +
    custom_theme(diag_text = TRUE, hide_legend = TRUE)
dev.off()

png(filename = "results/images/Figure_4/F4_TRIM9_barplot.png", width = 1600, height = 1400, res = 250)
ggplot(grouped_df, aes(x = CRISPR_type, y = TRIM9_mean, fill = type)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = TRIM9_mean - TRIM9_sd, ymax = TRIM9_mean + TRIM9_sd), width = 0.2) +
    ylim(-1, 7) +
    scale_fill_manual(values = c("#80AD3C")) +
    custom_theme(diag_text = TRUE, hide_legend = TRUE)
dev.off()

png(filename = "results/images/Figure_4/F4_EPHA4_barplot.png", width = 1600, height = 1400, res = 250)
ggplot(grouped_df, aes(x = CRISPR_type, y = EPHA4_mean, fill = type)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = EPHA4_mean - EPHA4_sd, ymax = EPHA4_mean + EPHA4_sd), width = 0.2) +
    ylim(-1, 7) +
    scale_fill_manual(values = c("#80AD3C")) +
    custom_theme(diag_text = TRUE, hide_legend = TRUE)
dev.off()

png(filename = "results/images/Figure_4/F4_ZFHX4_barplot.png", width = 1600, height = 1400, res = 250)
ggplot(grouped_df, aes(x = CRISPR_type, y = ZFHX4_mean, fill = type)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = ZFHX4_mean - ZFHX4_sd, ymax = ZFHX4_mean + ZFHX4_sd), width = 0.2) +
    ylim(-1, 7) +
    scale_fill_manual(values = c("#80AD3C")) +
    custom_theme(diag_text = TRUE, hide_legend = TRUE)
dev.off()

# Getting genes within SHH WGCNA module
SHH_cluster_genes_df <- read.csv("results/tables/Figure_3/SHH_cluster.csv")
View(SHH_cluster_genes_df)

corr_pos <- SHH_cluster_genes_df[which(SHH_cluster_genes_df$cor > 0), ]$gene
corr_neg <- SHH_cluster_genes_df[which(SHH_cluster_genes_df$cor < 0), ]$gene

# Looking at SHH co-expressed genes expression in ventral CRISPR samples
norm_scale <- norm - min(norm)

# Positively correlated genes
meta_pos <- meta
maximum_exp <- max(norm_scale[rownames(norm_scale) %in% union(corr_pos, corr_neg), ])
meta_pos$CRISPR_type <- paste(meta_pos$CRISPR, meta_pos$type, sep = "_")
for (gene in rownames(norm_scale)[rownames(norm_scale) %in% corr_pos]) {
    meta_pos$gene <- norm_scale[gene, ]
    grp_df_pos <- meta_pos %>%
        group_by(CRISPR) %>%
        summarize(
            gene_mean = mean(gene),
            gene_sd = sd(gene)
        )
    grp_df_pos$CRISPR_type <- factor(c("vAN +/+", "vAN +/-", "vAN -/-"), levels = c("vAN +/+", "vAN +/-", "vAN -/-"))
    ggplot(grp_df_pos, aes(x = CRISPR_type, y = gene_mean, fill = CRISPR)) +
        geom_bar(stat = "identity") +
        geom_errorbar(aes(ymin = gene_mean - gene_sd, ymax = gene_mean + gene_sd), width = 0.2) +
        ylim(-1, maximum_exp) +
        scale_fill_manual(values = c("#80AD3C", "#b9e27b", "#dfe981")) +
        custom_theme(diag_text = TRUE, hide_legend = TRUE) +
        ggtitle(paste0("Scaled normalized expression of ", gene))
    ggsave(filename = paste0("results/images/Figure_4/vAN_coex_genes_barplot/positive_correlation/barplot_", gene, ".png"), units = "px", width = 1800, height = 1400, dpi = 250)
}

# Negatively correlated genes
meta_neg <- meta
meta_neg$CRISPR_type <- paste(meta_neg$CRISPR, meta_neg$type, sep = "_")
for (gene in rownames(norm_scale)[rownames(norm_scale) %in% corr_neg]) {
    meta_neg$gene <- norm_scale[gene, ]
    grp_df_neg <- meta_neg %>%
        group_by(CRISPR) %>%
        summarize(
            gene_mean = mean(gene),
            gene_sd = sd(gene)
        )
    grp_df_neg$CRISPR_type <- factor(c("vAN +/+", "vAN +/-", "vAN -/-"), levels = c("vAN +/+", "vAN +/-", "vAN -/-"))
    ggplot(grp_df_neg, aes(x = CRISPR_type, y = gene_mean, fill = CRISPR)) +
        geom_bar(stat = "identity") +
        geom_errorbar(aes(ymin = gene_mean - gene_sd, ymax = gene_mean + gene_sd), width = 0.2) +
        ylim(-1, maximum_exp) +
        scale_fill_manual(values = c("#80AD3C", "#b9e27b", "#dfe981")) +
        custom_theme(diag_text = TRUE, hide_legend = TRUE) +
        ggtitle(paste0("Scaled normalized expression of ", gene))
    ggsave(filename = paste0("results/images/Figure_4/vAN_coex_genes_barplot/negative_correlation/barplot_", gene, ".png"), units = "px", width = 1800, height = 1400, dpi = 250)
}
