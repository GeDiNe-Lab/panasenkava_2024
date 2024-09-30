# Loading packages and functions
library(Matrix)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(VennDiagram)
library(gridExtra)
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
DE_vAN_het_vs_control <- dds_vAN %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("CRISPR", "hetero", "control")) %>%
    as.data.frame() %>%
    na.omit()
DE_vAN_het_vs_control$gene <- gene_converter(rownames(DE_vAN_het_vs_control), "ENSEMBL", "SYMBOL")
DE_vAN_het_vs_control$Thr1_5 <- DE_vAN_het_vs_control$padj < 0.05 & abs(DE_vAN_het_vs_control$log2FoldChange) >= 1.5
DE_vAN_het_vs_control$FCsign <- ifelse(DE_vAN_het_vs_control$log2FoldChange < 0, "neg", "pos")
DE_vAN_het_vs_control_f <- filter(DE_vAN_het_vs_control, !is.na(gene))
DE_vAN_het_vs_control_f <- filter(DE_vAN_het_vs_control_f, padj < 0.01, abs(log2FoldChange) >= 2)

# DEGs between control and homozygous for ventral, dorsal, and ventral+dorsal samples
DE_vAN_homo_vs_control <- dds_vAN %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("CRISPR", "homo", "control")) %>%
    as.data.frame() %>%
    na.omit()
DE_vAN_homo_vs_control$gene <- gene_converter(rownames(DE_vAN_homo_vs_control), "ENSEMBL", "SYMBOL") %>% as.vector()
DE_vAN_homo_vs_control$Thr1_5 <- DE_vAN_homo_vs_control$padj < 0.05 & abs(DE_vAN_homo_vs_control$log2FoldChange) >= 1.5
DE_vAN_homo_vs_control$FCsign <- ifelse(DE_vAN_homo_vs_control$log2FoldChange < 0, "neg", "pos")
DE_vAN_homo_vs_control_f <- filter(DE_vAN_homo_vs_control, !is.na(gene))
DE_vAN_homo_vs_control_f <- filter(DE_vAN_homo_vs_control_f, padj < 0.01, abs(log2FoldChange) >= 2)

# DEGs between heterozygous and homozygous for ventral, dorsal, and ventral+dorsal samples
DE_vAN_homo_vs_het <- dds_vAN %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("CRISPR", "homo", "hetero")) %>%
    as.data.frame() %>%
    na.omit()
DE_vAN_homo_vs_het$gene <- gene_converter(rownames(DE_vAN_homo_vs_het), "ENSEMBL", "SYMBOL")
DE_vAN_homo_vs_het$Thr1_5 <- DE_vAN_homo_vs_het$padj < 0.05 & abs(DE_vAN_homo_vs_het$log2FoldChange) >= 1.5
DE_vAN_homo_vs_het$FCsign <- ifelse(DE_vAN_homo_vs_het$log2FoldChange < 0, "neg", "pos")
DE_vAN_homo_vs_het_f <- filter(DE_vAN_homo_vs_het, !is.na(gene))
DE_vAN_homo_vs_het_f <- filter(DE_vAN_homo_vs_het_f, padj < 0.01, abs(log2FoldChange) >= 2)

#  Saving DE tables
write.csv(DE_vAN_het_vs_control, "results/tables/Figure_4/DE_vAN_het_vs_control.csv")
write.csv(DE_vAN_homo_vs_control, "results/tables/Figure_4/DE_vAN_homo_vs_control.csv")
write.csv(DE_vAN_homo_vs_het, "results/tables/Figure_4/DE_vAN_homo_vs_het.csv")

# Making volcanp plots
DE_CRISPR <- list(
    vAN_het_vs_control = DE_vAN_het_vs_control_f,
    vAN_homo_vs_control = DE_vAN_homo_vs_control_f,
    vAN_homo_vs_het = DE_vAN_homo_vs_het_f
)

for (contrast in names(DE_CRISPR)) {
    print(contrast)
    DE <- DE_CRISPR[[contrast]]
    print(max(-log10(DE$padj) + (max(-log10(DE$padj)) * 0.1)))
    ggplot(DE, aes(x = log2FoldChange, y = -log10(padj), label = gene)) +
        geom_text(size = 2) +
        custom_theme() +
        geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
        geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
        labs(x = "log2FoldChange", y = "-log10(padj)", title = paste0("DE: ", contrast), subtitle = "|log2FC| >= 2 & FDR < 0.01")
    ggsave(filename = paste0("results/images/Figure_4/Volcano_", contrast, ".png"), units = "px", width = 1800, height = 1400, dpi = 250)
}
View(DE_vAN_homo_vs_het_f)
filter(DE_vAN_het_vs_control_f, padj < 1.6857e-149)$gene
filter(DE_vAN_homo_vs_control_f, padj < 1.6703e-274)$gene
filter(DE_vAN_homo_vs_het_f, padj < 1.1168e-304)$gene




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



SHH_cluster_genes_df <- read.csv("results/tables/Figure_3/SHH_cluster.csv")
View(SHH_cluster_genes_df)

meta_CRISPR <- filter(rawmeta, diff == "diff12", type == "ventral")
meta_cyclo <- filter(rawmeta, diff == "diff9", type != "dorsal", sample != "L9C1_2")

c("EPHB1", "TRIM9", "EPHA4", "ZFHX4") %>% gene_converter("SYMBOL", "ENSEMBL")

markers <- c(SHH_cluster_genes_df$ENSEMBLE, c("EPHB1", "TRIM9", "EPHA4", "ZFHX4") %>% gene_converter("SYMBOL", "ENSEMBL"))
# Make sure to keep only markers in the matrix
markers <- intersect(markers, rownames(rawcounts))

# filtering out lowly expressed genes and keeping only seleted samples
retained_row_CRISPR <- rawcounts[rownames(rawcounts) %in% markers, meta_CRISPR$sample]
counts_CRISPR <- rawcounts[, meta_CRISPR$sample][rowSums(rawcounts[, meta_CRISPR$sample]) >= 25, ]

retained_row_cyclo <- rawcounts[rownames(rawcounts) %in% markers, meta_cyclo$sample]
counts_cyclo <- rawcounts[, meta_cyclo$sample][rowSums(rawcounts[, meta_cyclo$sample]) >= 25, ]

# putting back potentially filtered out markers (posterior markers for example as they should not be expressed)
if (length(which(!rownames(retained_row_CRISPR) %in% rownames(counts_CRISPR))) == 1) {
    counts_CRISPR <- rbind(counts_CRISPR, retained_row_CRISPR[which(!rownames(retained_row_CRISPR) %in% rownames(counts_CRISPR)), ])
    rownames(counts_CRISPR)[nrow(counts_CRISPR)] <- rownames(retained_row_CRISPR)[which(!rownames(retained_row_CRISPR) %in% rownames(counts_CRISPR))]
} else if (length(which(!rownames(retained_row_CRISPR) %in% rownames(counts_CRISPR))) > 1) {
    counts_CRISPR <- rbind(counts_CRISPR, retained_row_CRISPR[which(!rownames(retained_row_CRISPR) %in% rownames(counts_CRISPR)), ])
}

if (length(which(!rownames(retained_row_cyclo) %in% rownames(counts_cyclo))) == 1) {
    counts_cyclo <- rbind(counts_cyclo, retained_row_cyclo[which(!rownames(retained_row_cyclo) %in% rownames(counts_cyclo)), ])
    rownames(counts_cyclo)[nrow(counts_cyclo)] <- rownames(retained_row_cyclo)[which(!rownames(retained_row_cyclo) %in% rownames(counts_cyclo))]
} else if (length(which(!rownames(retained_row_cyclo) %in% rownames(counts_cyclo))) > 1) {
    counts_cyclo <- rbind(counts_cyclo, retained_row_cyclo[which(!rownames(retained_row_cyclo) %in% rownames(counts_cyclo)), ])
}


dds_cyclo <- DESeqDataSetFromMatrix(
    countData = counts_cyclo,
    colData = meta_cyclo,
    design = ~cyclo_dose_qual
)
vsd_cyclo <- vst(dds_cyclo, blind = FALSE)
dds_CRISPR <- DESeqDataSetFromMatrix(
    countData = counts_CRISPR,
    colData = meta_CRISPR,
    design = ~CRISPR
)
vsd_CRISPR <- vst(dds_CRISPR, blind = FALSE)

# Looking at SHH co-expressed genes expression in ventral CRISPR samples
scale_cyclo <- assay(vsd_cyclo) - min(assay(vsd_cyclo))
rownames(scale_cyclo) <- gene_converter(rownames(scale_cyclo), "ENSEMBL", "SYMBOL")
scale_CRISPR <- assay(vsd_CRISPR) - min(assay(vsd_CRISPR))
rownames(scale_CRISPR) <- gene_converter(rownames(scale_CRISPR), "ENSEMBL", "SYMBOL")

meta_bp_CRISPR <- meta_CRISPR
meta_bp_CRISPR$CRISPR_type <- paste(meta_bp_CRISPR$CRISPR, meta_bp_CRISPR$type, sep = "_")

meta_bp_cyclo <- meta_cyclo

for (gene in SHH_cluster_genes_df$gene) {
    meta_bp_CRISPR$gene <- scale_CRISPR[gene, ]
    red_df_CRISPR <- dplyr::select(meta_bp_CRISPR, c("CRISPR", "gene"))
    stat_df_CRISPR <- data.frame(control = filter(red_df_CRISPR, CRISPR == "control")$gene, hetero = filter(red_df_CRISPR, CRISPR == "hetero")$gene, homo = filter(red_df_CRISPR, CRISPR == "homo")$gene)

    grp_df_CRISPR <- meta_bp_CRISPR %>%
        group_by(CRISPR) %>%
        summarize(
            gene_mean = mean(gene),
            gene_sd = sd(gene)
        )

    grp_df_CRISPR$CRISPR <- factor(c("+/+", "+/-", "-/-"), levels = c("+/+", "+/-", "-/-"))
    c_het_p <- wilcox.test(stat_df_CRISPR$control, stat_df_CRISPR$hetero)$p.value %>% round(4)
    c_hom_p <- wilcox.test(stat_df_CRISPR$control, stat_df_CRISPR$homo)$p.value %>% round(4)
    het_hom_p <- wilcox.test(stat_df_CRISPR$hetero, stat_df_CRISPR$homo)$p.value %>% round(4)
    signif <- c(c_het_p, c_hom_p, het_hom_p) %>% sapply(function(pval) {
        if (is.nan(pval)) {
            return("***")
        } else if (pval < 0.001) {
            return("***")
        } else if (pval < 0.01) {
            return("**")
        } else if (pval < 0.05) {
            return("*")
        } else {
            return(" ")
        }
    })
    pval_data <- data.frame(
        group1 = c("+/+", "+/+", "+/-"),
        group2 = c("+/-", "-/-", "-/-"),
        p_value = paste(c(c_het_p, c_hom_p, het_hom_p), signif, sep = " "),
        y_position = c(max(scale_CRISPR) + max(scale_CRISPR) * 0.1, max(scale_CRISPR) + max(scale_CRISPR) * 0.2, max(scale_CRISPR) + max(scale_CRISPR) * 0.3) # Adjust these based on your plot's scale
    )
    bp_CRISPR <- ggbarplot(grp_df_CRISPR,
        x = "CRISPR", y = "gene_mean",
        fill = "CRISPR", palette = c("#80AD3C", "#b9e27b", "#dfe981"),
        add = "gene_sd",
        width = 0.9,
        ggtheme = custom_theme(diag_text = TRUE, hide_legend = TRUE)
    ) +
        ylim(-0.2, max(scale_CRISPR) + max(scale_CRISPR) * 0.3) +
        geom_errorbar(data = grp_df_CRISPR, aes(ymin = gene_mean - gene_sd, ymax = gene_mean + gene_sd), width = 0.2) +
        stat_pvalue_manual(
            pval_data,
            y.position = "y_position",
            label = "p_value"
        )


    meta_bp_cyclo$gene <- scale_cyclo[gene, ]
    red_df_cyclo <- dplyr::select(meta_bp_cyclo, c("cyclo_dose_qual", "gene"))

    grp_df_cyclo <- meta_bp_cyclo %>%
        group_by(cyclo_dose_qual) %>%
        summarize(
            gene_mean = mean(gene),
            gene_sd = sd(gene)
        )
    grp_df_cyclo <- grp_df_cyclo[c(3, 2, 1), ] # making sure everything is in the right order
    grp_df_cyclo$cyclo_dose_qual <- factor(c("none", "low", "high"), levels = c("none", "low", "high"))
    no_low_p <- wilcox.test(filter(red_df_cyclo, cyclo_dose_qual == "none")$gene, filter(red_df_cyclo, cyclo_dose_qual == "low")$gene)$p.value %>% round(4)
    no_high_p <- wilcox.test(filter(red_df_cyclo, cyclo_dose_qual == "none")$gene, filter(red_df_cyclo, cyclo_dose_qual == "high")$gene)$p.value %>% round(4)
    low_high_p <- wilcox.test(filter(red_df_cyclo, cyclo_dose_qual == "low")$gene, filter(red_df_cyclo, cyclo_dose_qual == "high")$gene)$p.value %>% round(4)
    signif <- c(no_low_p, no_high_p, low_high_p) %>% sapply(function(pval) {
        if (is.nan(pval)) {
            return("")
        } else if (pval < 0.001) {
            return("***")
        } else if (pval < 0.01) {
            return("**")
        } else if (pval < 0.05) {
            return("*")
        } else {
            return("")
        }
    })
    pval_data <- data.frame(
        group1 = c("none", "none", "low"),
        group2 = c("low", "high", "high"),
        p_value = paste(c(no_low_p, no_high_p, low_high_p), signif, sep = " "),
        y_position = c(max(scale_cyclo) + max(scale_cyclo) * 0.1, max(scale_cyclo) + max(scale_cyclo) * 0.2, max(scale_cyclo) + max(scale_cyclo) * 0.3) # Adjust these based on your plot's scale
    )
    bp_cyclo <- ggbarplot(grp_df_cyclo,
        x = "cyclo_dose_qual", y = "gene_mean",
        fill = "cyclo_dose_qual", palette = c("#80AD3C", "#f6a563", "#f86110"),
        add = "gene_sd",
        width = 0.9,
        ggtheme = custom_theme(diag_text = TRUE, hide_legend = TRUE)
    ) +
        ylim(-0.2, max(scale_cyclo) + max(scale_cyclo) * 0.3) +
        geom_errorbar(data = grp_df_cyclo, aes(ymin = gene_mean - gene_sd, ymax = gene_mean + gene_sd), width = 0.2) +
        stat_pvalue_manual(
            pval_data,
            y.position = "y_position",
            label = "p_value"
        )
    title <- textGrob(paste0("Scaled normalized expression of ", gene), gp = gpar(fontsize = 16, fontface = "bold"))

    # Arrange the two plots and add the title
    finalplot <- grid.arrange(
        title, # Add the title
        arrangeGrob(bp_CRISPR, bp_cyclo, ncol = 2), # Add the plots side by side
        nrow = 2, # The layout will have two rows, one for title and one for the plots
        heights = c(0.1, 1) # Relative heights (title takes 10% of the height)
    )
    if (SHH_cluster_genes_df[which(SHH_cluster_genes_df$gene == gene), ]$cor > 0) {
        ggsave(filename = paste0("results/images/Figure_4/vAN_coex_genes_barplot/positive_correlation/barplot_", gene, ".png"), plot = finalplot, units = "px", width = 1800, height = 1400, dpi = 250)
    } else {
        ggsave(filename = paste0("results/images/Figure_4/vAN_coex_genes_barplot/negative_correlation/barplot_", gene, ".png"), plot = finalplot, units = "px", width = 1800, height = 1400, dpi = 250)
    }
}


meta_CRISPRcyclo <- filter(rawmeta, diff == "diff12" & ((type != "dorsal" & CRISPR == "control") | (type == "ventral" & CRISPR != "control")))

dds_CRISPRcyclo <- DESeqDataSetFromMatrix(
    rawcounts[, meta_CRISPRcyclo$sample][rowSums(rawcounts[, meta_CRISPRcyclo$sample]) >= 25, ],
    colData = meta_CRISPRcyclo,
    design = ~CRISPR
)
vsd_CRISPRcyclo <- vst(dds_CRISPRcyclo, blind = FALSE)

max(assay(vsd_cyclo))
distance <- dist(t(assay(vsd_CRISPRcyclo)))
distance

# Step 2: Convert the distance matrix to a square matrix
dist_matrix <- as.matrix(distance)

# Step 3: Create a heatmap
png(filename = "results/images/Figure_4/CRISPR_cyclo_distances.png", width = 2000, height = 2000, res = 250)
heatmap(dist_matrix,
    main = "Distance Matrix Heatmap",
    xlab = "Samples",
    ylab = "Samples",
    col = heat.colors(256)
)
dev.off()




scaled_mat <- t(apply(assay(vsd_vAN)[SHH_cluster_genes_df$ENSEMBLE %in% rownames(assay(vsd_vAN)), ], 1, scale))
colnames(scaled_mat) <- colnames(assay(vsd_vAN))

# hierarchical clustering using euclidian distance and "complete" method
clustering <- hclust(dist(scaled_mat))
clusters <- cutree(clustering, k = 4)

# Subclustering of each cluster
sub_clusters_list <- unique(clusters) %>% lapply(function(cluster) {
    sub_mat <- scaled_mat[names(clusters[which(clusters == cluster)]), ]
    sub_clustering <- hclust(dist(sub_mat))
    return(cutree(sub_clustering, k = 5))
})
names(sub_clusters_list) <- paste0("cluster_", unique(clusters))
sub_clusters <- sub_clusters_list %>%
    unname() %>%
    unlist()
sub_clusters <- sub_clusters[names(clusters)]

sub_sub_clusters_list <- unique(sub_clusters) %>% lapply(function(sub_cluster) {
    sub_mat <- scaled_mat[names(sub_clusters[which(sub_clusters == sub_cluster)]), ]
    sub_sub_clustering <- hclust(dist(sub_mat))
    return(cutree(sub_sub_clustering, k = 5))
})
names(sub_sub_clusters_list) <- paste0("sub_cluster_", unique(sub_clusters))
sub_sub_clusters <- sub_sub_clusters_list %>%
    unname() %>%
    unlist()
sub_sub_clusters <- sub_sub_clusters[names(sub_clusters)]

# Heatmap gene annotation
clusters_ha <- rowAnnotation(
    cluster = as.character(clusters[clustering$order]),
    sub_cluster = as.character(sub_clusters[clustering$order]),
    sub_sub_cluster = as.character(sub_sub_clusters[clustering$order]),
    col = list(
        cluster = c(
            "1" = "#b16060",
            "2" = "#4d6da5",
            "3" = "#2f7439",
            "4" = "#613999"
        ),
        sub_cluster = c(
            "1" = "black",
            "2" = "pink",
            "3" = "yellow",
            "4" = "brown",
            "5" = "grey"
        ),
        sub_sub_cluster = c(
            "1" = "black",
            "2" = "pink",
            "3" = "yellow",
            "4" = "brown",
            "5" = "grey"
        )
    )
)



png(filename = "results/images/Figure_4/CRISPR_cycloWGCNA_genes_HM.png", width = 2400, height = 1600, res = 250)
Heatmap(
    scaled_mat[clustering$order, ],
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

cyclo_genes_df <- read.csv("results/tables/Figure_3/cyclo_genes_df.csv")
cyclo_genes_df %>% View()

cyclo_genes_df$CRISPR_cluster <- cyclo_genes_df$X %>% sapply(function(gene) {
    if (gene %in% names(clusters)) {
        return(clusters[gene])
    } else {
        return("NA")
    }
})

cyclo_genes_df$CRISPR_sub_cluster <- cyclo_genes_df$X %>% sapply(function(gene) {
    if (gene %in% names(sub_clusters)) {
        return(sub_clusters[gene])
    } else {
        return("NA")
    }
})

cyclo_genes_df$CRISPR_sub_sub_cluster <- cyclo_genes_df$X %>% sapply(function(gene) {
    if (gene %in% names(sub_sub_clusters)) {
        return(sub_sub_clusters[gene])
    } else {
        return("NA")
    }
})
View(cyclo_genes_df)
write.csv(cyclo_genes_df, "results/tables/Figure_4/cyclo_genes_df_CRISPR_cluster.csv")
