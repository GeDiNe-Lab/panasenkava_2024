# Loading packages and functions
library(Matrix)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(VennDiagram)
library(tibble)

rstudioapi::getSourceEditorContext()$path %>%
    str_split("/") %>%
    unlist() %>%
    head(-3) %>%
    str_c(collapse = "/") %>%
    str_c("/") %>%
    setwd()

source("R/custom_fct.R")
rawcounts <- readcounts("/home/jules/Documents/phd/Data/lab_RNAseq/IPSdiff12/IPSdiff12_counts.csv")
rawmeta <- read.table("/home/jules/Documents/phd/Data/lab_RNAseq/IPSdiff12/IPSdiff12_metadata.csv", sep = ",", header = T)
View(meta)
meta <- filter(rawmeta, type %in% c("cyclo", "ventral"))
counts <- rawcounts[which(rowSums(rawcounts) >= 50), meta$SAMPLE_NAME]
counts %>% colnames()
meta$cyclo_dose_qual <- meta$cyclo_dose %>% sapply(function(x) {
    if (x %in% c(0.125, 0.25)) {
        return("low")
    } else if (x %in% c(0.5, 1)) {
        return("high")
    } else {
        return("ventral")
    }
})
meta$cyclo_dose <- as.factor(meta$cyclo_dose)
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ cyclo_dose_qual + CRISPR
)

vsd <- varianceStabilizingTransformation(dds)

# PCA plot
pca.data <- plotPCA(vsd, intgroup = c("cyclo_dose", "cyclo_dose_qual", "CRISPR"), returnData = TRUE)
percentVar <- round(100 * attr(pca.data, "percentVar"))

png(filename = "results/images/Figure_4/F4_PCA_type.png", width = 1600, height = 1200, res = 250)
ggplot(pca.data, aes(PC1, PC2, color = cyclo_dose, shape = CRISPR)) +
    geom_point(size = 3, stroke = 1) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_color_manual(values = c("#80AD3C", "#dba800", "#fe8700", "#d55500", "#612000")) +
    custom_theme() +
    ggtitle("PCA of ventral and cyclo sample in CRISPR line")
dev.off()

dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~CRISPR
)

DE_control_vs_het <- dds %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("CRISPR", "hetero", "control")) %>%
    as.data.frame() %>%
    na.omit()
DE_control_vs_het_f <- filter(DE_control_vs_het, padj < 0.01, abs(log2FoldChange) >= 1)
DE_control_vs_het_f$gene <- gene_converter(rownames(DE_control_vs_het_f), "ENSEMBL", "SYMBOL")
DE_control_vs_het_f <- filter(DE_control_vs_het_f, !is.na(gene))

DE_control_vs_homo <- dds %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("CRISPR", "homo", "control")) %>%
    as.data.frame() %>%
    na.omit()
DE_control_vs_homo_f <- filter(DE_control_vs_homo, padj < 0.01, abs(log2FoldChange) >= 1)
DE_control_vs_homo_f$gene <- gene_converter(rownames(DE_control_vs_homo_f), "ENSEMBL", "SYMBOL")
DE_control_vs_homo_f <- filter(DE_control_vs_homo_f, !is.na(gene))

DE_het_vs_homo <- dds %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("CRISPR", "homo", "hetero")) %>%
    as.data.frame() %>%
    na.omit()
DE_het_vs_homo_f <- filter(DE_het_vs_homo, padj < 0.01, abs(log2FoldChange) >= 1)
DE_het_vs_homo_f$gene <- gene_converter(rownames(DE_het_vs_homo_f), "ENSEMBL", "SYMBOL")
DE_het_vs_homo_f <- filter(DE_het_vs_homo_f, !is.na(gene))

write.csv(DE_control_vs_het_f, "results/tables/Figure_4/DE_control_vs_het.csv")
write.csv(DE_control_vs_homo_f, "results/tables/Figure_4/DE_control_vs_homo.csv")
write.csv(DE_het_vs_homo_f, "results/tables/Figure_4/DE_het_vs_homo.csv")

DE_CRISPR <- list(
    C_vs_Het = DE_control_vs_het_f,
    C_vs_Hom = DE_control_vs_homo_f,
    Het_vs_Hom = DE_het_vs_homo_f
)
names(DE_CRISPR)

for (contrast in names(DE_CRISPR)) {
    print(contrast)
    DE <- DE_CRISPR[[contrast]]
    ggplot(DE, aes(x = log2FoldChange, y = -log10(padj), label = gene)) +
        ggrepel::geom_text_repel(box.padding = 0.001, size = 2.5, max.overlaps = 20) +
        custom_theme() +
        geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
        labs(x = "log2FoldChange", y = "-log10(padj)", title = paste0("DE: ", contrast))
    ggsave(filename = paste0("results/images/Figure_4/Volcano_", contrast, ".png"), units = "px", width = 1800, height = 1400, dpi = 250)
}

dds_ventral <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "ventral")$SAMPLE_NAME],
    colData = filter(meta, type == "ventral"),
    design = ~CRISPR
)

DE_vAN_control_vs_het <- dds_ventral %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("CRISPR", "hetero", "control")) %>%
    as.data.frame() %>%
    na.omit()
DE_vAN_control_vs_het_f <- filter(DE_vAN_control_vs_het, padj < 0.01, abs(log2FoldChange) >= 1)
DE_vAN_control_vs_het_f$gene <- gene_converter(rownames(DE_vAN_control_vs_het_f), "ENSEMBL", "SYMBOL")
DE_vAN_control_vs_het_f <- filter(DE_vAN_control_vs_het_f, !is.na(gene))

DE_vAN_control_vs_homo <- dds_ventral %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("CRISPR", "homo", "control")) %>%
    as.data.frame() %>%
    na.omit()
DE_vAN_control_vs_homo_f <- filter(DE_vAN_control_vs_homo, padj < 0.01, abs(log2FoldChange) >= 1)
DE_vAN_control_vs_homo_f$gene <- gene_converter(rownames(DE_vAN_control_vs_homo_f), "ENSEMBL", "SYMBOL")
DE_vAN_control_vs_homo_f <- filter(DE_vAN_control_vs_homo_f, !is.na(gene))

DE_vAN_het_vs_homo <- dds_ventral %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("CRISPR", "homo", "hetero")) %>%
    as.data.frame() %>%
    na.omit()
DE_vAN_het_vs_homo_f <- filter(DE_vAN_het_vs_homo, padj < 0.01, abs(log2FoldChange) >= 1)
DE_vAN_het_vs_homo_f$gene <- gene_converter(rownames(DE_vAN_het_vs_homo_f), "ENSEMBL", "SYMBOL")
DE_vAN_het_vs_homo_f <- filter(DE_vAN_het_vs_homo_f, !is.na(gene))

write.csv(DE_vAN_control_vs_het_f, "results/tables/Figure_4/DE_vAN_control_vs_het.csv")
write.csv(DE_vAN_control_vs_homo_f, "results/tables/Figure_4/DE_vAN_control_vs_homo.csv")
write.csv(DE_vAN_het_vs_homo_f, "results/tables/Figure_4/DE_vAN_het_vs_homo.csv")

DE_vAN_CRISPR <- list(
    C_vs_Het = DE_vAN_control_vs_het_f,
    C_vs_Hom = DE_vAN_control_vs_homo_f,
    Het_vs_Hom = DE_vAN_het_vs_homo_f
)
names(DE_vAN_CRISPR)

for (contrast in names(DE_vAN_CRISPR)) {
    print(contrast)
    DE <- DE_vAN_CRISPR[[contrast]]
    ggplot(DE, aes(x = log2FoldChange, y = -log10(padj), label = gene)) +
        ggrepel::geom_text_repel(box.padding = 0.001, size = 2.5, max.overlaps = 20) +
        custom_theme() +
        geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
        labs(x = "log2FoldChange", y = "-log10(padj)", title = paste0("DE vAN: ", contrast))
    ggsave(filename = paste0("results/images/Figure_4/Volcano_vAN_", contrast, ".png"), units = "px", width = 1800, height = 1400, dpi = 250)
}

dds_low_cyclo <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, cyclo_dose_qual == "low")$SAMPLE_NAME],
    colData = filter(meta, cyclo_dose_qual == "low"),
    design = ~CRISPR
)

DE_lowC_control_vs_het <- dds_low_cyclo %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("CRISPR", "hetero", "control")) %>%
    as.data.frame() %>%
    na.omit()
DE_lowC_control_vs_het_f <- filter(DE_lowC_control_vs_het, padj < 0.01, abs(log2FoldChange) >= 1)
DE_lowC_control_vs_het_f$gene <- gene_converter(rownames(DE_lowC_control_vs_het_f), "ENSEMBL", "SYMBOL")
DE_lowC_control_vs_het_f <- filter(DE_lowC_control_vs_het_f, !is.na(gene))

DE_lowC_control_vs_homo <- dds_low_cyclo %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("CRISPR", "homo", "control")) %>%
    as.data.frame() %>%
    na.omit()
DE_lowC_control_vs_homo_f <- filter(DE_lowC_control_vs_homo, padj < 0.01, abs(log2FoldChange) >= 1)
DE_lowC_control_vs_homo_f$gene <- gene_converter(rownames(DE_lowC_control_vs_homo_f), "ENSEMBL", "SYMBOL")
DE_lowC_control_vs_homo_f <- filter(DE_lowC_control_vs_homo_f, !is.na(gene))

DE_lowC_het_vs_homo <- dds_low_cyclo %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("CRISPR", "homo", "hetero")) %>%
    as.data.frame() %>%
    na.omit()
DE_lowC_het_vs_homo_f <- filter(DE_lowC_het_vs_homo, padj < 0.01, abs(log2FoldChange) >= 1)
DE_lowC_het_vs_homo_f$gene <- gene_converter(rownames(DE_lowC_het_vs_homo_f), "ENSEMBL", "SYMBOL")
DE_lowC_het_vs_homo_f <- filter(DE_lowC_het_vs_homo_f, !is.na(gene))

write.csv(DE_lowC_control_vs_het_f, "results/tables/Figure_4/DE_lowC_control_vs_het.csv")
write.csv(DE_lowC_control_vs_homo_f, "results/tables/Figure_4/DE_lowC_control_vs_homo.csv")
write.csv(DE_lowC_het_vs_homo_f, "results/tables/Figure_4/DE_lowC_het_vs_homo.csv")

DE_lowC_CRISPR <- list(
    C_vs_Het = DE_lowC_control_vs_het_f,
    C_vs_Hom = DE_lowC_control_vs_homo_f,
    Het_vs_Hom = DE_lowC_het_vs_homo_f
)
names(DE_lowC_CRISPR)

for (contrast in names(DE_lowC_CRISPR)) {
    print(contrast)
    DE <- DE_lowC_CRISPR[[contrast]]
    ggplot(DE, aes(x = log2FoldChange, y = -log10(padj), label = gene)) +
        ggrepel::geom_text_repel(box.padding = 0.001, size = 2.5, max.overlaps = 20) +
        custom_theme() +
        geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
        labs(x = "log2FoldChange", y = "-log10(padj)", title = paste0("DE lowC: ", contrast))
    ggsave(filename = paste0("results/images/Figure_4/Volcano_lowC_", contrast, ".png"), units = "px", width = 1800, height = 1400, dpi = 250)
}

dds_high_cyclo <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, cyclo_dose_qual == "high")$SAMPLE_NAME],
    colData = filter(meta, cyclo_dose_qual == "high"),
    design = ~CRISPR
)

DE_highC_control_vs_het <- dds_high_cyclo %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("CRISPR", "hetero", "control")) %>%
    as.data.frame() %>%
    na.omit()
DE_highC_control_vs_het_f <- filter(DE_highC_control_vs_het, padj < 0.01, abs(log2FoldChange) >= 1)
DE_highC_control_vs_het_f$gene <- gene_converter(rownames(DE_highC_control_vs_het_f), "ENSEMBL", "SYMBOL")
DE_highC_control_vs_het_f <- filter(DE_highC_control_vs_het_f, !is.na(gene))

DE_highC_control_vs_homo <- dds_high_cyclo %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("CRISPR", "homo", "control")) %>%
    as.data.frame() %>%
    na.omit()
DE_highC_control_vs_homo_f <- filter(DE_highC_control_vs_homo, padj < 0.01, abs(log2FoldChange) >= 1)
DE_highC_control_vs_homo_f$gene <- gene_converter(rownames(DE_highC_control_vs_homo_f), "ENSEMBL", "SYMBOL")
DE_highC_control_vs_homo_f <- filter(DE_highC_control_vs_homo_f, !is.na(gene))

DE_highC_het_vs_homo <- dds_high_cyclo %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("CRISPR", "homo", "hetero")) %>%
    as.data.frame() %>%
    na.omit()
DE_highC_het_vs_homo_f <- filter(DE_highC_het_vs_homo, padj < 0.01, abs(log2FoldChange) >= 1)
DE_highC_het_vs_homo_f$gene <- gene_converter(rownames(DE_highC_het_vs_homo_f), "ENSEMBL", "SYMBOL")
DE_highC_het_vs_homo_f <- filter(DE_highC_het_vs_homo_f, !is.na(gene))

write.csv(DE_highC_control_vs_het_f, "results/tables/Figure_4/DE_highC_control_vs_het.csv")
write.csv(DE_highC_control_vs_homo_f, "results/tables/Figure_4/DE_highC_control_vs_homo.csv")
write.csv(DE_highC_het_vs_homo_f, "results/tables/Figure_4/DE_highC_het_vs_homo.csv")

DE_highC_CRISPR <- list(
    C_vs_Het = DE_highC_control_vs_het_f,
    C_vs_Hom = DE_highC_control_vs_homo_f,
    Het_vs_Hom = DE_highC_het_vs_homo_f
)
names(DE_highC_CRISPR)

for (contrast in names(DE_highC_CRISPR)) {
    print(contrast)
    DE <- DE_highC_CRISPR[[contrast]]
    ggplot(DE, aes(x = log2FoldChange, y = -log10(padj), label = gene)) +
        ggrepel::geom_text_repel(box.padding = 0.001, size = 2.5, max.overlaps = 20) +
        custom_theme() +
        geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
        labs(x = "log2FoldChange", y = "-log10(padj)", title = paste0("DE highC: ", contrast))
    ggsave(filename = paste0("results/images/Figure_4/Volcano_highC_", contrast, ".png"), units = "px", width = 1800, height = 1400, dpi = 250)
}

meta <- filter(rawmeta, type %in% c("ventral", "dorsal"))
counts <- rawcounts[which(rowSums(rawcounts) >= 25), meta$SAMPLE_NAME]
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ CRISPR + type
)

vsd <- varianceStabilizingTransformation(dds)

norm <- assay(vsd)
rownames(norm) <- gene_converter(rownames(norm), "ENSEMBL", "SYMBOL")

meta$SHH <- norm["SHH", ] - min(norm)
meta$NKX21 <- norm["NKX2-1", ] - min(norm)
meta$PAX6 <- norm["PAX6", ] - min(norm)
meta$TBR1 <- norm["TBR1", ] - min(norm)
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
        TBR1_sd = sd(TBR1)
    )
grouped_df <- grouped_df[order(grouped_df$type, decreasing = TRUE), ]
grouped_df$CRISPR_type <- factor(c("vAN +/+", "vAN +/-", "vAN -/-", "dAN +/+", "dAN +/-", "dAN -/-"), levels = c("vAN +/+", "vAN +/-", "vAN -/-", "dAN +/+", "dAN +/-", "dAN -/-"))

png(filename = "results/images/Figure_4/F4_SHH_barplot.png", width = 1600, height = 1400, res = 250)
ggplot(grouped_df, aes(x = CRISPR_type, y = SHH_mean, fill = type)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = SHH_mean - SHH_sd, ymax = SHH_mean + SHH_sd), width = 0.2) +
    ylim(-1, 7) +
    scale_fill_manual(values = c("#A1A1DE", "#80AD3C")) +
    custom_theme(diag_text = TRUE)
dev.off()

png(filename = "results/images/Figure_4/F4_NKX21_barplot.png", width = 1600, height = 1400, res = 250)
ggplot(grouped_df, aes(x = CRISPR_type, y = NKX21_mean, fill = type)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = NKX21_mean - NKX21_sd, ymax = NKX21_mean + NKX21_sd), width = 0.2) +
    ylim(-1, 7) +
    scale_fill_manual(values = c("#A1A1DE", "#80AD3C")) +
    custom_theme(diag_text = TRUE)
dev.off()

png(filename = "results/images/Figure_4/F4_PAX6_barplot.png", width = 1600, height = 1400, res = 250)
ggplot(grouped_df, aes(x = CRISPR_type, y = PAX6_mean, fill = type)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = PAX6_mean - PAX6_sd, ymax = PAX6_mean + PAX6_sd), width = 0.2) +
    ylim(-1, 7) +
    scale_fill_manual(values = c("#A1A1DE", "#80AD3C")) +
    custom_theme(diag_text = TRUE)
dev.off()

png(filename = "results/images/Figure_4/F4_TBR1_barplot.png", width = 1600, height = 1400, res = 250)
ggplot(grouped_df, aes(x = CRISPR_type, y = TBR1_mean, fill = type)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = TBR1_mean - TBR1_sd, ymax = TBR1_mean + TBR1_sd), width = 0.2) +
    ylim(-1, 7) +
    scale_fill_manual(values = c("#A1A1DE", "#80AD3C")) +
    custom_theme(diag_text = TRUE)
dev.off()

corr_pos <- read.csv("results/tables/Figure_3/cytoscape_SHH_pos_0_80.csv")
corr_pos <- filter(corr_pos, target_info %in% c("selected", "known"))
corr_neg <- read.csv("results/tables/Figure_3/cytoscape_SHH_neg_0_80.csv")
corr_neg <- filter(corr_neg, target_info %in% c("selected", "known"))


meta <- filter(rawmeta, type %in% c("ventral"))
counts <- rawcounts[which(rowSums(rawcounts) >= 25), meta$SAMPLE_NAME]
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~CRISPR
)

vsd <- varianceStabilizingTransformation(dds)

norm <- assay(vsd)
norm <- norm - min(norm)
rownames(norm) <- gene_converter(rownames(norm), "ENSEMBL", "SYMBOL")
setdiff(corr_pos$target, rownames(norm))
maximum_exp <- max(norm[rownames(norm) %in% union(corr_pos$target, corr_neg$target), ])
maximum_exp
max(norm)
meta_pos <- meta
meta_pos$CRISPR_type <- paste(meta_pos$CRISPR, meta_pos$type, sep = "_")
rownames(norm)[rownames(norm) %in% corr_pos$target]
for (gene in rownames(norm)[rownames(norm) %in% corr_pos$target]) {
    meta_pos$gene <- norm[gene, ]
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
        custom_theme(diag_text = TRUE, hide_legend = TRUE)
    ggtitle(paste0("Scaled normalized expression of ", gene))
    ggsave(filename = paste0("results/images/Figure_4/coex_genes_barplot/positive_correlation/barplot_", gene, ".png"), units = "px", width = 1800, height = 1400, dpi = 250)
}

meta_neg <- meta
maximum_exp <- max(norm[rownames(norm) %in% union(corr_neg$target, corr_neg$target), ])
meta_neg$CRISPR_type <- paste(meta_neg$CRISPR, meta_neg$type, sep = "_")
rownames(norm)[rownames(norm) %in% corr_neg$target]
for (gene in rownames(norm)[rownames(norm) %in% corr_neg$target]) {
    meta_neg$gene <- norm[gene, ]
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
        custom_theme(diag_text = TRUE, hide_legend = TRUE)
    ggtitle(paste0("Scaled normalized expression of ", gene))
    ggsave(filename = paste0("results/images/Figure_4/coex_genes_barplot/negative_correlation/barplot_", gene, ".png"), units = "px", width = 1800, height = 1400, dpi = 250)
}
