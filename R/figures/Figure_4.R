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
meta <- read.table("/home/jules/Documents/phd/Data/lab_RNAseq/IPSdiff12/IPSdiff12_metadata.csv", sep = ",", header = T)
View(meta)
meta <- filter(meta, type %in% c("cyclo", "ventral"))
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
