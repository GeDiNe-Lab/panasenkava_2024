# Loading packages and functions
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(DESeq2)
library(ggrepel)

rstudioapi::getSourceEditorContext()$path %>%
    str_split("/") %>%
    unlist() %>%
    head(-3) %>%
    str_c(collapse = "/") %>%
    str_c("/") %>%
    setwd()

source("R/custom_fct.R")
# Loading data (path to change later)
rawcounts <- readcounts("/home/jules/Documents/phd/Data/lab_RNAseq/diff13/diff13_counts.csv", sep = ",", header = TRUE)
meta <- read.table("/home/jules/Documents/phd/Data/lab_RNAseq/diff13/diff13_meta.csv", sep = ",", header = TRUE)

# LON71_D12_2 is a biiiiiig outlier
meta <- filter(meta, sample != "LON71_D12_2", type %in% c("ventral", "dorsal"), line %in% c("LON71", "WTC"))
counts <- rawcounts[which(rowSums(rawcounts) >= 50), meta$sample]

# DEGs dorsal VS ventral at day02
DEGs_DV_day02 <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, day == "day02")$sample],
    colData = filter(meta, day == "day02"),
    design = ~ line + type
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "ventral", "dorsal")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_DV_day02_f <- filter(DEGs_DV_day02, padj < 0.01, abs(log2FoldChange) >= 1)
DEGs_DV_day02_f$gene <- gene_converter(rownames(DEGs_DV_day02_f), "ENSEMBL", "SYMBOL")
DEGs_DV_day02_f <- filter(DEGs_DV_day02_f, !is.na(gene))

write.csv(DEGs_DV_day02_f, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2A/DEGs_DV_day02.csv", row.names = FALSE)
# if (nrow(filter(DEGs_DV_day02_f, padj < 0.01, abs(log2FoldChange) >= 1)) > 25) {
#     DEGs_DV_day02_f$volcano <- rep(NA, nrow(DEGs_DV_day02_f))
#     FC_top25 <- sort(abs(filter(DEGs_DV_day02_f, padj < 0.01)$log2FoldChange), decreasing = TRUE)[25]
#     DEGs_DV_day02_f$volcano[abs(DEGs_DV_day02_f$log2FoldChange) >= FC_top25 & DEGs_DV_day02_f$padj < 0.01] <- DEGs_DV_day02_f$gene[abs(DEGs_DV_day02_f$log2FoldChange) >= FC_top25 & DEGs_DV_day02_f$padj < 0.01]
# } else {
#     DEGs_DV_day02_f$volcano <- DEGs_DV_day02_f$gene
# }
# png(filename = "results/images/Figure_2A/volcano_plots/DEGs_DV_day02.png", width = 1600, height = 1200, res = 250)
# ggplot(DEGs_DV_day02_f, aes(x = log2FoldChange, y = -log10(padj), label = volcano)) +
#     geom_point(size = 1) +
#     geom_text() +
#     custom_theme() +
#     geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
#     geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
#     labs(x = "log2FoldChange", y = "-log10(padj)", title = "DEGs dorsal VS ventral at day02")
# dev.off()

# DEGs dorsal VS ventral at day04
DEGs_DV_day04 <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, day == "day04")$sample],
    colData = filter(meta, day == "day04"),
    design = ~ line + type
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "ventral", "dorsal")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_DV_day04_f <- filter(DEGs_DV_day04, padj < 0.01, abs(log2FoldChange) >= 1)
DEGs_DV_day04_f$gene <- gene_converter(rownames(DEGs_DV_day04_f), "ENSEMBL", "SYMBOL")
DEGs_DV_day04_f <- filter(DEGs_DV_day04_f, !is.na(gene))

write.csv(DEGs_DV_day04_f, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2A/DEGs_DV_day04.csv", row.names = FALSE)

# DEGs dorsal VS ventral at day06
DEGs_DV_day06 <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, day == "day06")$sample],
    colData = filter(meta, day == "day06"),
    design = ~ line + type
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "ventral", "dorsal")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_DV_day06_f <- filter(DEGs_DV_day06, padj < 0.01, abs(log2FoldChange) >= 1)
DEGs_DV_day06_f$gene <- gene_converter(rownames(DEGs_DV_day06_f), "ENSEMBL", "SYMBOL")
DEGs_DV_day06_f <- filter(DEGs_DV_day06_f, !is.na(gene))

write.csv(DEGs_DV_day06_f, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2A/DEGs_DV_day06.csv", row.names = FALSE)


# DEGs dorsal VS ventral at day08
DEGs_DV_day08 <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, day == "day08")$sample],
    colData = filter(meta, day == "day08"),
    design = ~ line + type
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "ventral", "dorsal")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_DV_day08_f <- filter(DEGs_DV_day08, padj < 0.01, abs(log2FoldChange) >= 1)
DEGs_DV_day08_f$gene <- gene_converter(rownames(DEGs_DV_day08_f), "ENSEMBL", "SYMBOL")
DEGs_DV_day08_f <- filter(DEGs_DV_day08_f, !is.na(gene))

write.csv(DEGs_DV_day08_f, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2A/DEGs_DV_day08.csv", row.names = FALSE)

# DEGs dorsal VS ventral at day10
DEGs_DV_day10 <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, day == "day10")$sample],
    colData = filter(meta, day == "day10"),
    design = ~ line + type
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "ventral", "dorsal")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_DV_day10_f <- filter(DEGs_DV_day10, padj < 0.01, abs(log2FoldChange) >= 1)
DEGs_DV_day10_f$gene <- gene_converter(rownames(DEGs_DV_day10_f), "ENSEMBL", "SYMBOL")
DEGs_DV_day10_f <- filter(DEGs_DV_day10_f, !is.na(gene))

write.csv(DEGs_DV_day10_f, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2A/DEGs_DV_day10.csv", row.names = FALSE)

# DEGs dorsal VS ventral at day12
DEGs_DV_day12 <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, day == "day12")$sample],
    colData = filter(meta, day == "day12"),
    design = ~ line + type
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "ventral", "dorsal")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_DV_day12_f <- filter(DEGs_DV_day12, padj < 0.01, abs(log2FoldChange) >= 1)
DEGs_DV_day12_f$gene <- gene_converter(rownames(DEGs_DV_day12_f), "ENSEMBL", "SYMBOL")
DEGs_DV_day12_f <- filter(DEGs_DV_day12_f, !is.na(gene))

write.csv(DEGs_DV_day12_f, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2A/DEGs_DV_day12.csv", row.names = FALSE)



# DEGs dorsal VS ventral at day12
DEGs_day02_04_dorsal <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "dorsal")$sample],
    colData = filter(meta, type == "dorsal"),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day02", "day04")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day02_04_dorsal_f <- filter(DEGs_day02_04_dorsal, padj < 0.01, abs(log2FoldChange) >= 1)
DEGs_day02_04_dorsal_f$gene <- gene_converter(rownames(DEGs_day02_04_dorsal_f), "ENSEMBL", "SYMBOL")
DEGs_day02_04_dorsal_f <- filter(DEGs_day02_04_dorsal_f, !is.na(gene))

write.csv(DEGs_day02_04_dorsal_f, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2A/DEGs_day02_04_dorsal.csv", row.names = FALSE)


DEGs_day02_04_ventral <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "ventral")$sample],
    colData = filter(meta, type == "ventral"),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day02", "day04")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day02_04_ventral_f <- filter(DEGs_day02_04_ventral, padj < 0.01, abs(log2FoldChange) >= 1)
DEGs_day02_04_ventral_f$gene <- gene_converter(rownames(DEGs_day02_04_ventral_f), "ENSEMBL", "SYMBOL")
DEGs_day02_04_ventral_f <- filter(DEGs_day02_04_ventral_f, !is.na(gene))

write.csv(DEGs_day02_04_ventral_f, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2A/DEGs_day02_04_ventral.csv", row.names = FALSE)


# day 4 to 8
DEGs_day04_06_dorsal <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "dorsal")$sample],
    colData = filter(meta, type == "dorsal"),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day04", "day06")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day04_06_dorsal_f <- filter(DEGs_day04_06_dorsal, padj < 0.01, abs(log2FoldChange) >= 1)
DEGs_day04_06_dorsal_f$gene <- gene_converter(rownames(DEGs_day04_06_dorsal_f), "ENSEMBL", "SYMBOL")
DEGs_day04_06_dorsal_f <- filter(DEGs_day04_06_dorsal_f, !is.na(gene))

write.csv(DEGs_day04_06_dorsal_f, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2A/DEGs_day04_06_dorsal.csv", row.names = FALSE)

DEGs_day04_06_ventral <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "ventral")$sample],
    colData = filter(meta, type == "ventral"),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day04", "day06")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day04_06_ventral_f <- filter(DEGs_day04_06_ventral, padj < 0.01, abs(log2FoldChange) >= 1)
DEGs_day04_06_ventral_f$gene <- gene_converter(rownames(DEGs_day04_06_ventral_f), "ENSEMBL", "SYMBOL")
DEGs_day04_06_ventral_f <- filter(DEGs_day04_06_ventral_f, !is.na(gene))

write.csv(DEGs_day04_06_ventral_f, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2A/DEGs_day04_06_ventral.csv", row.names = FALSE)

# day 6 to 8
DEGs_day06_08_dorsal <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "dorsal")$sample],
    colData = filter(meta, type == "dorsal"),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day06", "day08")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day06_08_dorsal_f <- filter(DEGs_day06_08_dorsal, padj < 0.01, abs(log2FoldChange) >= 1)
DEGs_day06_08_dorsal_f$gene <- gene_converter(rownames(DEGs_day06_08_dorsal_f), "ENSEMBL", "SYMBOL")
DEGs_day06_08_dorsal_f <- filter(DEGs_day06_08_dorsal_f, !is.na(gene))

write.csv(DEGs_day06_08_dorsal_f, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2A/DEGs_day06_08_dorsal.csv", row.names = FALSE)

DEGs_day06_08_ventral <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "ventral")$sample],
    colData = filter(meta, type == "ventral"),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day06", "day08")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day06_08_ventral_f <- filter(DEGs_day06_08_ventral, padj < 0.01, abs(log2FoldChange) >= 1)
DEGs_day06_08_ventral_f$gene <- gene_converter(rownames(DEGs_day06_08_ventral_f), "ENSEMBL", "SYMBOL")
DEGs_day06_08_ventral_f <- filter(DEGs_day06_08_ventral_f, !is.na(gene))

write.csv(DEGs_day06_08_ventral_f, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2A/DEGs_day06_08_ventral.csv", row.names = FALSE)

# day 6 to 8
DEGs_day08_10_dorsal <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "dorsal")$sample],
    colData = filter(meta, type == "dorsal"),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day08", "day10")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day08_10_dorsal_f <- filter(DEGs_day08_10_dorsal, padj < 0.01, abs(log2FoldChange) >= 1)
DEGs_day08_10_dorsal_f$gene <- gene_converter(rownames(DEGs_day08_10_dorsal_f), "ENSEMBL", "SYMBOL")
DEGs_day08_10_dorsal_f <- filter(DEGs_day08_10_dorsal_f, !is.na(gene))

write.csv(DEGs_day08_10_dorsal_f, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2A/DEGs_day08_10_dorsal.csv", row.names = FALSE)

DEGs_day08_10_ventral <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "ventral")$sample],
    colData = filter(meta, type == "ventral"),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day08", "day10")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day08_10_ventral_f <- filter(DEGs_day08_10_ventral, padj < 0.01, abs(log2FoldChange) >= 1)
DEGs_day08_10_ventral_f$gene <- gene_converter(rownames(DEGs_day08_10_ventral_f), "ENSEMBL", "SYMBOL")
DEGs_day08_10_ventral_f <- filter(DEGs_day08_10_ventral_f, !is.na(gene))

write.csv(DEGs_day08_10_ventral_f, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2A/DEGs_day08_10_ventral.csv", row.names = FALSE)

# day 6 to 8
DEGs_day10_12_dorsal <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "dorsal")$sample],
    colData = filter(meta, type == "dorsal"),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day10", "day12")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day10_12_dorsal_f <- filter(DEGs_day10_12_dorsal, padj < 0.01, abs(log2FoldChange) >= 1)
DEGs_day10_12_dorsal_f$gene <- gene_converter(rownames(DEGs_day10_12_dorsal_f), "ENSEMBL", "SYMBOL")
DEGs_day10_12_dorsal_f <- filter(DEGs_day10_12_dorsal_f, !is.na(gene))

write.csv(DEGs_day10_12_dorsal_f, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2A/DEGs_day10_12_dorsal.csv", row.names = FALSE)

DEGs_day10_12_ventral <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "ventral")$sample],
    colData = filter(meta, type == "ventral"),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day10", "day12")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day10_12_ventral_f <- filter(DEGs_day10_12_ventral, padj < 0.01, abs(log2FoldChange) >= 1)
DEGs_day10_12_ventral_f$gene <- gene_converter(rownames(DEGs_day10_12_ventral_f), "ENSEMBL", "SYMBOL")
DEGs_day10_12_ventral_f <- filter(DEGs_day10_12_ventral_f, !is.na(gene))

write.csv(DEGs_day10_12_ventral_f, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2A/DEGs_day10_12_ventral.csv", row.names = FALSE)
