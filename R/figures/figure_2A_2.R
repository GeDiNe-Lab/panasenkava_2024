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
DEGs_dorso_ventral <- DESeqDataSetFromMatrix(
    countData = counts[, meta$sample],
    colData = meta,
    design = ~ line + type
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "ventral", "dorsal")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_dorso_ventral_f <- filter(DEGs_dorso_ventral, padj < 0.01, abs(log2FoldChange) >= 1)
DEGs_dorso_ventral_f$gene <- gene_converter(rownames(DEGs_dorso_ventral_f), "ENSEMBL", "SYMBOL")
DEGs_dorso_ventral_f <- filter(DEGs_dorso_ventral_f, !is.na(gene))

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

DE_days_ventral <- list(
    day02_04 = DEGs_day02_04_ventral_f,
    day04_06 = DEGs_day04_06_ventral_f,
    day06_08 = DEGs_day06_08_ventral_f,
    day08_10 = DEGs_day08_10_ventral_f,
    day10_12 = DEGs_day10_12_ventral_f
)
names(DE_days_ventral)

for (dayrange in names(DE_days_ventral)) {
    print(dayrange)
    DE <- DE_days_ventral[[dayrange]]
    ggplot(DE, aes(x = log2FoldChange, y = -log10(padj), label = gene)) +
        ggrepel::geom_text_repel(box.padding = 0.001, size = 2.5, max.overlaps = 20) +
        custom_theme() +
        geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
        labs(x = "log2FoldChange", y = "-log10(padj)", title = paste0("DEGs for ventral samples between ", dayrange))
    ggsave(filename = paste0("results/images/Figure_2A/volcano_plots/DEGs_ventral_", dayrange, ".png"), units = "px", width = 1800, height = 1400, dpi = 250)
}


DE_days_dorsal <- list(
    day02_04 = DEGs_day02_04_dorsal_f,
    day04_06 = DEGs_day04_06_dorsal_f,
    day06_08 = DEGs_day06_08_dorsal_f,
    day08_10 = DEGs_day08_10_dorsal_f,
    day10_12 = DEGs_day10_12_dorsal_f
)
names(DE_days_dorsal)

for (dayrange in names(DE_days_dorsal)) {
    print(dayrange)
    DE <- DE_days_dorsal[[dayrange]]
    ggplot(DE, aes(x = log2FoldChange, y = -log10(padj), label = gene)) +
        ggrepel::geom_text_repel(box.padding = 0.001, size = 2.5, max.overlaps = 20) +
        custom_theme() +
        geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
        labs(x = "log2FoldChange", y = "-log10(padj)", title = paste0("DEGs for dorsal samples between ", dayrange))
    ggsave(filename = paste0("results/images/Figure_2A/volcano_plots/DEGs_dorsal_", dayrange, ".png"), units = "px", width = 1800, height = 1400, dpi = 250)
}
