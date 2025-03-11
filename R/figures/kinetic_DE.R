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

# Loading custom functions
source("R/custom_fct.R")
# Loading data (path to change later)
# Loading data (path to change later)
rawcounts <- readcounts("data/rawcounts.csv", sep = ",", header = TRUE)
rawmeta <- read.table("data/meta.csv", sep = ",", header = TRUE)

# LON71_D12_2 does not have any reads in the count file
# though, the fastQC report shows that the sample is good
meta <- filter(rawmeta, sample != "LON71_D12_2", diff == "diff13", line %in% c("LON71", "WTC"), ((manip == "veranika" & day != "day12") | (manip == "lauryane" & day == "day12")))

# filtering out lowly expressed genes
counts <- rawcounts[, meta$sample][which(rowSums(rawcounts[, meta$sample]) >= 25), ]

# DEGs ventral VS dorsal with all samples
DEGs_vAN_vs_dAN <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ line + type
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "ventral", "dorsal")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_vAN_vs_dAN$gene <- gene_converter(rownames(DEGs_vAN_vs_dAN), "ENSEMBL", "SYMBOL")
DEGs_vAN_vs_dAN_f <- filter(DEGs_vAN_vs_dAN, !is.na(gene))
DEGs_vAN_vs_dAN_f <- filter(DEGs_vAN_vs_dAN_f, padj < 0.01, abs(log2FoldChange) >= 2)

DEGs_vAN_vs_dAN$classic_threshold <- ifelse(DEGs_vAN_vs_dAN$padj < 0.01 & abs(DEGs_vAN_vs_dAN$log2FoldChange) >= 1, "yes", "no")
DEGs_vAN_vs_dAN$high_threshold <- ifelse(DEGs_vAN_vs_dAN$padj < 0.01 & abs(DEGs_vAN_vs_dAN$log2FoldChange) >= 2, "yes", "no")
write.csv(DEGs_vAN_vs_dAN, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/DEGs_vAN_vs_dAN.csv")

# DEGs ventral VS dorsal at day02
DEGs_vAN_vs_dAN_day02 <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, day == "day02")$sample],
    colData = filter(meta, day == "day02"),
    design = ~ line + type
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "ventral", "dorsal")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_vAN_vs_dAN_day02$gene <- gene_converter(rownames(DEGs_vAN_vs_dAN_day02), "ENSEMBL", "SYMBOL")
DEGs_vAN_vs_dAN_day02_f <- filter(DEGs_vAN_vs_dAN_day02, !is.na(gene))
DEGs_vAN_vs_dAN_day02_f <- filter(DEGs_vAN_vs_dAN_day02_f, padj < 0.01, abs(log2FoldChange) >= 2)

DEGs_vAN_vs_dAN_day02$classic_threshold <- ifelse(DEGs_vAN_vs_dAN_day02$padj < 0.01 & abs(DEGs_vAN_vs_dAN_day02$log2FoldChange) >= 1, "yes", "no")
DEGs_vAN_vs_dAN_day02$high_threshold <- ifelse(DEGs_vAN_vs_dAN_day02$padj < 0.01 & abs(DEGs_vAN_vs_dAN_day02$log2FoldChange) >= 2, "yes", "no")
write.csv(DEGs_vAN_vs_dAN_day02, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/DEGs_vAN_vs_dAN_day02.csv")


# DEGs ventral VS dorsal at day04
DEGs_vAN_vs_dAN_day04 <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, day == "day04")$sample],
    colData = filter(meta, day == "day04"),
    design = ~ line + type
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "ventral", "dorsal")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_vAN_vs_dAN_day04$gene <- gene_converter(rownames(DEGs_vAN_vs_dAN_day04), "ENSEMBL", "SYMBOL")
DEGs_vAN_vs_dAN_day04_f <- filter(DEGs_vAN_vs_dAN_day04, !is.na(gene))
DEGs_vAN_vs_dAN_day04_f <- filter(DEGs_vAN_vs_dAN_day04_f, padj < 0.01, abs(log2FoldChange) >= 2)

DEGs_vAN_vs_dAN_day04$classic_threshold <- ifelse(DEGs_vAN_vs_dAN_day04$padj < 0.01 & abs(DEGs_vAN_vs_dAN_day04$log2FoldChange) >= 1, "yes", "no")
DEGs_vAN_vs_dAN_day04$high_threshold <- ifelse(DEGs_vAN_vs_dAN_day04$padj < 0.01 & abs(DEGs_vAN_vs_dAN_day04$log2FoldChange) >= 2, "yes", "no")
write.csv(DEGs_vAN_vs_dAN_day04, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/DEGs_vAN_vs_dAN_day04.csv")

# DEGs ventral VS dorsal at day06
DEGs_vAN_vs_dAN_day06 <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, day == "day06")$sample],
    colData = filter(meta, day == "day06"),
    design = ~ line + type
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "ventral", "dorsal")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_vAN_vs_dAN_day06$gene <- gene_converter(rownames(DEGs_vAN_vs_dAN_day06), "ENSEMBL", "SYMBOL")
DEGs_vAN_vs_dAN_day06_f <- filter(DEGs_vAN_vs_dAN_day06, !is.na(gene))
DEGs_vAN_vs_dAN_day06_f <- filter(DEGs_vAN_vs_dAN_day06_f, padj < 0.01, abs(log2FoldChange) >= 2)

DEGs_vAN_vs_dAN_day06$classic_threshold <- ifelse(DEGs_vAN_vs_dAN_day06$padj < 0.01 & abs(DEGs_vAN_vs_dAN_day06$log2FoldChange) >= 1, "yes", "no")
DEGs_vAN_vs_dAN_day06$high_threshold <- ifelse(DEGs_vAN_vs_dAN_day06$padj < 0.01 & abs(DEGs_vAN_vs_dAN_day06$log2FoldChange) >= 2, "yes", "no")
write.csv(DEGs_vAN_vs_dAN_day06, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/DEGs_vAN_vs_dAN_day06.csv")


# DEGs ventral VS dorsal at day08
DEGs_vAN_vs_dAN_day08 <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, day == "day08")$sample],
    colData = filter(meta, day == "day08"),
    design = ~ line + type
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "ventral", "dorsal")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_vAN_vs_dAN_day08$gene <- gene_converter(rownames(DEGs_vAN_vs_dAN_day08), "ENSEMBL", "SYMBOL")
DEGs_vAN_vs_dAN_day08_f <- filter(DEGs_vAN_vs_dAN_day08, !is.na(gene))
DEGs_vAN_vs_dAN_day08_f <- filter(DEGs_vAN_vs_dAN_day08_f, padj < 0.01, abs(log2FoldChange) >= 2)

DEGs_vAN_vs_dAN_day08$classic_threshold <- ifelse(DEGs_vAN_vs_dAN_day08$padj < 0.01 & abs(DEGs_vAN_vs_dAN_day08$log2FoldChange) >= 1, "yes", "no")
DEGs_vAN_vs_dAN_day08$high_threshold <- ifelse(DEGs_vAN_vs_dAN_day08$padj < 0.01 & abs(DEGs_vAN_vs_dAN_day08$log2FoldChange) >= 2, "yes", "no")
write.csv(DEGs_vAN_vs_dAN_day08, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/DEGs_vAN_vs_dAN_day08.csv")

# DEGs ventral VS dorsal at day10
DEGs_vAN_vs_dAN_day10 <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, day == "day10")$sample],
    colData = filter(meta, day == "day10"),
    design = ~ line + type
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "ventral", "dorsal")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_vAN_vs_dAN_day10$gene <- gene_converter(rownames(DEGs_vAN_vs_dAN_day10), "ENSEMBL", "SYMBOL")
DEGs_vAN_vs_dAN_day10_f <- filter(DEGs_vAN_vs_dAN_day10, !is.na(gene))
DEGs_vAN_vs_dAN_day10_f <- filter(DEGs_vAN_vs_dAN_day10_f, padj < 0.01, abs(log2FoldChange) >= 2)

DEGs_vAN_vs_dAN_day10$classic_threshold <- ifelse(DEGs_vAN_vs_dAN_day10$padj < 0.01 & abs(DEGs_vAN_vs_dAN_day10$log2FoldChange) >= 1, "yes", "no")
DEGs_vAN_vs_dAN_day10$high_threshold <- ifelse(DEGs_vAN_vs_dAN_day10$padj < 0.01 & abs(DEGs_vAN_vs_dAN_day10$log2FoldChange) >= 2, "yes", "no")
write.csv(DEGs_vAN_vs_dAN_day10, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/DEGs_vAN_vs_dAN_day10.csv")

# DEGs ventral VS dorsal at day12
DEGs_vAN_vs_dAN_day12 <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, day == "day12")$sample],
    colData = filter(meta, day == "day12"),
    design = ~ line + type
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "ventral", "dorsal")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_vAN_vs_dAN_day12$gene <- gene_converter(rownames(DEGs_vAN_vs_dAN_day12), "ENSEMBL", "SYMBOL")
DEGs_vAN_vs_dAN_day12_f <- filter(DEGs_vAN_vs_dAN_day12, !is.na(gene))
DEGs_vAN_vs_dAN_day12_f <- filter(DEGs_vAN_vs_dAN_day12_f, padj < 0.01, abs(log2FoldChange) >= 2)

DEGs_vAN_vs_dAN_day12$classic_threshold <- ifelse(DEGs_vAN_vs_dAN_day12$padj < 0.01 & abs(DEGs_vAN_vs_dAN_day12$log2FoldChange) >= 1, "yes", "no")
DEGs_vAN_vs_dAN_day12$high_threshold <- ifelse(DEGs_vAN_vs_dAN_day12$padj < 0.01 & abs(DEGs_vAN_vs_dAN_day12$log2FoldChange) >= 2, "yes", "no")
write.csv(DEGs_vAN_vs_dAN_day12, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/DEGs_vAN_vs_dAN_day12.csv")


genes_ventral <- list(
    day_02 = filter(DEGs_vAN_vs_dAN_day02, log2FoldChange >= 2 & padj < 0.01 & !is.na(gene))$gene,
    day_04 = filter(DEGs_vAN_vs_dAN_day04, log2FoldChange >= 2 & padj < 0.01 & !is.na(gene))$gene,
    day_06 = filter(DEGs_vAN_vs_dAN_day06, log2FoldChange >= 2 & padj < 0.01 & !is.na(gene))$gene,
    day_08 = filter(DEGs_vAN_vs_dAN_day08, log2FoldChange >= 2 & padj < 0.01 & !is.na(gene))$gene,
    day_10 = filter(DEGs_vAN_vs_dAN_day10, log2FoldChange >= 2 & padj < 0.01 & !is.na(gene))$gene,
    day_12 = filter(DEGs_vAN_vs_dAN_day12, log2FoldChange >= 2 & padj < 0.01 & !is.na(gene))$gene
) %>% Reduce(union, .)
DE_df_ventral <- data.frame(
    gene = genes_ventral,
    day02 = ifelse(genes_ventral %in% filter(DEGs_vAN_vs_dAN_day02, log2FoldChange >= 2 & padj < 0.01 & !is.na(gene))$gene, "YES", "no"),
    day04 = ifelse(genes_ventral %in% filter(DEGs_vAN_vs_dAN_day04, log2FoldChange >= 2 & padj < 0.01 & !is.na(gene))$gene, "YES", "no"),
    day06 = ifelse(genes_ventral %in% filter(DEGs_vAN_vs_dAN_day06, log2FoldChange >= 2 & padj < 0.01 & !is.na(gene))$gene, "YES", "no"),
    day08 = ifelse(genes_ventral %in% filter(DEGs_vAN_vs_dAN_day08, log2FoldChange >= 2 & padj < 0.01 & !is.na(gene))$gene, "YES", "no"),
    day10 = ifelse(genes_ventral %in% filter(DEGs_vAN_vs_dAN_day10, log2FoldChange >= 2 & padj < 0.01 & !is.na(gene))$gene, "YES", "no"),
    day12 = ifelse(genes_ventral %in% filter(DEGs_vAN_vs_dAN_day12, log2FoldChange >= 2 & padj < 0.01 & !is.na(gene))$gene, "YES", "no")
)
write.csv(DE_df_ventral, file = "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/highly_DEG_by_day_ventral.csv", quote = FALSE)


genes_dorsal <- list(
    day_02 = filter(DEGs_vAN_vs_dAN_day02, log2FoldChange <= -2 & padj < 0.01 & !is.na(gene))$gene,
    day_04 = filter(DEGs_vAN_vs_dAN_day04, log2FoldChange <= -2 & padj < 0.01 & !is.na(gene))$gene,
    day_06 = filter(DEGs_vAN_vs_dAN_day06, log2FoldChange <= -2 & padj < 0.01 & !is.na(gene))$gene,
    day_08 = filter(DEGs_vAN_vs_dAN_day08, log2FoldChange <= -2 & padj < 0.01 & !is.na(gene))$gene,
    day_10 = filter(DEGs_vAN_vs_dAN_day10, log2FoldChange <= -2 & padj < 0.01 & !is.na(gene))$gene,
    day_12 = filter(DEGs_vAN_vs_dAN_day12, log2FoldChange <= -2 & padj < 0.01 & !is.na(gene))$gene
) %>% Reduce(union, .)
DE_df_dorsal <- data.frame(
    gene = genes_dorsal,
    day02 = ifelse(genes_dorsal %in% filter(DEGs_vAN_vs_dAN_day02, log2FoldChange <= -2 & padj < 0.01 & !is.na(gene))$gene, "YES", "no"),
    day04 = ifelse(genes_dorsal %in% filter(DEGs_vAN_vs_dAN_day04, log2FoldChange <= -2 & padj < 0.01 & !is.na(gene))$gene, "YES", "no"),
    day06 = ifelse(genes_dorsal %in% filter(DEGs_vAN_vs_dAN_day06, log2FoldChange <= -2 & padj < 0.01 & !is.na(gene))$gene, "YES", "no"),
    day08 = ifelse(genes_dorsal %in% filter(DEGs_vAN_vs_dAN_day08, log2FoldChange <= -2 & padj < 0.01 & !is.na(gene))$gene, "YES", "no"),
    day10 = ifelse(genes_dorsal %in% filter(DEGs_vAN_vs_dAN_day10, log2FoldChange <= -2 & padj < 0.01 & !is.na(gene))$gene, "YES", "no"),
    day12 = ifelse(genes_dorsal %in% filter(DEGs_vAN_vs_dAN_day12, log2FoldChange <= -2 & padj < 0.01 & !is.na(gene))$gene, "YES", "no")
)
write.csv(DE_df_dorsal, file = "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/highly_DEG_by_day_dorsal.csv", quote = FALSE)


# DEGs day04 vs day02 for dorsal samples
DEGs_day_04_vs_02_dorsal <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "dorsal")$sample],
    colData = filter(meta, type == "dorsal"),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day04", "day02")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day_04_vs_02_dorsal$gene <- gene_converter(rownames(DEGs_day_04_vs_02_dorsal), "ENSEMBL", "SYMBOL")
DEGs_day_04_vs_02_dorsal_f <- filter(DEGs_day_04_vs_02_dorsal, !is.na(gene))
DEGs_day_04_vs_02_dorsal_f <- filter(DEGs_day_04_vs_02_dorsal_f, padj < 0.01, abs(log2FoldChange) >= 2)

DEGs_day_04_vs_02_dorsal$classic_threshold <- ifelse(DEGs_day_04_vs_02_dorsal$padj < 0.01 & abs(DEGs_day_04_vs_02_dorsal$log2FoldChange) >= 1, "yes", "no")
DEGs_day_04_vs_02_dorsal$high_threshold <- ifelse(DEGs_day_04_vs_02_dorsal$padj < 0.01 & abs(DEGs_day_04_vs_02_dorsal$log2FoldChange) >= 2, "yes", "no")
write.csv(DEGs_day_04_vs_02_dorsal, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/DEGs_day_04_vs_02_dorsal.csv")

# DEGs day04 vs day02 for ventral samples
DEGs_day_04_vs_02_ventral <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "ventral")$sample],
    colData = filter(meta, type == "ventral"),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day04", "day02")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day_04_vs_02_ventral$gene <- gene_converter(rownames(DEGs_day_04_vs_02_ventral), "ENSEMBL", "SYMBOL")
DEGs_day_04_vs_02_ventral_f <- filter(DEGs_day_04_vs_02_ventral, !is.na(gene))
DEGs_day_04_vs_02_ventral_f <- filter(DEGs_day_04_vs_02_ventral_f, padj < 0.01, abs(log2FoldChange) >= 2)

DEGs_day_04_vs_02_ventral$classic_threshold <- ifelse(DEGs_day_04_vs_02_ventral$padj < 0.01 & abs(DEGs_day_04_vs_02_ventral$log2FoldChange) >= 1, "yes", "no")
DEGs_day_04_vs_02_ventral$high_threshold <- ifelse(DEGs_day_04_vs_02_ventral$padj < 0.01 & abs(DEGs_day_04_vs_02_ventral$log2FoldChange) >= 2, "yes", "no")
write.csv(DEGs_day_04_vs_02_ventral, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/DEGs_day_04_vs_02_ventral.csv")

# DEGs day06 vs day04 for dorsal samples
DEGs_day_06_vs_04_dorsal <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "dorsal")$sample],
    colData = filter(meta, type == "dorsal"),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day06", "day04")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day_06_vs_04_dorsal$gene <- gene_converter(rownames(DEGs_day_06_vs_04_dorsal), "ENSEMBL", "SYMBOL")
DEGs_day_06_vs_04_dorsal_f <- filter(DEGs_day_06_vs_04_dorsal, !is.na(gene))
DEGs_day_06_vs_04_dorsal_f <- filter(DEGs_day_06_vs_04_dorsal_f, padj < 0.01, abs(log2FoldChange) >= 2)

DEGs_day_06_vs_04_dorsal$classic_threshold <- ifelse(DEGs_day_06_vs_04_dorsal$padj < 0.01 & abs(DEGs_day_06_vs_04_dorsal$log2FoldChange) >= 1, "yes", "no")
DEGs_day_06_vs_04_dorsal$high_threshold <- ifelse(DEGs_day_06_vs_04_dorsal$padj < 0.01 & abs(DEGs_day_06_vs_04_dorsal$log2FoldChange) >= 2, "yes", "no")
write.csv(DEGs_day_06_vs_04_dorsal, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/DEGs_day_06_vs_04_dorsal.csv")

# DEGs day06 vs day04 for ventral samples
DEGs_day_06_vs_04_ventral <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "ventral")$sample],
    colData = filter(meta, type == "ventral"),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day06", "day04")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day_06_vs_04_ventral$gene <- gene_converter(rownames(DEGs_day_06_vs_04_ventral), "ENSEMBL", "SYMBOL")
DEGs_day_06_vs_04_ventral_f <- filter(DEGs_day_06_vs_04_ventral, !is.na(gene))
DEGs_day_06_vs_04_ventral_f <- filter(DEGs_day_06_vs_04_ventral_f, padj < 0.01, abs(log2FoldChange) >= 2)

DEGs_day_06_vs_04_ventral$classic_threshold <- ifelse(DEGs_day_06_vs_04_ventral$padj < 0.01 & abs(DEGs_day_06_vs_04_ventral$log2FoldChange) >= 1, "yes", "no")
DEGs_day_06_vs_04_ventral$high_threshold <- ifelse(DEGs_day_06_vs_04_ventral$padj < 0.01 & abs(DEGs_day_06_vs_04_ventral$log2FoldChange) >= 2, "yes", "no")
write.csv(DEGs_day_06_vs_04_ventral, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/DEGs_day_06_vs_04_ventral.csv")

# DEGs day08 vs day06 for dorsal samples
DEGs_day_08_vs_06_dorsal <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "dorsal")$sample],
    colData = filter(meta, type == "dorsal"),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day08", "day06")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day_08_vs_06_dorsal$gene <- gene_converter(rownames(DEGs_day_08_vs_06_dorsal), "ENSEMBL", "SYMBOL")
DEGs_day_08_vs_06_dorsal_f <- filter(DEGs_day_08_vs_06_dorsal, !is.na(gene))
DEGs_day_08_vs_06_dorsal_f <- filter(DEGs_day_08_vs_06_dorsal_f, padj < 0.01, abs(log2FoldChange) >= 2)

DEGs_day_08_vs_06_dorsal$classic_threshold <- ifelse(DEGs_day_08_vs_06_dorsal$padj < 0.01 & abs(DEGs_day_08_vs_06_dorsal$log2FoldChange) >= 1, "yes", "no")
DEGs_day_08_vs_06_dorsal$high_threshold <- ifelse(DEGs_day_08_vs_06_dorsal$padj < 0.01 & abs(DEGs_day_08_vs_06_dorsal$log2FoldChange) >= 2, "yes", "no")
write.csv(DEGs_day_08_vs_06_dorsal, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/DEGs_day_08_vs_06_dorsal.csv")

# DEGs day08 vs day06 for ventral samples
DEGs_day_08_vs_06_ventral <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "ventral")$sample],
    colData = filter(meta, type == "ventral"),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day08", "day06")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day_08_vs_06_ventral$gene <- gene_converter(rownames(DEGs_day_08_vs_06_ventral), "ENSEMBL", "SYMBOL")
DEGs_day_08_vs_06_ventral_f <- filter(DEGs_day_08_vs_06_ventral, !is.na(gene))
DEGs_day_08_vs_06_ventral_f <- filter(DEGs_day_08_vs_06_ventral_f, padj < 0.01, abs(log2FoldChange) >= 2)

DEGs_day_08_vs_06_ventral$classic_threshold <- ifelse(DEGs_day_08_vs_06_ventral$padj < 0.01 & abs(DEGs_day_08_vs_06_ventral$log2FoldChange) >= 1, "yes", "no")
DEGs_day_08_vs_06_ventral$high_threshold <- ifelse(DEGs_day_08_vs_06_ventral$padj < 0.01 & abs(DEGs_day_08_vs_06_ventral$log2FoldChange) >= 2, "yes", "no")
write.csv(DEGs_day_08_vs_06_ventral, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/DEGs_day_08_vs_06_ventral.csv")

# DEGs day10 vs day08 for dorsal samples
DEGs_day_10_vs_08_dorsal <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "dorsal")$sample],
    colData = filter(meta, type == "dorsal"),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day10", "day08")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day_10_vs_08_dorsal$gene <- gene_converter(rownames(DEGs_day_10_vs_08_dorsal), "ENSEMBL", "SYMBOL")
DEGs_day_10_vs_08_dorsal_f <- filter(DEGs_day_10_vs_08_dorsal, !is.na(gene))
DEGs_day_10_vs_08_dorsal_f <- filter(DEGs_day_10_vs_08_dorsal_f, padj < 0.01, abs(log2FoldChange) >= 2)

DEGs_day_10_vs_08_dorsal$classic_threshold <- ifelse(DEGs_day_10_vs_08_dorsal$padj < 0.01 & abs(DEGs_day_10_vs_08_dorsal$log2FoldChange) >= 1, "yes", "no")
DEGs_day_10_vs_08_dorsal$high_threshold <- ifelse(DEGs_day_10_vs_08_dorsal$padj < 0.01 & abs(DEGs_day_10_vs_08_dorsal$log2FoldChange) >= 2, "yes", "no")
write.csv(DEGs_day_10_vs_08_dorsal, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/DEGs_day_10_vs_08_dorsal.csv")

# DEGs day10 vs day08 for ventral samples
DEGs_day_10_vs_08_ventral <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "ventral")$sample],
    colData = filter(meta, type == "ventral"),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day10", "day08")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day_10_vs_08_ventral$gene <- gene_converter(rownames(DEGs_day_10_vs_08_ventral), "ENSEMBL", "SYMBOL")
DEGs_day_10_vs_08_ventral_f <- filter(DEGs_day_10_vs_08_ventral, !is.na(gene))
DEGs_day_10_vs_08_ventral_f <- filter(DEGs_day_10_vs_08_ventral_f, padj < 0.01, abs(log2FoldChange) >= 2)

DEGs_day_10_vs_08_ventral$classic_threshold <- ifelse(DEGs_day_10_vs_08_ventral$padj < 0.01 & abs(DEGs_day_10_vs_08_ventral$log2FoldChange) >= 1, "yes", "no")
DEGs_day_10_vs_08_ventral$high_threshold <- ifelse(DEGs_day_10_vs_08_ventral$padj < 0.01 & abs(DEGs_day_10_vs_08_ventral$log2FoldChange) >= 2, "yes", "no")
write.csv(DEGs_day_10_vs_08_ventral, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/DEGs_day_10_vs_08_ventral.csv")

# DEGs day12 vs day10 for dorsal samples
DEGs_day_12_vs_10_dorsal <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "dorsal")$sample],
    colData = filter(meta, type == "dorsal"),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day12", "day10")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day_12_vs_10_dorsal$gene <- gene_converter(rownames(DEGs_day_12_vs_10_dorsal), "ENSEMBL", "SYMBOL")
DEGs_day_12_vs_10_dorsal_f <- filter(DEGs_day_12_vs_10_dorsal, !is.na(gene))
DEGs_day_12_vs_10_dorsal_f <- filter(DEGs_day_12_vs_10_dorsal_f, padj < 0.01, abs(log2FoldChange) >= 2)

DEGs_day_12_vs_10_dorsal$classic_threshold <- ifelse(DEGs_day_12_vs_10_dorsal$padj < 0.01 & abs(DEGs_day_12_vs_10_dorsal$log2FoldChange) >= 1, "yes", "no")
DEGs_day_12_vs_10_dorsal$high_threshold <- ifelse(DEGs_day_12_vs_10_dorsal$padj < 0.01 & abs(DEGs_day_12_vs_10_dorsal$log2FoldChange) >= 2, "yes", "no")
write.csv(DEGs_day_12_vs_10_dorsal, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/DEGs_day_12_vs_10_dorsal.csv")

# DEGs day12 vs day10 for ventral samples
DEGs_day_12_vs_10_ventral <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "ventral")$sample],
    colData = filter(meta, type == "ventral"),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day12", "day10")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day_12_vs_10_ventral$gene <- gene_converter(rownames(DEGs_day_12_vs_10_ventral), "ENSEMBL", "SYMBOL")
DEGs_day_12_vs_10_ventral_f <- filter(DEGs_day_12_vs_10_ventral, !is.na(gene))
DEGs_day_12_vs_10_ventral_f <- filter(DEGs_day_12_vs_10_ventral_f, padj < 0.01, abs(log2FoldChange) >= 2)

DEGs_day_12_vs_10_ventral$classic_threshold <- ifelse(DEGs_day_12_vs_10_ventral$padj < 0.01 & abs(DEGs_day_12_vs_10_ventral$log2FoldChange) >= 1, "yes", "no")
DEGs_day_12_vs_10_ventral$high_threshold <- ifelse(DEGs_day_12_vs_10_ventral$padj < 0.01 & abs(DEGs_day_12_vs_10_ventral$log2FoldChange) >= 2, "yes", "no")
write.csv(DEGs_day_12_vs_10_ventral, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/DEGs_day_12_vs_10_ventral.csv")

# DEGs day12 vs day6 for ventral samples
DEGs_day_12_vs_06_ventral <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "ventral")$sample],
    colData = filter(meta, type == "ventral"),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day12", "day06")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day_12_vs_06_ventral$gene <- gene_converter(rownames(DEGs_day_12_vs_06_ventral), "ENSEMBL", "SYMBOL")
DEGs_day_12_vs_06_ventral_f <- filter(DEGs_day_12_vs_06_ventral, !is.na(gene))
DEGs_day_12_vs_06_ventral_f <- filter(DEGs_day_12_vs_06_ventral_f, padj < 0.01, abs(log2FoldChange) >= 2)

DEGs_day_12_vs_06_ventral$classic_threshold <- ifelse(DEGs_day_12_vs_06_ventral$padj < 0.01 & abs(DEGs_day_12_vs_06_ventral$log2FoldChange) >= 1, "yes", "no")
DEGs_day_12_vs_06_ventral$high_threshold <- ifelse(DEGs_day_12_vs_06_ventral$padj < 0.01 & abs(DEGs_day_12_vs_06_ventral$log2FoldChange) >= 2, "yes", "no")
write.csv(DEGs_day_12_vs_06_ventral, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/DEGs_day_12_vs_06_ventral.csv")

# DEGs day12 vs day6 for dorsal samples
DEGs_day_12_vs_06_dorsal <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "dorsal")$sample],
    colData = filter(meta, type == "dorsal"),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day12", "day06")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day_12_vs_06_dorsal$gene <- gene_converter(rownames(DEGs_day_12_vs_06_dorsal), "ENSEMBL", "SYMBOL")
DEGs_day_12_vs_06_dorsal_f <- filter(DEGs_day_12_vs_06_dorsal, !is.na(gene))
DEGs_day_12_vs_06_dorsal_f <- filter(DEGs_day_12_vs_06_dorsal_f, padj < 0.01, abs(log2FoldChange) >= 2)

DEGs_day_12_vs_06_dorsal$classic_threshold <- ifelse(DEGs_day_12_vs_06_dorsal$padj < 0.01 & abs(DEGs_day_12_vs_06_dorsal$log2FoldChange) >= 1, "yes", "no")
DEGs_day_12_vs_06_dorsal$high_threshold <- ifelse(DEGs_day_12_vs_06_dorsal$padj < 0.01 & abs(DEGs_day_12_vs_06_dorsal$log2FoldChange) >= 2, "yes", "no")
write.csv(DEGs_day_12_vs_06_dorsal, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/DEGs_day_12_vs_06_dorsal.csv")

# DEGs day06 vs day2 for ventral samples
DEGs_day_06_vs_02_ventral <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "ventral")$sample],
    colData = filter(meta, type == "ventral"),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day06", "day02")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day_06_vs_02_ventral$gene <- gene_converter(rownames(DEGs_day_06_vs_02_ventral), "ENSEMBL", "SYMBOL")
DEGs_day_06_vs_02_ventral_f <- filter(DEGs_day_06_vs_02_ventral, !is.na(gene))
DEGs_day_06_vs_02_ventral_f <- filter(DEGs_day_06_vs_02_ventral_f, padj < 0.01, abs(log2FoldChange) >= 2)

DEGs_day_06_vs_02_ventral$classic_threshold <- ifelse(DEGs_day_06_vs_02_ventral$padj < 0.01 & abs(DEGs_day_06_vs_02_ventral$log2FoldChange) >= 1, "yes", "no")
DEGs_day_06_vs_02_ventral$high_threshold <- ifelse(DEGs_day_06_vs_02_ventral$padj < 0.01 & abs(DEGs_day_06_vs_02_ventral$log2FoldChange) >= 2, "yes", "no")
write.csv(DEGs_day_06_vs_02_ventral, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/DEGs_day_06_vs_02_ventral.csv")

# DEGs day06 vs day2 for dorsal samples
DEGs_day_06_vs_02_dorsal <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "dorsal")$sample],
    colData = filter(meta, type == "dorsal"),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day06", "day02")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day_06_vs_02_dorsal$gene <- gene_converter(rownames(DEGs_day_06_vs_02_dorsal), "ENSEMBL", "SYMBOL")
DEGs_day_06_vs_02_dorsal_f <- filter(DEGs_day_06_vs_02_dorsal, !is.na(gene))
DEGs_day_06_vs_02_dorsal_f <- filter(DEGs_day_06_vs_02_dorsal_f, padj < 0.01, abs(log2FoldChange) >= 2)

DEGs_day_06_vs_02_dorsal$classic_threshold <- ifelse(DEGs_day_06_vs_02_dorsal$padj < 0.01 & abs(DEGs_day_06_vs_02_dorsal$log2FoldChange) >= 1, "yes", "no")
DEGs_day_06_vs_02_dorsal$high_threshold <- ifelse(DEGs_day_06_vs_02_dorsal$padj < 0.01 & abs(DEGs_day_06_vs_02_dorsal$log2FoldChange) >= 2, "yes", "no")
write.csv(DEGs_day_06_vs_02_dorsal, "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/DEGs_day_06_vs_02_dorsal.csv")

genes_ventral <- list(
    day_04_02 = filter(DEGs_day_04_vs_02_ventral, abs(log2FoldChange) >= 1 & padj < 0.01 & !is.na(gene))$gene,
    day_06_04 = filter(DEGs_day_06_vs_04_ventral, abs(log2FoldChange) >= 1 & padj < 0.01 & !is.na(gene))$gene,
    day_08_06 = filter(DEGs_day_08_vs_06_ventral, abs(log2FoldChange) >= 1 & padj < 0.01 & !is.na(gene))$gene,
    day_10_08 = filter(DEGs_day_10_vs_08_ventral, abs(log2FoldChange) >= 1 & padj < 0.01 & !is.na(gene))$gene,
    day_12_10 = filter(DEGs_day_12_vs_10_ventral, abs(log2FoldChange) >= 1 & padj < 0.01 & !is.na(gene))$gene,
    day_12_06 = filter(DEGs_day_12_vs_06_ventral, abs(log2FoldChange) >= 1 & padj < 0.01 & !is.na(gene))$gene,
    day_06_02 = filter(DEGs_day_06_vs_02_ventral, abs(log2FoldChange) >= 1 & padj < 0.01 & !is.na(gene))$gene
) %>% Reduce(union, .)
DE_df_ventral <- data.frame(
    gene = genes_ventral,
    day_04_02 = ifelse(
        genes_ventral %in% filter(DEGs_day_04_vs_02_ventral, log2FoldChange >= 1 & padj < 0.01 & !is.na(gene))$gene, "UP",
        ifelse(
            genes_ventral %in% filter(DEGs_day_04_vs_02_ventral, log2FoldChange <= -1 & padj < 0.01 & !is.na(gene))$gene, "DOWN", "no"
        )
    ),
    day_06_04 = ifelse(
        genes_ventral %in% filter(DEGs_day_06_vs_04_ventral, log2FoldChange >= 1 & padj < 0.01 & !is.na(gene))$gene, "UP",
        ifelse(
            genes_ventral %in% filter(DEGs_day_06_vs_04_ventral, log2FoldChange <= -1 & padj < 0.01 & !is.na(gene))$gene, "DOWN", "no"
        )
    ),
    day_08_06 = ifelse(
        genes_ventral %in% filter(DEGs_day_08_vs_06_ventral, log2FoldChange >= 1 & padj < 0.01 & !is.na(gene))$gene, "UP",
        ifelse(
            genes_ventral %in% filter(DEGs_day_08_vs_06_ventral, log2FoldChange <= -1 & padj < 0.01 & !is.na(gene))$gene, "DOWN", "no"
        )
    ),
    day_10_08 = ifelse(
        genes_ventral %in% filter(DEGs_day_10_vs_08_ventral, log2FoldChange >= 1 & padj < 0.01 & !is.na(gene))$gene, "UP",
        ifelse(
            genes_ventral %in% filter(DEGs_day_10_vs_08_ventral, log2FoldChange <= -1 & padj < 0.01 & !is.na(gene))$gene, "DOWN", "no"
        )
    ),
    day_12_10 = ifelse(
        genes_ventral %in% filter(DEGs_day_12_vs_10_ventral, log2FoldChange >= 1 & padj < 0.01 & !is.na(gene))$gene, "UP",
        ifelse(
            genes_ventral %in% filter(DEGs_day_12_vs_10_ventral, log2FoldChange <= -1 & padj < 0.01 & !is.na(gene))$gene, "DOWN", "no"
        )
    ),
    day_12_06 = ifelse(
        genes_ventral %in% filter(DEGs_day_12_vs_06_ventral, log2FoldChange >= 1 & padj < 0.01 & !is.na(gene))$gene, "UP",
        ifelse(
            genes_ventral %in% filter(DEGs_day_12_vs_06_ventral, log2FoldChange <= -1 & padj < 0.01 & !is.na(gene))$gene, "DOWN", "no"
        )
    ),
    day_06_02 = ifelse(
        genes_ventral %in% filter(DEGs_day_06_vs_02_ventral, log2FoldChange >= 1 & padj < 0.01 & !is.na(gene))$gene, "UP",
        ifelse(
            genes_ventral %in% filter(DEGs_day_06_vs_02_ventral, log2FoldChange <= -1 & padj < 0.01 & !is.na(gene))$gene, "DOWN", "no"
        )
    )
)
write.csv(DE_df_ventral, file = "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/Highly_DEG_dbd_ventral.csv", quote = FALSE)

genes_dorsal <- list(
    day_04_02 = filter(DEGs_day_04_vs_02_dorsal, abs(log2FoldChange) >= 1 & padj < 0.01 & !is.na(gene))$gene,
    day_06_04 = filter(DEGs_day_06_vs_04_dorsal, abs(log2FoldChange) >= 1 & padj < 0.01 & !is.na(gene))$gene,
    day_08_06 = filter(DEGs_day_08_vs_06_dorsal, abs(log2FoldChange) >= 1 & padj < 0.01 & !is.na(gene))$gene,
    day_10_08 = filter(DEGs_day_10_vs_08_dorsal, abs(log2FoldChange) >= 1 & padj < 0.01 & !is.na(gene))$gene,
    day_12_10 = filter(DEGs_day_12_vs_10_dorsal, abs(log2FoldChange) >= 1 & padj < 0.01 & !is.na(gene))$gene,
    day_12_06 = filter(DEGs_day_12_vs_06_dorsal, abs(log2FoldChange) >= 1 & padj < 0.01 & !is.na(gene))$gene,
    day_06_02 = filter(DEGs_day_06_vs_02_dorsal, abs(log2FoldChange) >= 1 & padj < 0.01 & !is.na(gene))$gene
) %>% Reduce(union, .)

DE_df_dorsal <- data.frame(
    gene = genes_dorsal,
    day_04_02 = ifelse(
        genes_dorsal %in% filter(DEGs_day_04_vs_02_dorsal, log2FoldChange >= 1 & padj < 0.01 & !is.na(gene))$gene, "UP",
        ifelse(
            genes_dorsal %in% filter(DEGs_day_04_vs_02_dorsal, log2FoldChange <= -1 & padj < 0.01 & !is.na(gene))$gene, "DOWN", "no"
        )
    ),
    day_06_04 = ifelse(
        genes_dorsal %in% filter(DEGs_day_06_vs_04_dorsal, log2FoldChange >= 1 & padj < 0.01 & !is.na(gene))$gene, "UP",
        ifelse(
            genes_dorsal %in% filter(DEGs_day_06_vs_04_dorsal, log2FoldChange <= -1 & padj < 0.01 & !is.na(gene))$gene, "DOWN", "no"
        )
    ),
    day_08_06 = ifelse(
        genes_dorsal %in% filter(DEGs_day_08_vs_06_dorsal, log2FoldChange >= 1 & padj < 0.01 & !is.na(gene))$gene, "UP",
        ifelse(
            genes_dorsal %in% filter(DEGs_day_08_vs_06_dorsal, log2FoldChange <= -1 & padj < 0.01 & !is.na(gene))$gene, "DOWN", "no"
        )
    ),
    day_10_08 = ifelse(
        genes_dorsal %in% filter(DEGs_day_10_vs_08_dorsal, log2FoldChange >= 1 & padj < 0.01 & !is.na(gene))$gene, "UP",
        ifelse(
            genes_dorsal %in% filter(DEGs_day_10_vs_08_dorsal, log2FoldChange <= -1 & padj < 0.01 & !is.na(gene))$gene, "DOWN", "no"
        )
    ),
    day_12_10 = ifelse(
        genes_dorsal %in% filter(DEGs_day_12_vs_10_dorsal, log2FoldChange >= 1 & padj < 0.01 & !is.na(gene))$gene, "UP",
        ifelse(
            genes_dorsal %in% filter(DEGs_day_12_vs_10_dorsal, log2FoldChange <= -1 & padj < 0.01 & !is.na(gene))$gene, "DOWN", "no"
        )
    ),
    day_12_06 = ifelse(
        genes_dorsal %in% filter(DEGs_day_12_vs_06_dorsal, log2FoldChange >= 1 & padj < 0.01 & !is.na(gene))$gene, "UP",
        ifelse(
            genes_dorsal %in% filter(DEGs_day_12_vs_06_dorsal, log2FoldChange <= -1 & padj < 0.01 & !is.na(gene))$gene, "DOWN", "no"
        )
    ),
    day_06_02 = ifelse(
        genes_dorsal %in% filter(DEGs_day_06_vs_02_dorsal, log2FoldChange >= 1 & padj < 0.01 & !is.na(gene))$gene, "UP",
        ifelse(
            genes_dorsal %in% filter(DEGs_day_06_vs_02_dorsal, log2FoldChange <= -1 & padj < 0.01 & !is.na(gene))$gene, "DOWN", "no"
        )
    )
)
write.csv(DE_df_dorsal, file = "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/Highly_DEG_dbd_dorsal.csv", quote = FALSE)
