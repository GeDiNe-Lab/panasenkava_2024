# Loading packages and functions
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(DESeq2)

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
meta <- filter(meta, sample != "LON71_D12_2", type %in% c("ventral", "dorsal"))
counts <- rawcounts[which(rowSums(rawcounts) >= 50), meta$sample]

dds <- DESeqDataSetFromMatrix(
    countData = counts[, meta$sample],
    colData = meta,
    design = ~ type + day
)

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# PCA plot
pca.data <- plotPCA(vsd, intgroup = c("type", "day", "line"), returnData = TRUE)
percentVar <- round(100 * attr(pca.data, "percentVar"))

png(filename = "results/images/Figure_2A/F2A_1_PCA_type.png", width = 1600, height = 1200, res = 250)
ggplot(pca.data, aes(PC1, PC2, color = type, shape = day)) +
    geom_point(size = 2, stroke = 1) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_color_manual(values = c("#A1A1DE", "#80AD3C")) +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6)) +
    custom_theme() +
    ggtitle("PCA of dorsal and ventral kinetics")
dev.off()

png(filename = "results/images/Figure_2A/F2A_1_PCA_line.png", width = 1600, height = 1200, res = 250)
ggplot(pca.data, aes(PC1, PC2, color = line, shape = day)) +
    geom_point(size = 2, stroke = 1) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6)) +
    custom_theme() +
    ggtitle("PCA of dorsal and ventral kinetics")
dev.off()

norm <- vsd@assays@data[[1]]

DEGs_DV <- DESeq(dds) %>%
    results(alpha = 0.05, contrast = c("type", "ventral", "dorsal")) %>%
    as.data.frame() %>%
    na.omit()

DEGs_DV_f <- filter(DEGs_DV, padj < 0.05, abs(log2FoldChange) >= 1)
DEGs_DV_f$gene <- gene_converter(rownames(DEGs_DV_f), "ENSEMBL", "SYMBOL")
DEGs_DV_f <- filter(DEGs_DV_f, !is.na(gene))


vstnorm <- vst(dds, blind = FALSE)
mat <- assay(vstnorm)[rownames(DEGs_DV_f), ]
# scaling the matrix
scaled_mat <- t(apply(mat, 1, scale))
colnames(scaled_mat) <- colnames(mat)

# png(filename = "results/images/F1_3_DE_HM.png", width = 1600, height = 1600, res = 250)


sample_order <- c(
    filter(meta, type == "dorsal")$sample[order(filter(meta, type == "dorsal")$day, decreasing = TRUE)],
    filter(meta, type == "ventral")$sample[order(filter(meta, type == "ventral")$day)]
)

Heatmap(
    scaled_mat[, sample_order],
    name = "Normalized expression",
    column_names_gp = gpar(fontsize = 6),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    row_names_side = "left",
    show_column_names = TRUE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    width = ncol(scaled_mat) * unit(2, "mm"),
    # height = nrow(mat) * unit(5, "mm"),
    col = colorRampPalette(c(
        "blue",
        "white",
        "red"
    ))(1000),
)
# dev.off()


norm <- assay(vsd)
rownames(norm) <- rownames(norm) %>% gene_converter("ENSEMBL", "SYMBOL")
norm <- norm[!is.na(rownames(norm)), ]

IF <- c(
    "PAX6",
    "TBR1",
    "EOMES",
    "EMX2",
    "NKX2-1",
    "SHH",
    "FOXA2",
    "FOXA1",
    "MEIS2",
    "OLIG2",
    "OTX2"
)

meta_IF <- cbind(meta, t(norm[IF, ]))

meta_IF <- meta_IF[order(meta_IF$day), ]
meta_IF %>% colnames()

norm %>% min()
norm %>% max()
meta_IF_merged <- meta_IF %>%
    group_by(type, line, day, sexe) %>%
    summarise(
        PAX6 = mean(PAX6, na.rm = TRUE),
        TBR1 = mean(TBR1, na.rm = TRUE),
        EOMES = mean(EOMES, na.rm = TRUE),
        EMX2 = mean(EMX2, na.rm = TRUE),
        `NKX2-1` = mean(`NKX2-1`, na.rm = TRUE),
        SHH = mean(SHH, na.rm = TRUE),
        FOXA2 = mean(FOXA2, na.rm = TRUE),
        FOXA1 = mean(FOXA1, na.rm = TRUE),
        MEIS2 = mean(MEIS2, na.rm = TRUE),
        OLIG2 = mean(OLIG2, na.rm = TRUE),
        OTX2 = mean(OTX2, na.rm = TRUE),
    )
meta_IF_merged
meta_IF_merged$type_line <- paste(meta_IF_merged$type, meta_IF_merged$line, sep = "_")
filter(meta_IF_merged, type == "ventral")

for (gene in IF) {
    print(paste0("results/images/check_IF/", gene, ".png"))
    meta_IF_merged$gene <- meta_IF_merged[[gene]]
    ggplot(data = meta_IF_merged, aes(x = day, y = gene, color = line, linetype = type, group = type_line)) +
        geom_point() +
        geom_line(stat = "identity") +
        ylab(gene) +
        geom_hline(yintercept = max(norm[IF, ]), linetype = "dotted") +
        geom_hline(yintercept = max(norm), linetype = "dashed") +
        ylim(min(norm), max(norm)) +
        custom_theme()
    ggsave(paste0("results/images/check_IF/", gene, ".png"))
}
