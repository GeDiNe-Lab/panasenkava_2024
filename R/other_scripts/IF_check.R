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
