# Loading packages and functions
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(DESeq2)
library(Seurat)

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
# Loading data (path to change later)
rawcounts <- readcounts("data/rawcounts.csv", sep = ",", header = TRUE)
rawmeta <- read.table("data/meta.csv", sep = ",", header = TRUE)

# LON71_D12_2 does not have any reads in the count file
# though, the fastQC report shows that the sample is good
meta <- filter(rawmeta, sample != "LON71_D12_2", diff == "diff13", line %in% c("LON71", "WTC"))
View(meta)
# filtering out lowly expressed genes
counts <- rawcounts[, meta$sample][which(rowSums(rawcounts[, meta$sample]) >= 25), ]

# making DESeq object with lineage,days and type as covariates
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ line + day + type
)

# Normalization by variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)
rownames(vsd) <- gene_converter(rownames(vsd), "ENSEMBL", "SYMBOL")
rownames(vsd)
# loading single cell data from Zeng et al from week3,week4 and week5
sc_counts <- readMM("scdata/week345_wholebody.mtx") %>% t()

# keeping only cells related to Central Nervous System (SNC)
cell_ids <- read.table("scdata/indices_week345_wholebody.csv", header = TRUE, sep = ",")
genes <- read.table("scdata/genes_week345_wholebody.tsv", header = FALSE)$V1


rownames(sc_counts) <- genes

formated_ids <- c(1:nrow(cell_ids)) %>% sapply(function(i) {
    c(str_split(cell_ids[i, ]$index, "-")[[1]][1:2], cell_ids[i, ]$week_stage) %>%
        paste(collapse = "-") %>%
        return()
})
colnames(sc_counts) <- formated_ids

authors_meta <- read.csv("scdata/week345_wholebody_metadata.csv")
authors_meta_f <- filter(authors_meta, celltype_region_num %in% c(1:8))


sc_counts_f_int <- sc_counts[, authors_meta_f$barcode]
gene_filter <- apply(sc_counts_f_int, 1, function(row) sum(row > 0) >= 20)
sc_counts_f <- sc_counts_f_int[names(which(gene_filter == TRUE)), ]

rm(sc_counts_f_int)
rm(sc_counts)


genes <- c("FGF10", "DDC", "LRRK2", "PLCL1", "DRC1", "GADL1", "SLIT2", "NKX2-2", "NKX2-1", "SHH", "LRP2", "FRZB", "FOXA2", "PTCH1", "CAPN6", "SIX3", "LINC00261")
other_genes <- c("CAPN6", "PLCL1", "FRZB", "LINC00261", "DDC", "SLIT2", "LRRK2")
other_genes %in% rownames(sc_counts_f)

sc_violin <- sc_counts_f[genes, ] %>%
    t() %>%
    as.data.frame()
sc_violin$week <- authors_meta_f$week_stage
sc_violin$celltype <- authors_meta_f$celltype_region

sc_violin_long <- sc_violin %>%
    pivot_longer(cols = -c(week, celltype), names_to = "gene", values_to = "expression") %>%
    as.data.frame()

table(authors_meta_f$celltype_region, authors_meta_f$week_stage)

sc_violin_long$gene <- factor(sc_violin_long$gene, levels = genes)
ggplot(filter(sc_violin_long, expression != 0, celltype %in% c("FB", "Brain Neuron")), aes(x = expression, y = gene, fill = gene)) +
    geom_boxplot() +
    xlim(0, 30) +
    theme_minimal() +
    theme(
        legend.position = "none",
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        strip.text = element_text(size = 20, face = "bold"),
    ) +
    facet_grid(~week, scales = "fixed")

seurat_obj <- CreateSeuratObject(counts = sc_counts_f, meta.data = authors_meta_f)
seurat_obj <- NormalizeData(seurat_obj)
Idents(seurat_obj) <- "week_stage"

DE_weeks <- FindAllMarkers(seurat_obj, logfc.threshold = 2) %>% filter(p_val_adj < 0.01)
DE_weeks$abs_logFC <- abs(DE_weeks$avg_log2FC)
DE_weeks$diff_pct <- DE_weeks$pct.1 - DE_weeks$pct.2


DE_genes <- list(
    filter(DE_weeks, cluster == "W3-1")$gene[order(filter(DE_weeks, cluster == "W3-1")$diff_pct, decreasing = TRUE)] %>% head(20),
    filter(DE_weeks, cluster == "W4-1")$gene[order(filter(DE_weeks, cluster == "W3-1")$diff_pct, decreasing = TRUE)] %>% head(20),
    filter(DE_weeks, cluster == "W5-1")$gene[order(filter(DE_weeks, cluster == "W3-1")$diff_pct, decreasing = TRUE)] %>% head(20)
) %>%
    unlist() %>%
    unique()
DE_genes












DE_weeks %>%
    filter(cluster == "W4-1") %>%
    View()
DE_weeks %>%
    filter(cluster == "W5-1") %>%
    View()

sc_norm <- seurat_obj@assays$RNA$data

common_genes <- intersect(rownames(vsd), rownames(sc_norm))

scDEGs <- read.csv("scdata/DEGs_week345.csv")
scDEGs <- filter(scDEGs, celltype_region %in% c(
    "Spinal Cord Motor Neuron 1",
    "Spinal Cord Motor Neuron 2",
    "Spinal Cord Motor Neuron 3",
    "Forebrain progenitor",
    "Midbrain progenitor",
    "Hindbrain progenitor",
    "Spinal cord progenitor",
    "Inhibitory Neuron"
) & Gene %in% common_genes)
scDEGs %>% View()
scDEGs$celltype_region %>% table()

genes1 <- filter(scDEGs, celltype_region %in% c("Forebrain progenitor"))$Gene
genes2 <- filter(scDEGs, celltype_region %in% c("Midbrain progenitor"))$Gene

test1 <- sc_norm[genes2, filter(authors_meta_f, week_stage == "W3-1" & celltype_region == "MB")$barcode]
testrank1 <- percent_rank(rowSums(test1))

test2 <- assay(vsd)[genes2, filter(meta, type == "ventral")$sample]
testrank2 <- percent_rank(rowSums(test2))

cor(testrank1, testrank2, method = "pearson")
