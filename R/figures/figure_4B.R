# Loading packages and functions
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(DESeq2)
library(Seurat)
library(reshape2)

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
authors_meta$celltype_region %>% unique()
authors_meta_f <- filter(
    authors_meta,
    (celltype_region == "FB" & week_stage == "W3-1") |
        (celltype_region %in% c("Brain Neuron", "FB") & week_stage == "W4-1") |
        (celltype_region == "Brain Neuron" & week_stage == "W5-1")
)

######################
######################
# Percentages of cells expresing WGCNA modules genes in human single cell data

SHH_cluster <- read.csv("results/tables/Figure_4/SHH_cluster.csv", header = TRUE)

setdiff(SHH_cluster$gene, rownames(sc_counts))
genes <- intersect(rownames(sc_counts), SHH_cluster$gene)
gene_long <- melt(as.matrix(sc_counts[genes, authors_meta_f$barcode]))

colnames(gene_long) <- c("gene", "barcode", "expression")

# Merge with week_stage information
merged_data <- merge(gene_long, dplyr::select(authors_meta_f, c("barcode", "week_stage")), by = "barcode")

# Calculate percentage of cells with expression > 0 per week_stage and gene
week_gene_df <- merged_data %>%
    group_by(week_stage, gene) %>%
    summarise(percent_expressed = mean(expression > 0) * 100, .groups = "drop") %>%
    tidyr::spread(key = week_stage, value = percent_expressed) %>%
    as.data.frame()

write.csv(week_gene_df, "results/tables/Figure_4/sc_percent.csv", row.names = FALSE)

######################
######################
# Figure 4C : Expression in human single cell for key genes from WGCNA

#  SPON1 just not in it
ventral_genes <- c(
    "GLI1",
    "FOXA2",
    "NKX2-2",
    "NKX2-1",
    "HEY1",
    "SHH",
    "FREM1",
    "PTCH1",
    "NES",
    "IFT88",
    "CAMK2B",
    "FGFR3",
    "PTCH2",
    "NCAM1",
    "FOXB1"
)
dorsal_genes <- c(
    "PAX6",
    "CDON",
    "ZIC2",
    "GLI2",
    "GLI3",
    "LHX2",
    "FGF2",
    "BOC",
    "ZIC1",
    "ZIC3",
    "OTX2",
    "SP8",
    "FOXG1",
    "FOS",
    "GAS1",
    "SOX1",
    "KLF7",
    "SMOC1"
)
selected_genes <- c(ventral_genes, dorsal_genes)
sc_counts_genes <- sc_counts[selected_genes[selected_genes %in% rownames(sc_counts)], authors_meta_f$barcode]
sc_meta_genes <- authors_meta_f %>% filter(barcode %in% colnames(sc_counts_genes))
sc_meta_genes <- cbind(sc_meta_genes, t(sc_counts_genes))

summary_by_week_region <- sc_meta_genes %>%
    group_by(week_stage) %>%
    summarise(across(all_of(selected_genes),
        list(
            expression = ~ mean(.),
            percent = ~ mean(. > 0) * 100
        ),
        .names = "{.col}_{.fn}"
    )) %>%
    pivot_longer(
        cols = -week_stage,
        names_to = c("gene", ".value"),
        names_sep = "_"
    ) %>%
    as.data.frame()

ventral_genes <- ggplot(
    filter(summary_by_week_region, gene %in% ventral_genes & percent > 0),
    aes(x = week_stage, y = gene, size = percent, color = expression)
) +
    geom_point() +
    scale_size_continuous(range = c(1, 10), limits = c(0, 100)) + # Ensure percent scale from 0 to 100
    scale_color_gradientn(colors = c("grey", "#cd7dff", "#312eff")) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
    ) +
    labs(size = "Percent\nexpresed", color = "Average\nexpression") # Add legend labels

dorsal_genes <- ggplot(
    filter(summary_by_week_region, gene %in% dorsal_genes & percent > 0),
    aes(x = week_stage, y = gene, size = percent, color = expression)
) +
    geom_point() +
    scale_size_continuous(range = c(1, 10), limits = c(0, 100)) + # Ensure percent scale from 0 to 100
    scale_color_gradientn(colors = c("grey", "#cd7dff", "#312eff")) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
    ) +
    labs(size = "Percent\nexpresed", color = "Average\nexpression") # Add legend labels

ggsave("results/images/Figure_4/Figure4C_v.png", ventral_genes, width = 3.5, height = 5.5)
ggsave("results/images/Figure_4/Figure4C_d.png", dorsal_genes, width = 3.5, height = 6.5)
