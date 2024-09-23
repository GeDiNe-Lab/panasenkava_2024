# Loading packages and functions
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(DESeq2)

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

# Keeping only necessary samples
meta <- filter(rawmeta, type %in% c("dorsal", "ventral"), diff %in% c("diff9", "diff12"), CRISPR %in% c("no", "control"), sample != "L9C1_2")

C_Het <- read.csv("results/tables/Figure_4/DE_vAN_control_vs_het.csv")
C_Hom <- read.csv("results/tables/Figure_4/DE_vAN_control_vs_homo.csv")
Het_Hom <- read.csv("results/tables/Figure_4/DE_vAN_het_vs_homo.csv")

load("results/tables/Figure_3/cyclo_genes_df.RData")
cyclopamin <- cyclo_genes_df
View(cyclopamin)
View(C_Het)
# Loading the VennDiagram package
library(VennDiagram)

# CRISPR DEGs
C_Het_DE <- filter(C_Het, log2FoldChange >= 1, padj < 0.01, !is.na(gene))$gene
C_Hom_DE <- filter(C_Hom, log2FoldChange >= 1, padj < 0.01, !is.na(gene))$gene
Het_Hom_DE <- filter(Het_Hom, log2FoldChange >= 1, padj < 0.01, !is.na(gene))$gene

# Cyclopamin DEGs
High_No_DE <- filter(cyclopamin, HvsN_thr == TRUE, !is.na(genes))$genes
High_Low_DE <- filter(cyclopamin, HvsL_thr == TRUE, !is.na(genes))$genes
Low_No_DE <- filter(cyclopamin, LvsN_thr == TRUE, !is.na(genes))$genes

venn.diagram(
    x = list(C_Het_DE, High_No_DE),
    disable.logging = TRUE,
    category.names = c("C_Het_DE", "High_No_DE"),
    main = "C_Het_DE    High_No_DE",
    filename = "results/venn_C_Het_DE_vs_High_No_DE.png"
)

venn.diagram(
    x = list(C_Het_DE, High_Low_DE),
    category.names = c("C_Het_DE", "High_Low_DE"),
    disable.logging = TRUE,
    main = "C_Het_DE    High_Low_DE",
    filename = "results/venn_C_Het_DE_vs_High_Low_DE.png"
)

venn.diagram(
    x = list(C_Het_DE, Low_No_DE),
    category.names = c("C_Het_DE", "Low_No_DE"),
    disable.logging = TRUE,
    main = "C_Het_DE    Low_No_DE",
    filename = "results/venn_C_Het_DE_vs_Low_No_DE.png"
)

# C_Hom_DE vs cyclopamin DEGs
venn.diagram(
    x = list(C_Hom_DE, High_No_DE),
    category.names = c("C_Hom_DE", "High_No_DE"),
    disable.logging = TRUE,
    main = "C_Hom_DE    High_No_DE",
    filename = "results/venn_C_Hom_DE_vs_High_No_DE.png"
)

venn.diagram(
    x = list(C_Hom_DE, High_Low_DE),
    category.names = c("C_Hom_DE", "High_Low_DE"),
    disable.logging = TRUE,
    main = "C_Hom_DE    High_Low_DE",
    filename = "results/venn_C_Hom_DE_vs_High_Low_DE.png"
)

venn.diagram(
    x = list(C_Hom_DE, Low_No_DE),
    category.names = c("C_Hom_DE", "Low_No_DE"),
    disable.logging = TRUE,
    main = "C_Hom_DE    Low_No_DE",
    filename = "results/venn_C_Hom_DE_vs_Low_No_DE.png"
)

# Het_Hom_DE vs cyclopamin DEGs
venn.diagram(
    x = list(Het_Hom_DE, High_No_DE),
    category.names = c("Het_Hom_DE", "High_No_DE"),
    disable.logging = TRUE,
    main = "Het_Hom_DE    High_No_DE",
    filename = "results/venn_Het_Hom_DE_vs_High_No_DE.png"
)

venn.diagram(
    x = list(Het_Hom_DE, High_Low_DE),
    category.names = c("Het_Hom_DE", "High_Low_DE"),
    disable.logging = TRUE,
    main = "Het_Hom_DE    High_Low_DE",
    filename = "results/venn_Het_Hom_DE_vs_High_Low_DE.png"
)

venn.diagram(
    x = list(Het_Hom_DE, Low_No_DE),
    category.names = c("Het_Hom_DE", "Low_No_DE"),
    disable.logging = TRUE,
    main = "Het_Hom_DE    Low_No_DE",
    filename = "results/venn_Het_Hom_DE_vs_Low_No_DE.png"
)
