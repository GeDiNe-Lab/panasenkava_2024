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
C_Het_DE <- filter(C_Het, abs(log2FoldChange) >= 1, padj < 0.01, !is.na(gene))$gene
C_Hom_DE <- filter(C_Hom, abs(log2FoldChange) >= 1, padj < 0.01, !is.na(gene))$gene
Het_Hom_DE <- filter(Het_Hom, abs(log2FoldChange) >= 1, padj < 0.01, !is.na(gene))$gene

# Cyclopamin DEGs
High_No_DE <- filter(cyclopamin, HvsN_thr == TRUE, !is.na(genes))$genes
High_Low_DE <- filter(cyclopamin, HvsL_thr == TRUE, !is.na(genes))$genes
Low_No_DE <- filter(cyclopamin, LvsN_thr == TRUE, !is.na(genes))$genes

# C_Het_DE vs High_No_DE
venn.diagram(
    x = list(C_Het_DE, High_No_DE),
    disable.logging = TRUE,
    category.names = c("C_Het_DE", "High_No_DE"),
    main = "C_Het_DE    High_No_DE",
    filename = "results/images/Figure_4/venn_C_Het_DE_vs_High_No_DE.png"
)
C_Het_DE_vs_High_No_DE <- data.frame(
    genes = union(C_Het_DE, High_No_DE)
)
C_Het_DE_vs_High_No_DE$state <- C_Het_DE_vs_High_No_DE$genes %>% sapply(function(gene) {
    if (gene %in% C_Het_DE & gene %in% High_No_DE) {
        return("common")
    } else if (gene %in% C_Het_DE) {
        return("C_Het_DE")
    } else {
        return("High_No_DE")
    }
})
write.csv(C_Het_DE_vs_High_No_DE, "results/tables/Figure_4/venn_C_Het_DE_vs_High_No_DE.csv")


# C_Het_DE vs High_Low_DE
venn.diagram(
    x = list(C_Het_DE, High_Low_DE),
    category.names = c("C_Het_DE", "High_Low_DE"),
    disable.logging = TRUE,
    main = "C_Het_DE    High_Low_DE",
    filename = "results/images/Figure_4/venn_C_Het_DE_vs_High_Low_DE.png"
)
C_Het_DE_vs_High_Low_DE <- data.frame(
    genes = union(C_Het_DE, High_Low_DE)
)
C_Het_DE_vs_High_Low_DE$state <- C_Het_DE_vs_High_Low_DE$genes %>% sapply(function(gene) {
    if (gene %in% C_Het_DE & gene %in% High_Low_DE) {
        return("common")
    } else if (gene %in% C_Het_DE) {
        return("C_Het_DE")
    } else {
        return("High_Low_DE")
    }
})
write.csv(C_Het_DE_vs_High_Low_DE, "results/tables/Figure_4/venn_C_Het_DE_vs_High_Low_DE.csv")


# C_Het_DE vs Low_No_DE
venn.diagram(
    x = list(C_Het_DE, Low_No_DE),
    category.names = c("C_Het_DE", "Low_No_DE"),
    disable.logging = TRUE,
    main = "C_Het_DE    Low_No_DE",
    filename = "results/images/Figure_4/venn_C_Het_DE_vs_Low_No_DE.png"
)
C_Het_DE_vs_Low_No_DE <- data.frame(
    genes = union(C_Het_DE, Low_No_DE)
)
C_Het_DE_vs_Low_No_DE$state <- C_Het_DE_vs_Low_No_DE$genes %>% sapply(function(gene) {
    if (gene %in% C_Het_DE & gene %in% Low_No_DE) {
        return("common")
    } else if (gene %in% C_Het_DE) {
        return("C_Het_DE")
    } else {
        return("Low_No_DE")
    }
})
write.csv(C_Het_DE_vs_Low_No_DE, "results/tables/Figure_4/venn_C_Het_DE_vs_Low_No_DE.csv")


# C_Hom_DE vs High_No_DE
venn.diagram(
    x = list(C_Hom_DE, High_No_DE),
    category.names = c("C_Hom_DE", "High_No_DE"),
    disable.logging = TRUE,
    main = "C_Hom_DE    High_No_DE",
    filename = "results/images/Figure_4/venn_C_Hom_DE_vs_High_No_DE.png"
)
C_Hom_DE_vs_High_No_DE <- data.frame(
    genes = union(C_Hom_DE, High_No_DE)
)
C_Hom_DE_vs_High_No_DE$state <- C_Hom_DE_vs_High_No_DE$genes %>% sapply(function(gene) {
    if (gene %in% C_Hom_DE & gene %in% High_No_DE) {
        return("common")
    } else if (gene %in% C_Hom_DE) {
        return("C_Hom_DE")
    } else {
        return("High_No_DE")
    }
})
write.csv(C_Hom_DE_vs_High_No_DE, "results/tables/Figure_4/venn_C_Hom_DE_vs_High_No_DE.csv")


# C_Hom_DE vs High_Low_DE
venn.diagram(
    x = list(C_Hom_DE, High_Low_DE),
    category.names = c("C_Hom_DE", "High_Low_DE"),
    disable.logging = TRUE,
    main = "C_Hom_DE    High_Low_DE",
    filename = "results/images/Figure_4/venn_C_Hom_DE_vs_High_Low_DE.png"
)
C_Hom_DE_vs_High_Low_DE <- data.frame(
    genes = union(C_Hom_DE, High_Low_DE)
)
C_Hom_DE_vs_High_Low_DE$state <- C_Hom_DE_vs_High_Low_DE$genes %>% sapply(function(gene) {
    if (gene %in% C_Hom_DE & gene %in% High_Low_DE) {
        return("common")
    } else if (gene %in% C_Hom_DE) {
        return("C_Hom_DE")
    } else {
        return("High_Low_DE")
    }
})
write.csv(C_Hom_DE_vs_High_Low_DE, "results/tables/Figure_4/venn_C_Hom_DE_vs_High_Low_DE.csv")


# C_Hom_DE vs Low_No_DE
venn.diagram(
    x = list(C_Hom_DE, Low_No_DE),
    category.names = c("C_Hom_DE", "Low_No_DE"),
    disable.logging = TRUE,
    main = "C_Hom_DE    Low_No_DE",
    filename = "results/images/Figure_4/venn_C_Hom_DE_vs_Low_No_DE.png"
)
C_Hom_DE_vs_Low_No_DE <- data.frame(
    genes = union(C_Hom_DE, Low_No_DE)
)
C_Hom_DE_vs_Low_No_DE$state <- C_Hom_DE_vs_Low_No_DE$genes %>% sapply(function(gene) {
    if (gene %in% C_Hom_DE & gene %in% Low_No_DE) {
        return("common")
    } else if (gene %in% C_Hom_DE) {
        return("C_Hom_DE")
    } else {
        return("Low_No_DE")
    }
})
write.csv(C_Hom_DE_vs_Low_No_DE, "results/tables/Figure_4/venn_C_Hom_DE_vs_Low_No_DE.csv")


# Het_Hom_DE vs High_No_DE
venn.diagram(
    x = list(Het_Hom_DE, High_No_DE),
    category.names = c("Het_Hom_DE", "High_No_DE"),
    disable.logging = TRUE,
    main = "Het_Hom_DE    High_No_DE",
    filename = "results/images/Figure_4/venn_Het_Hom_DE_vs_High_No_DE.png"
)
Het_Hom_DE_vs_High_No_DE <- data.frame(
    genes = union(Het_Hom_DE, High_No_DE)
)
Het_Hom_DE_vs_High_No_DE$state <- Het_Hom_DE_vs_High_No_DE$genes %>% sapply(function(gene) {
    if (gene %in% Het_Hom_DE & gene %in% High_No_DE) {
        return("common")
    } else if (gene %in% Het_Hom_DE) {
        return("Het_Hom_DE")
    } else {
        return("High_No_DE")
    }
})
write.csv(Het_Hom_DE_vs_High_No_DE, "results/tables/Figure_4/venn_Het_Hom_DE_vs_High_No_DE.csv")


# Het_Hom_DE vs High_Low_DE
venn.diagram(
    x = list(Het_Hom_DE, High_Low_DE),
    category.names = c("Het_Hom_DE", "High_Low_DE"),
    disable.logging = TRUE,
    main = "Het_Hom_DE    High_Low_DE",
    filename = "results/images/Figure_4/venn_Het_Hom_DE_vs_High_Low_DE.png"
)
Het_Hom_DE_vs_High_Low_DE <- data.frame(
    genes = union(Het_Hom_DE, High_Low_DE)
)
Het_Hom_DE_vs_High_Low_DE$state <- Het_Hom_DE_vs_High_Low_DE$genes %>% sapply(function(gene) {
    if (gene %in% Het_Hom_DE & gene %in% High_Low_DE) {
        return("common")
    } else if (gene %in% Het_Hom_DE) {
        return("Het_Hom_DE")
    } else {
        return("High_Low_DE")
    }
})
write.csv(Het_Hom_DE_vs_High_Low_DE, "results/tables/Figure_4/venn_Het_Hom_DE_vs_High_Low_DE.csv")

# Het_Hom_DE vs Low_No_DE
venn.diagram(
    x = list(Het_Hom_DE, Low_No_DE),
    category.names = c("Het_Hom_DE", "Low_No_DE"),
    disable.logging = TRUE,
    main = "Het_Hom_DE    Low_No_DE",
    filename = "results/images/Figure_4/venn_Het_Hom_DE_vs_Low_No_DE.png"
)
Het_Hom_DE_vs_Low_No_DE <- data.frame(
    genes = union(Het_Hom_DE, Low_No_DE)
)
Het_Hom_DE_vs_Low_No_DE$state <- Het_Hom_DE_vs_Low_No_DE$genes %>% sapply(function(gene) {
    if (gene %in% Het_Hom_DE & gene %in% Low_No_DE) {
        return("common")
    } else if (gene %in% Het_Hom_DE) {
        return("Het_Hom_DE")
    } else {
        return("Low_No_DE")
    }
})
write.csv(Het_Hom_DE_vs_Low_No_DE, "results/tables/Figure_4/venn_Het_Hom_DE_vs_Low_No_DE.csv")


# missing barplots (genes tested in mice not in WGCNA cluster)
