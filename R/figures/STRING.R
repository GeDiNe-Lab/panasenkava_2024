# Loading packages and functions
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(DESeq2)
library(Seurat)
library(reshape2)
library(networkD3)

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


##############
##############
# Formating STRINGdb network outpout for Cytoscape


WGCNA_genes <- read.csv("results/tables/Figure_4/SHH_cluster.csv", header = TRUE)
rownames(WGCNA_genes) <- WGCNA_genes$gene
View(WGCNA_genes)
figure_genes <- read.csv("results/tables/Figure_genelist.csv", header = TRUE)
#  Got a doublet symbol for CKMT1B
figure_genes$gene[which(figure_genes$X == "ENSG00000223572")] <- "CKMT1B_b"
figure_genes <- filter(figure_genes, !is.na(gene))
rownames(figure_genes) <- figure_genes$gene

STRING_res <- read.table("results/tables/Figure_4/string_interactions_short.tsv", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")




STRING_formated <- STRING_res %>% dplyr::select(X.node1, node2, combined_score)
colnames(STRING_formated) <- c("gene1", "gene2", "STRING_score")
STRING_formated$double <- rep("keep", nrow(STRING_formated))
STRING_formated$double
STRING_formated_rev <- data.frame(gene1 = STRING_formated$gene2, gene2 = STRING_formated$gene1, STRING_score = STRING_formated$STRING_score, double = rep("exclude", nrow(STRING_formated)))
STRING_formated <- rbind(STRING_formated, STRING_formated_rev)

# Retaking mouse genes without any STRING_connection (they are linked to themselves) in the network
retake_genes <- setdiff(filter(figure_genes, mouse == "tested")$gene, c(STRING_formated$gene1, STRING_formated$gene2) %>% unique()) %>%
    lapply(function(gene) {
        return(data.frame(gene1 = gene, gene2 = gene, STRING_score = 0.4, double = "exclude"))
    }) %>%
    Reduce(rbind, .)
STRING_formated <- rbind(STRING_formated, retake_genes)

STRING_formated$gene_mouse <- STRING_formated$gene1 %>% sapply(function(gene) {
    return(ifelse(gene == "SHH", "SHH", figure_genes[gene, "mouse"]))
})

STRING_formated$scRNA <- STRING_formated$gene1 %>% sapply(function(gene) {
    return(figure_genes[gene, "sc_w3_w4"])
})

STRING_formated$gene_SHH_link <- c(1:nrow(STRING_formated)) %>% sapply(function(i) {
    return(ifelse(STRING_formated$gene2[i] == "SHH", TRUE, FALSE))
})

# to modify
mouse_genes_cor <- read.csv("results/tables/Figure_4/mouse_genes_cor.csv", header = TRUE)
rownames(mouse_genes_cor) <- mouse_genes_cor$X %>% gene_converter("ENSEMBL", "SYMBOL")

STRING_formated$gene1_SHH_corr <- STRING_formated$gene1 %>% sapply(function(gene) {
    if (gene %in% rownames(mouse_genes_cor)) {
        return(mouse_genes_cor[gene, "cor"])
    } else {
        return(WGCNA_genes[gene, "cor"])
    }
})

STRING_formated$gene2_SHH_corr <- STRING_formated$gene2 %>% sapply(function(gene) {
    if (gene %in% rownames(mouse_genes_cor)) {
        return(mouse_genes_cor[gene, "cor"])
    } else {
        return(WGCNA_genes[gene, "cor"])
    }
})



write.csv(
    filter(
        STRING_formated,
        (gene1_SHH_corr > 0 & gene2_SHH_corr > 0)
    ),
    file = "results/tables/Figure_4/STRING_formated_ventral.csv", row.names = FALSE, quote = FALSE
)

write.csv(
    filter(
        STRING_formated,
        (
            (gene1_SHH_corr < 0 & gene2_SHH_corr < 0) |
                (gene1_SHH_corr < 0 & gene2 == "SHH") |
                (gene1 == "SHH" & gene2_SHH_corr < 0)
        )
    ),
    file = "results/tables/Figure_4/STRING_formated_dorsal.csv", row.names = FALSE, quote = FALSE
)
