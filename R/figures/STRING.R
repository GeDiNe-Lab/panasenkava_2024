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

figure_genes <- read.csv("results/tables/Figure_genelist.csv", header = TRUE)
#  Got a doublet symbol for CKMT1B
figure_genes$gene[which(figure_genes$X == "ENSG00000223572")] <- "CKMT1B_b"
figure_genes <- filter(figure_genes, !is.na(gene))
rownames(figure_genes) <- figure_genes$gene

STRING_res <- read.table("results/tables/Figure_4/STRING_short_tab_output.tsv", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")

c(STRING_res$X.node1, STRING_res$node2) %>%
    unique() %>%
    length()


STRING_formated <- STRING_res %>% dplyr::select(X.node1, node2, combined_score)
colnames(STRING_formated) <- c("gene1", "gene2", "STRING_score")
STRING_formated$double <- rep("keep", nrow(STRING_formated))

STRING_formated_rev <- data.frame(gene1 = STRING_formated$gene2, gene2 = STRING_formated$gene1, STRING_score = STRING_formated$STRING_score, double = rep("exclude", nrow(STRING_formated)))
STRING_formated <- rbind(STRING_formated, STRING_formated_rev)

STRING_formated$gene_mouse <- STRING_formated$gene1 %>% sapply(function(gene) {
    return(ifelse(gene == "SHH", "SHH", figure_genes[gene, "mouse"]))
})

STRING_formated$gene_SHH_link <- c(1:nrow(STRING_formated)) %>% sapply(function(i) {
    return(ifelse(STRING_formated$gene2[i] == "SHH", TRUE, FALSE))
})

STRING_formated$gene_SHH_corr <- STRING_formated$gene1 %>% sapply(function(gene) {
    return(WGCNA_genes[gene, "cor"])
})

write.csv(
    filter(
        STRING_formated,
        ((gene1 %in% filter(WGCNA_genes, cor > 0)$gene & gene2 %in% filter(WGCNA_genes, cor > 0)$gene) & (gene1 %in% filter(figure_genes, sc_w3_w4 == "expressed")$gene & gene2 %in% filter(figure_genes, sc_w3_w4 == "expressed")$gene))
    ),
    file = "results/tables/Figure_4/STRING_formated_ventral.csv", row.names = FALSE, quote = FALSE
)

write.csv(
    filter(
        STRING_formated,
        (((gene1 %in% filter(WGCNA_genes, cor < 0)$gene & gene2 %in% filter(WGCNA_genes, cor < 0)$gene) | (gene1 == "SHH" & gene2 %in% filter(WGCNA_genes, cor < 0)$gene) | (gene1 %in% filter(WGCNA_genes, cor < 0)$gene & gene2 == "SHH")) & (gene1 %in% filter(figure_genes, sc_w3_w4 == "expressed")$gene & gene2 %in% filter(figure_genes, sc_w3_w4 == "expressed")$gene))
    ),
    file = "results/tables/Figure_4/STRING_formated_dorsal.csv", row.names = FALSE, quote = FALSE
)