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

#  get WGCNA gene infos table and figures genes info table
WGCNA_genes <- read.csv("results/tables/Figure_4/SHH_cluster.csv", header = TRUE)
rownames(WGCNA_genes) <- WGCNA_genes$gene
figure_genes <- read.csv("results/tables/Figure_genelist.csv", header = TRUE)
#  Got a doublet symbol for CKMT1B
figure_genes$gene[which(figure_genes$X == "ENSG00000223572")] <- "CKMT1B_b"
figure_genes <- filter(figure_genes, !is.na(gene))
rownames(figure_genes) <- figure_genes$gene

#  Get STRINGdb interactions results table and format it
STRING_res <- read.table("results/tables/Figure_4/string_interactions_short.tsv", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
STRING_formated <- STRING_res %>% dplyr::select(X.node1, node2, combined_score)
colnames(STRING_formated) <- c("gene1", "gene2", "STRING_score")
STRING_formated$double <- rep("keep", nrow(STRING_formated))
STRING_formated_rev <- data.frame(gene1 = STRING_formated$gene2, gene2 = STRING_formated$gene1, STRING_score = STRING_formated$STRING_score, double = rep("exclude", nrow(STRING_formated)))
STRING_formated <- rbind(STRING_formated, STRING_formated_rev)

# Retaking mouse genes without any STRING_connection (they are linked to themselves in the network)
retake_genes <- setdiff(filter(figure_genes, mouse == "tested")$gene, c(STRING_formated$gene1, STRING_formated$gene2) %>% unique()) %>%
    lapply(function(gene) {
        return(data.frame(gene1 = gene, gene2 = gene, STRING_score = 0.4, double = "exclude"))
    }) %>%
    Reduce(rbind, .)
STRING_formated <- rbind(STRING_formated, retake_genes)

# Tagging mice tested genes
STRING_formated$gene_mouse <- STRING_formated$gene1 %>% sapply(function(gene) {
    return(ifelse(gene == "SHH", "SHH", figure_genes[gene, "mouse"]))
})

#  Tagging scRNA expressed genes
STRING_formated$scRNA <- STRING_formated$gene1 %>% sapply(function(gene) {
    return(figure_genes[gene, "sc_w3_w4"])
})

# Tagging genes directly linked to SHH
STRING_formated$gene_SHH_link <- c(1:nrow(STRING_formated)) %>% sapply(function(i) {
    return(ifelse(STRING_formated$gene2[i] == "SHH", TRUE, FALSE))
})

# Getting correlation to SHH info for mice tested genes
mouse_genes_cor <- read.csv("results/tables/Figure_4/mouse_genes_cor.csv", header = TRUE)
rownames(mouse_genes_cor) <- mouse_genes_cor$X %>% gene_converter("ENSEMBL", "SYMBOL")

#  Getting correlation to SHH info for gene1 column genes
STRING_formated$gene1_SHH_corr <- STRING_formated$gene1 %>% sapply(function(gene) {
    if (gene %in% rownames(mouse_genes_cor)) {
        return(mouse_genes_cor[gene, "cor"])
    } else {
        return(WGCNA_genes[gene, "cor"])
    }
})

#  Getting correlation to SHH info for gene2 column genes
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



fgl <- read.csv("results/tables/Figure_genelist.csv")

nd1 <- read.table("results/tables/string_node_degrees.tsv", header = FALSE)
nd2 <- read.table("results/tables/string_node_degrees(1).tsv", header = FALSE)
nd3 <- read.table("results/tables/string_node_degrees(2).tsv", header = FALSE)

nd <- rbind(nd1, nd2, nd3)

STRING_gene <- unique(nd$V1)

sn1 <- read.table("results/tables/string_interactions_short.tsv", header = FALSE)
sn2 <- read.table("results/tables/string_interactions_short(1).tsv", header = FALSE)
sn3 <- read.table("results/tables/string_interactions_short(2).tsv", header = FALSE)

sn <- rbind(sn1, sn2, sn3)

sn_SHH <- filter(sn, V1 == "SHH" | V2 == "SHH")


fgl$STRINGdb <- fgl$genes %>% sapply(function(gene) {
    if (is.na(gene)) {
        return("non_coding")
    } else if (gene == "SHH") {
        return("SHH")
    } else if (gene %in% sn$V1 | gene %in% sn$V2) {
        return("coding_SHH_linked")
    } else if (gene %in% STRING_gene) {
        return("coding")
    } else {
        return("non_coding")
    }
})

write.csv(fgl, file = "results/tables/Figure_genelist_complete.csv", row.names = FALSE, quote = FALSE)
