# Loading packages and functions
library(Matrix)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(VennDiagram)
library(tibble)

rstudioapi::getSourceEditorContext()$path %>%
    str_split("/") %>%
    unlist() %>%
    head(-3) %>%
    str_c(collapse = "/") %>%
    str_c("/") %>%
    setwd()

source("R/custom_fct.R")
rawcounts <- readcounts("/home/jules/Documents/phd/Data/lab_RNAseq/manip4/manip4_counts.csv")
meta <- read.table("/home/jules/Documents/phd/Data/lab_RNAseq/manip4/manip4_metadata.csv", sep = ",", header = T)

meta <- filter(meta, type %in% c("cyclo", "ventral"))
counts <- rawcounts[which(rowSums(rawcounts) >= 50), meta$samples]
norm <- varianceStabilizingTransformation(counts)
meta %>% View()
meta$cyclo_dose_qual <- meta$cyclo_dose %>% sapply(function(x) {
    if (x %in% c(0.125, 0.25)) {
        return("low")
    } else if (x %in% c(0.5, 1)) {
        return("high")
    } else {
        return("no_cyclo")
    }
})

DE_no_high <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~cyclo_dose_qual
) %>%
    DESeq() %>%
    results(alpha = 0.01, contrast = c("cyclo_dose_qual", "no_cyclo", "high")) %>%
    as.data.frame()
DE_no_high$gene <- gene_converter(rownames(DE_no_high), "ENSEMBL", "SYMBOL")
DE_no_high <- filter(DE_no_high, padj < 0.05 & abs(log2FoldChange) > 1 & !is.na(gene))

DE_low_high <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~cyclo_dose_qual
) %>%
    DESeq() %>%
    results(alpha = 0.01, contrast = c("cyclo_dose_qual", "low", "high")) %>%
    as.data.frame()
DE_low_high$gene <- gene_converter(rownames(DE_low_high), "ENSEMBL", "SYMBOL")
DE_low_high <- filter(DE_low_high, padj < 0.05 & abs(log2FoldChange) > 1 & !is.na(gene))

DE_no_low <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~cyclo_dose_qual
) %>%
    DESeq() %>%
    results(alpha = 0.01, contrast = c("cyclo_dose_qual", "no_cyclo", "low")) %>%
    as.data.frame()
DE_no_low$gene <- gene_converter(rownames(DE_no_low), "ENSEMBL", "SYMBOL")
DE_no_low <- filter(DE_no_low, padj < 0.05 & abs(log2FoldChange) > 1 & !is.na(gene))



scaled_mat <- t(apply(norm[rownames(DE_no_high), ], 1, scale))
colnames(scaled_mat) <- colnames(norm)

clustering <- hclust(dist(scaled_mat))

clusters_1 <- cutree(clustering, k = 2)
clusters_2 <- cutree(clustering, k = 5)
clusters_3 <- cutree(clustering, k = 7)

clusters_ha <- rowAnnotation(
    cluster_1 = as.character(clusters_1[clustering$order]),
    cluster_2 = as.character(clusters_2[clustering$order]),
    cluster_3 = as.character(clusters_3[clustering$order]),
    col = list(
        cluster_1 = c(
            "1" = "red",
            "2" = "blue"
        ),
        cluster_2 = c(
            "1" = "red",
            "2" = "blue",
            "3" = "green",
            "4" = "purple",
            "5" = "orange"
        ),
        cluster_3 = c(
            "1" = "red",
            "2" = "blue",
            "3" = "green",
            "4" = "purple",
            "5" = "orange",
            "6" = "black",
            "7" = "pink"
        )
    )
)



for (cluster in unique(clusters_1)) {
    print(cluster)
    GO_enrichment <- clusterProfiler::enrichGO(names(clusters_1[which(clusters_1 == cluster)]),
        OrgDb = "org.Hs.eg.db",
        keyType = "ENSEMBL",
        ont = "BP"
    )
    write.csv(GO_enrichment, paste0("results/tables/Figure_3/GO_enrichment_cluster_", cluster, ".csv"))
    goplot <- clusterProfiler::dotplot(GO_enrichment,
        title = paste0("GO enrichment on cluster", cluster, " (biological processes only)"),
        showCategory = 15
    )
    ggsave(paste0("results/images/Figure_3/GO_enrichment_cluster_", cluster, ".png"), goplot, width = 8, height = 10)
}


png(filename = "results/images/Figure_3/F3_cyclo_DE_HM.png", width = 2400, height = 1600, res = 250)
Heatmap(
    scaled_mat[clustering$order, order(meta$cyclo_dose)],
    name = "Normalized expression",
    column_names_gp = gpar(fontsize = 6),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    left_annotation = clusters_ha,
    show_row_names = FALSE,
    row_names_side = "left",
    show_column_names = TRUE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    width = ncol(scaled_mat) * unit(4, "mm"),
    # height = nrow(mat) * unit(5, "mm"),
    col = colorRampPalette(c(
        "black",
        "purple",
        "orange",
        "yellow"
    ))(1000),
)
dev.off()

DE_no_high$cluster_1 <- clusters_1[rownames(DE_no_high)]
DE_no_high$cluster_2 <- clusters_2[rownames(DE_no_high)]
DE_no_high$cluster_3 <- clusters_3[rownames(DE_no_high)]

write.csv(DE_no_high, "results/tables/Figure_3/DEG_cyclo_vAN_VS_high.csv")
write.csv(DE_low_high, "results/tables/Figure_3/DEG_cyclo_low_VS_high.csv")
write.csv(DE_no_low, "results/tables/Figure_3/DEG_cyclo_vAN_VS_low.csv")

ventral <- c(
    "AFF2", "ATP2C2", "AUTS2", "BAHCC1", "CAPN6", "CNTN6", "EDNRA", "FOXA1", "FOXA2", "FRZB", "GPM6B", "GRIK3", "HTR1D", "LDB2", "LINC00261", "MBIP", "MPPED1", "MYRF", "NAALAD2", "NACC2", "NKX2-1", "NKX2-1-AS1", "NKX2-2", "NTN1", "PDZRN3", "PLCL1", "PNMA2", "PNRC2", "PPM1L", "PTCH1", "QKI", "RGMA", "RORA", "RPS6KA6", "RXRA", "SERPINF1", "SERPINI1", "SFRP1", "SHH", "SLC38A2", "SLC38A4", "SLIT1", "SMIM32", "SPON1", "SPTSSB", "TMTC2", "TRIM9", "USP2"
)
dorsal <- c("PAX6", "ADD3", "ATP2B1", "CNN3", "COLGALT2", "EPHA4", "FZD3", "GLI2", "GLI3", "HOMER1", "NLGN1", "NUAK2", "OPTN", "PALLD", "PDP1", "PLK2", "PRDX6", "SLC3A2", "VCL", "ZIC2", "ZIC5", "ZNF385B", "ZFHX4")
known_genes <- c("GLI2", "GLI3", "ZIC2", "FOXA1", "FOXA2", "NKX2-1", "PAX6", "PTCH1")
# co-expression :
counts_coex <- rawcounts[, meta$samples]
counts_coex <- counts_coex[which(rowSums(counts_coex) >= 100), ]
norm_coex <- varianceStabilizingTransformation(counts_coex)
corr <- WGCNA::cor(t(norm_coex))

corr_df <- as.data.frame(corr[, c("ENSG00000164690", "ENSG00000136352", "ENSG00000007372")])
corr_df$abs_max_cor <- apply(corr_df, 1, function(x) max(abs(x)))
corr_df$gene <- gene_converter(rownames(corr_df), "ENSEMBL", "SYMBOL")
corr_df <- filter(corr_df, !is.na(gene))

genelist <- filter(corr_df, abs_max_cor >= 0.8)$gene
lit_gene_list <- intersect(genelist, union(ventral, dorsal))
corr_df_f <- filter(corr_df, abs_max_cor >= 0.8)



cytoscape <- corr[rownames(filter(corr_df, gene %in% lit_gene_list)), rownames(filter(corr_df, gene %in% lit_gene_list))]
cytoscape[lower.tri(cytoscape, diag = TRUE)] <- NA
cytoscape <- data.table::melt(cytoscape)

cytoscape <- filter(cytoscape, is.na(value) == FALSE)

cytoscape$abs_cor <- cytoscape$value %>% abs()
cytoscape$Var1 <- cytoscape$Var1 %>%
    as.vector() %>%
    gene_converter("ENSEMBL", "SYMBOL")
cytoscape$Var2 <- cytoscape$Var2 %>%
    as.vector() %>%
    gene_converter("ENSEMBL", "SYMBOL")
cytoscape$link <- ifelse(cytoscape$value > 0, "pos", "neg")
cytoscape$knownV1 <- c(1:nrow(cytoscape)) %>% sapply(function(x) {
    if (cytoscape$Var1[x] %in% known_genes) {
        return("known")
    } else if (cytoscape$Var1[x] == "SHH") {
        return("SHH")
    } else {
        return("no")
    }
})
cytoscape$knownV2 <- c(1:nrow(cytoscape)) %>% sapply(function(x) {
    if (cytoscape$Var2[x] %in% known_genes) {
        return("known")
    } else if (cytoscape$Var2[x] == "SHH") {
        return("SHH")
    } else {
        return("no")
    }
})

# Add a column for node degree
cytoscape$degree <- sapply(cytoscape$Var1, function(node) sum(cytoscape$Var1 == node | cytoscape_all$Var2 == node))
View(filter(cytoscape, abs_cor >= 0.9))
write.csv(filter(cytoscape, abs_cor >= 0.9), "results/tables/Figure_3/cytoscape_all.csv")
write.csv(genelist, "results/tables/Figure_3/coex_genelist.csv")

SHH_pos <- filter(cytoscape, value >= 0.8 & (Var1 == "ENSG00000164690" | Var2 == "ENSG00000164690"))
SHH_pos$Var1 <- SHH_pos$Var1 %>%
    as.vector() %>%
    gene_converter("ENSEMBL", "SYMBOL")
SHH_pos$Var2 <- SHH_pos$Var2 %>%
    as.vector() %>%
    gene_converter("ENSEMBL", "SYMBOL")

SHH_neg <- filter(cytoscape, value <= -0.8 & (Var1 == "ENSG00000164690" | Var2 == "ENSG00000164690"))
SHH_neg$Var1 <- SHH_neg$Var1 %>%
    as.vector() %>%
    gene_converter("ENSEMBL", "SYMBOL")
SHH_neg$Var2 <- SHH_neg$Var2 %>%
    as.vector() %>%
    gene_converter("ENSEMBL", "SYMBOL")


write.csv(SHH_pos, "results/tables/Figure_3/SHH_coex_pos.csv")
write.csv(SHH_neg, "results/tables/Figure_3/SHH_coex_neg.csv")
