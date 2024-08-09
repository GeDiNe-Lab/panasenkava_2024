# Loading packages and functions
library(Matrix)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(VennDiagram)
library(tibble)

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
rawcounts <- readcounts("/home/jules/Documents/phd/Data/lab_RNAseq/manip4/manip4_counts.csv")
meta <- read.table("/home/jules/Documents/phd/Data/lab_RNAseq/manip4/manip4_metadata.csv", sep = ",", header = T)

# L9C1_2 is an outlier and is removed
meta <- filter(meta, type %in% c("cyclo", "ventral") & samples != "L9C1_2")
counts <- rawcounts[which(rowSums(rawcounts) >= 50), meta$samples]

meta$cyclo_dose_qual <- meta$cyclo_dose %>% sapply(function(x) {
    if (x %in% c(0.125, 0.25)) {
        return("low")
    } else if (x %in% c(0.5, 1)) {
        return("high")
    } else {
        return("no_cyclo")
    }
})

# making DESeq object
dds <- DESeqDataSetFromMatrix(
    countData = counts[, meta$samples],
    colData = meta,
    design = ~cyclo_dose_qual
)

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

pca.data <- plotPCA.DESeqTransform(vsd, intgroup = c("type", "cyclo_dose_qual"), returnData = TRUE)
percentVar <- round(100 * attr(pca.data, "percentVar"))

png(filename = "results/images/Figure_2A/F3_PCA_1_2.png", width = 1600, height = 1200, res = 250)
ggplot(pca.data, aes(PC1, PC2, color = type, shape = cyclo_dose_qual)) +
    geom_point(size = 2, stroke = 1) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_color_manual(values = c("#ecb039", "#80AD3C")) +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6)) +
    custom_theme() +
    ggtitle("First and second PCs of dorsal and ventral kinetics all genes")
dev.off()

meta_bin <- meta %>%
    dplyr::select(c("type", "cyclo_dose")) %>%
    apply(2, function(x) {
        return(as.numeric(factor(x)) - 1)
    }) %>%
    as.matrix()

PC_covariate_cor <- cor(pca.data[, 1:5], meta_bin) %>% abs()
rownames(PC_covariate_cor) <- paste0(rownames(PC_covariate_cor), " (", percentVar[1:5], "%)")
PC_covariate_cor

png(filename = "results/images/Figure_2A/F3_PC_covariate_correlation.png", width = 2000, height = 1800, res = 250)
Heatmap(
    PC_covariate_cor,
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.2f", PC_covariate_cor[i, j]), x, y, gp = gpar(fontsize = 10, fontface = "bold", col = "#646464"))
    },
    name = "Absolute pearson correlation",
    row_title_gp = gpar(fontsize = 20, fontface = "bold"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    column_names_side = "top",
    column_names_rot = 0,
    column_names_centered = TRUE,
    show_column_names = TRUE,
    show_heatmap_legend = TRUE,
    width = ncol(PC_covariate_cor) * unit(2.5, "cm"),
    height = nrow(PC_covariate_cor) * unit(1.5, "cm"),
    col = colorRampPalette(c(
        "lightblue",
        "darkblue"
    ))(1000),
)
dev.off()


# getting genes correlated to cyclopamine dose
cyclo_genes <- cor(t(assay(vsd)[, meta$samples[order(meta$cyclo_dose)]]), sort(meta$cyclo_dose)) %>% as.data.frame()
colnames(cyclo_genes) <- c("cor")
cyclo_genes$genes <- rownames(cyclo_genes) %>% gene_converter("ENSEMBL", "SYMBOL")
cyclo_genes$abscor <- abs(cyclo_genes$cor)
cyclo_genes_f <- filter(cyclo_genes, abscor > 0.5)

# DEGs high vs no cyclo, high vs low cyclo, low vs no cyclo
DE_high_vs_no <- dds %>%
    DESeq() %>%
    results(alpha = 0.01, contrast = c("cyclo_dose_qual", "high", "no_cyclo")) %>%
    as.data.frame()
DE_high_vs_no$gene <- gene_converter(rownames(DE_high_vs_no), "ENSEMBL", "SYMBOL")
DE_high_vs_no_f <- filter(DE_high_vs_no, padj < 0.05 & abs(log2FoldChange) >= 2 & !is.na(gene))

DE_high_vs_low <- dds %>%
    DESeq() %>%
    results(alpha = 0.01, contrast = c("cyclo_dose_qual", "high", "low")) %>%
    as.data.frame()
DE_high_vs_low$gene <- gene_converter(rownames(DE_high_vs_low), "ENSEMBL", "SYMBOL")
DE_high_vs_low_f <- filter(DE_high_vs_low, padj < 0.05 & abs(log2FoldChange) >= 2 & !is.na(gene))

DE_low_vs_no <- dds %>%
    DESeq() %>%
    results(alpha = 0.01, contrast = c("cyclo_dose_qual", "low", "no_cyclo")) %>%
    as.data.frame()
DE_low_vs_no$gene <- gene_converter(rownames(DE_low_vs_no), "ENSEMBL", "SYMBOL")
DE_low_vs_no_f <- filter(DE_low_vs_no, padj < 0.05 & abs(log2FoldChange) >= 2 & !is.na(gene))

# Preparation for heatmap, clustering and GO enrichment
sample_order <- meta$samples[order(meta$cyclo_dose)]
scaled_mat <- t(apply(assay(vsd)[rownames(cyclo_genes_f), sample_order], 1, scale))
colnames(scaled_mat) <- colnames(assay(vsd)[, sample_order])

clustering <- hclust(dist(scaled_mat))
clusters <- cutree(clustering, k = 2)

sub_clusters_list <- unique(clusters) %>% lapply(function(cluster) {
    sub_mat <- scaled_mat[names(clusters[which(clusters == cluster)]), sample_order]
    sub_clustering <- hclust(dist(sub_mat))
    return(cutree(sub_clustering, k = 5))
})
names(sub_clusters_list) <- paste0("cluster_", unique(clusters))
sub_clusters <- sub_clusters_list %>%
    unname() %>%
    unlist()
sub_clusters <- sub_clusters[names(clusters)]

clusters_ha <- rowAnnotation(
    cluster = as.character(clusters[clustering$order]),
    sub_cluster = as.character(sub_clusters[clustering$order]),
    col = list(
        cluster = c(
            "1" = "red",
            "2" = "blue"
        ),
        sub_cluster = c(
            "1" = "black",
            "2" = "pink",
            "3" = "yellow",
            "4" = "brown",
            "5" = "grey"
        )
    )
)

for (cluster in unique(clusters)) {
    print(cluster)
    GO_enrichment <- clusterProfiler::enrichGO(names(clusters[which(clusters == cluster)]),
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


png(filename = "results/images/Figure_3/F3_cyclo_genes_HM.png", width = 2400, height = 1600, res = 250)
Heatmap(
    scaled_mat[clustering$order, sample_order],
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

# Building heatmap genes annotation (DE, correlation with cyclo, cluster, subcluster)
colnames(DE_high_vs_no) <- paste0("HvsN_", colnames(DE_high_vs_no))
colnames(DE_high_vs_low) <- paste0("HvsL_", colnames(DE_high_vs_low))
colnames(DE_low_vs_no) <- paste0("LvsN_", colnames(DE_low_vs_no))
cyclo_genes_df <- Reduce(cbind, list(
    dplyr::select(DE_high_vs_no, c("HvsN_padj", "HvsN_log2FoldChange")),
    dplyr::select(DE_high_vs_low, c("HvsL_padj", "HvsL_log2FoldChange")),
    dplyr::select(DE_low_vs_no, c("LvsN_padj", "LvsN_log2FoldChange"))
))
cyclo_genes_df$HvsN_thr <- ifelse(abs(cyclo_genes_df$HvsN_log2FoldChange) >= 2 & cyclo_genes_df$HvsN_padj < 0.01, "DE", "no")
cyclo_genes_df$HvsL_thr <- ifelse(abs(cyclo_genes_df$HvsL_log2FoldChange) >= 2 & cyclo_genes_df$HvsL_padj < 0.01, "DE", "no")
cyclo_genes_df$LvsN_thr <- ifelse(abs(cyclo_genes_df$LvsN_log2FoldChange) >= 2 & cyclo_genes_df$LvsN_padj < 0.01, "DE", "no")

cyclo_genes_df$cyclo_cor <- ifelse(rownames(cyclo_genes_df) %in% rownames(cyclo_genes_f), "yes", "no")
cyclo_genes_df$cluster <- rownames(cyclo_genes_df) %>% sapply(function(gene) {
    if (gene %in% names(clusters)) {
        return(clusters[gene])
    } else {
        return("no")
    }
})
cyclo_genes_df$sub_cluster <- rownames(cyclo_genes_df) %>% sapply(function(gene) {
    if (gene %in% names(sub_clusters)) {
        return(sub_clusters[gene])
    } else {
        return("no")
    }
})
cyclo_genes_df$genes <- gene_converter(rownames(cyclo_genes_df), "ENSEMBL", "SYMBOL")
View(cyclo_genes_df)

write.csv(cyclo_genes_df, "results/tables/Figure_3/cyclo_genes_df.csv")

# Genes tested in mices
ventral <- c(
    "AFF2", "ATP2C2", "AUTS2", "BAHCC1", "CAPN6", "CNTN6", "EDNRA", "FOXA1", "FOXA2", "FRZB", "GPM6B", "GRIK3", "HTR1D", "LDB2", "LINC00261", "MBIP", "MPPED1", "MYRF", "NAALAD2", "NACC2", "NKX2-1", "NKX2-1-AS1", "NKX2-2", "NTN1", "PDZRN3", "PLCL1", "PNMA2", "PNRC2", "PPM1L", "PTCH1", "QKI", "RGMA", "RORA", "RPS6KA6", "RXRA", "SERPINF1", "SERPINI1", "SFRP1", "SHH", "SLC38A2", "SLC38A4", "SLIT1", "SMIM32", "SPON1", "SPTSSB", "TMTC2", "TRIM9", "USP2"
)
dorsal <- c("PAX6", "ADD3", "ATP2B1", "CNN3", "COLGALT2", "EPHA4", "FZD3", "GLI2", "GLI3", "HOMER1", "NLGN1", "NUAK2", "OPTN", "PALLD", "PDP1", "PLK2", "PRDX6", "SLC3A2", "VCL", "ZIC2", "ZIC5", "ZNF385B", "ZFHX4")

# Known SHH related genes
known_genes <- c("GLI2", "GLI3", "ZIC2", "FOXA1", "FOXA2", "NKX2-1", "PAX6", "PTCH1")

# co-expression

# getting normalized counts for co-expression
counts_coex <- rawcounts[, meta$samples]
counts_coex <- counts_coex[which(rowSums(counts_coex) >= 200), ]
dds_coex <- DESeqDataSetFromMatrix(
    countData = counts_coex,
    colData = meta,
    design = ~cyclo_dose_qual
)
vsd_coex <- vst(dds_coex, blind = FALSE)

# getting co-expression matrix (Pearson correlation)
corr <- WGCNA::cor(t(assay(vsd_coex)))

# getting genes correlation with SHH, NKX2-1 and PAX6
corr_df <- as.data.frame(corr[, c("ENSG00000164690", "ENSG00000136352", "ENSG00000007372")])
corr_df$abs_max_cor <- apply(corr_df, 1, function(x) max(abs(x)))
corr_df$gene <- gene_converter(rownames(corr_df), "ENSEMBL", "SYMBOL")
corr_df <- filter(corr_df, !is.na(gene))

# Getting genes passing 2 correlation thresholds
corr_df_f1 <- filter(corr_df, abs_max_cor >= 0.8)
corr_df_f2 <- filter(corr_df, abs_max_cor >= 0.85)

SHH_pos_1 <- data.frame(
    cor = filter(corr_df_f1, ENSG00000164690 > 0)$abs_max_cor,
    link = rep(1, nrow(filter(corr_df_f1, ENSG00000164690 > 0))),
    target = filter(corr_df_f1, ENSG00000164690 > 0)$gene,
    source = rep("SHH", nrow(filter(corr_df_f1, ENSG00000164690 > 0)))
)
SHH_pos_1$target_info <- SHH_pos_1$target %>% sapply(function(x) {
    if (x %in% known_genes) {
        return("known")
    } else if (x %in% union(dorsal, ventral)) {
        return("selected")
    } else {
        return("no")
    }
})
SHH_pos_1$ttarget_info <- SHH_pos_1$target_info

SHH_neg_1 <- data.frame(
    cor = filter(corr_df_f1, ENSG00000164690 < 0)$abs_max_cor,
    link = rep(1, nrow(filter(corr_df_f1, ENSG00000164690 < 0))),
    target = filter(corr_df_f1, ENSG00000164690 < 0)$gene,
    source = rep("SHH", nrow(filter(corr_df_f1, ENSG00000164690 < 0)))
)
SHH_neg_1$target_info <- SHH_neg_1$target %>% sapply(function(x) {
    if (x %in% known_genes) {
        return("known")
    } else if (x %in% union(dorsal, ventral)) {
        return("selected")
    } else {
        return("no")
    }
})
SHH_neg_1$ttarget_info <- SHH_neg_1$target_info

SHH_pos_2 <- data.frame(
    cor = filter(corr_df_f2, ENSG00000164690 > 0)$abs_max_cor,
    link = rep(1, nrow(filter(corr_df_f2, ENSG00000164690 > 0))),
    target = filter(corr_df_f2, ENSG00000164690 > 0)$gene,
    source = rep("SHH", nrow(filter(corr_df_f2, ENSG00000164690 > 0)))
)
SHH_pos_2$target_info <- SHH_pos_2$target %>% sapply(function(x) {
    if (x %in% known_genes) {
        return("known")
    } else if (x %in% union(dorsal, ventral)) {
        return("selected")
    } else {
        return("no")
    }
})
SHH_pos_2$ttarget_info <- SHH_pos_2$target_info

SHH_neg_2 <- data.frame(
    cor = filter(corr_df_f2, ENSG00000164690 < 0)$abs_max_cor,
    link = rep(1, nrow(filter(corr_df_f2, ENSG00000164690 < 0))),
    target = filter(corr_df_f2, ENSG00000164690 < 0)$gene,
    source = rep("SHH", nrow(filter(corr_df_f2, ENSG00000164690 < 0)))
)
SHH_neg_2$target_info <- SHH_neg_2$target %>% sapply(function(x) {
    if (x %in% known_genes) {
        return("known")
    } else if (x %in% union(dorsal, ventral)) {
        return("selected")
    } else {
        return("no")
    }
})
SHH_neg_2$ttarget_info <- SHH_neg_2$target_info

write.csv(SHH_pos_1, file = "results/tables/Figure_3/cytoscape_SHH_pos_0_80.csv", row.names = FALSE)
write.csv(SHH_neg_1, file = "results/tables/Figure_3/cytoscape_SHH_neg_0_80.csv", row.names = FALSE)
write.csv(SHH_pos_2, file = "results/tables/Figure_3/cytoscape_SHH_pos_0_85.csv", row.names = FALSE)
write.csv(SHH_neg_2, file = "results/tables/Figure_3/cytoscape_SHH_neg_0_85.csv", row.names = FALSE)

# Making volcano plots
DE_cyclo <- list(
    high_vs_no = DE_high_vs_no_f,
    high_vs_low = DE_high_vs_low_f,
    low_vs_no = DE_low_vs_no_f
)

for (contrast in names(DE_cyclo)) {
    print(contrast)
    DE <- DE_cyclo[[contrast]]
    ggplot(DE, aes(x = log2FoldChange, y = -log10(padj), label = gene)) +
        ggrepel::geom_text_repel(box.padding = 0.001, size = 2.5, max.overlaps = 20) +
        custom_theme() +
        geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
        labs(x = "log2FoldChange", y = "-log10(padj)", title = paste0("DE: ", contrast, " ", nrow(DE), " DE genes total"))
    ggsave(filename = paste0("results/images/Figure_3/Volcano_", contrast, ".png"), units = "px", width = 1800, height = 1400, dpi = 250)
}
