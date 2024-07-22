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
meta <- filter(meta, sample != "LON71_D12_2", type %in% c("ventral", "dorsal"), line %in% c("LON71", "WTC"))
counts <- rawcounts[which(rowSums(rawcounts) >= 50), meta$sample]

dds <- DESeqDataSetFromMatrix(
    countData = counts[, meta$sample],
    colData = meta,
    design = ~ type + day
)

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# PCA plot
pca.data <- plotPCA(vsd, intgroup = c("type", "day", "line"), returnData = TRUE)
percentVar <- round(100 * attr(pca.data, "percentVar"))

png(filename = "results/images/Figure_2A/F2A_1_PCA_type.png", width = 1600, height = 1200, res = 250)
ggplot(pca.data, aes(PC1, PC2, color = type, shape = day)) +
    geom_point(size = 2, stroke = 1) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_color_manual(values = c("#A1A1DE", "#80AD3C")) +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6)) +
    custom_theme() +
    ggtitle("PCA of dorsal and ventral kinetics")
dev.off()

png(filename = "results/images/Figure_2A/F2A_1_PCA_line.png", width = 1600, height = 1200, res = 250)
ggplot(pca.data, aes(PC1, PC2, color = line, shape = day)) +
    geom_point(size = 2, stroke = 1) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6)) +
    custom_theme() +
    ggtitle("PCA of dorsal and ventral kinetics")
dev.off()

DEGs_day_ventral <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "ventral" & day %in% c("day02", "day12"))$sample],
    colData = filter(meta, type == "ventral" & day %in% c("day02", "day12")),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day12", "day02")) %>%
    as.data.frame() %>%
    na.omit()

DEGs_day_ventral_f <- filter(DEGs_day_ventral, padj < 0.01, abs(log2FoldChange) >= 1)
DEGs_day_ventral_f$gene <- gene_converter(rownames(DEGs_day_ventral_f), "ENSEMBL", "SYMBOL")
DEGs_day_ventral_f <- filter(DEGs_day_ventral_f, !is.na(gene))
DEGs_day_ventral_f$type <- rep("ventral", nrow(DEGs_day_ventral_f))


DEGs_day_dorsal <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, type == "dorsal" & day %in% c("day02", "day12"))$sample],
    colData = filter(meta, type == "dorsal" & day %in% c("day02", "day12")),
    design = ~ line + day
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("day", "day12", "day02")) %>%
    as.data.frame() %>%
    na.omit()

DEGs_day_dorsal_f <- filter(DEGs_day_dorsal, padj < 0.01, abs(log2FoldChange) >= 1)
DEGs_day_dorsal_f$gene <- gene_converter(rownames(DEGs_day_dorsal_f), "ENSEMBL", "SYMBOL")
DEGs_day_dorsal_f <- filter(DEGs_day_dorsal_f, !is.na(gene))
DEGs_day_dorsal_f$type <- rep("dorsal", nrow(DEGs_day_dorsal_f))

degs_total <- union(rownames(DEGs_day_ventral_f), rownames(DEGs_day_dorsal_f))




vstnorm <- vst(dds, blind = FALSE)
mat <- assay(vstnorm)[degs_total, ]
# scaling the matrix
scaled_mat <- t(apply(mat, 1, scale))
colnames(scaled_mat) <- colnames(mat)

sample_order <- c(
    filter(meta, type == "ventral")$sample[order(filter(meta, type == "ventral")$day, decreasing = TRUE)],
    filter(meta, type == "dorsal")$sample[order(filter(meta, type == "dorsal")$day)]
)

clustering <- hclust(dist(scaled_mat[, sample_order]))

clusters_1 <- cutree(clustering, k = 5)
clusters_2 <- cutree(clustering, k = 7)
clusters_3 <- cutree(clustering, k = 10)

clusters_ha <- rowAnnotation(
    cluster_1 = as.character(clusters_1[clustering$order]),
    cluster_2 = as.character(clusters_2[clustering$order]),
    cluster_3 = as.character(clusters_3[clustering$order]),
    col = list(
        cluster_1 = c(
            "1" = "red",
            "2" = "blue",
            "3" = "green",
            "4" = "purple",
            "5" = "orange"
        ),
        cluster_2 = c(
            "1" = "red",
            "2" = "blue",
            "3" = "green",
            "4" = "purple",
            "5" = "orange",
            "6" = "black",
            "7" = "pink"
        ),
        cluster_3 = c(
            "1" = "red",
            "2" = "blue",
            "3" = "green",
            "4" = "purple",
            "5" = "orange",
            "6" = "black",
            "7" = "pink",
            "8" = "yellow",
            "9" = "brown",
            "10" = "grey"
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
    write.csv(GO_enrichment, paste0("results/tables/Figure_2A/GO_enrichment_cluster_", cluster, ".csv"))
    goplot <- clusterProfiler::dotplot(GO_enrichment,
        title = paste0("GO enrichment on cluster", cluster, " (biological processes only)"),
        showCategory = 15
    )
    ggsave(paste0("results/images/Figure_2A/F2A_DE_GO_clust", cluster, ".png"), goplot, width = 8, height = 10)
}

png(filename = "results/images/Figure_2A/F2A_DE_HM.png", width = 2400, height = 1600, res = 250)
Heatmap(
    scaled_mat[clustering$order, sample_order],
    name = "Normalized expression",
    column_names_gp = gpar(fontsize = 6),
    cluster_rows = FALSE,
    left_annotation = clusters_ha,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    row_names_side = "left",
    show_column_names = TRUE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    width = ncol(scaled_mat) * unit(2, "mm"),
    # height = nrow(mat) * unit(5, "mm"),
    col = colorRampPalette(c(
        "black",
        "purple",
        "orange",
        "yellow"
    ))(1000),
)
dev.off()

sample_order_LON <- c(
    filter(meta, type == "ventral" & line == "LON71")$sample[order(filter(meta, type == "ventral" & line == "LON71")$day, decreasing = TRUE)],
    filter(meta, type == "dorsal" & line == "LON71")$sample[order(filter(meta, type == "dorsal" & line == "LON71")$day)]
)
png(filename = "results/images/Figure_2A/F2A_DE_HM_LON.png", width = 2400, height = 1600, res = 250)
Heatmap(
    scaled_mat[clustering$order, sample_order_LON],
    name = "Normalized expression",
    column_names_gp = gpar(fontsize = 6),
    cluster_rows = FALSE,
    left_annotation = clusters_ha,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    row_names_side = "left",
    show_column_names = TRUE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    width = ncol(scaled_mat[, sample_order_LON]) * unit(2, "mm"),
    # height = nrow(mat) * unit(5, "mm"),
    col = colorRampPalette(c(
        "black",
        "purple",
        "orange",
        "yellow"
    ))(1000),
)
dev.off()


sample_order_WTC <- c(
    filter(meta, type == "ventral" & line == "WTC")$sample[order(filter(meta, type == "ventral" & line == "WTC")$day, decreasing = TRUE)],
    filter(meta, type == "dorsal" & line == "WTC")$sample[order(filter(meta, type == "dorsal" & line == "WTC")$day)]
)
png(filename = "results/images/Figure_2A/F2A_DE_HM_WTC.png", width = 2400, height = 1600, res = 250)
Heatmap(
    scaled_mat[clustering$order, sample_order_WTC],
    name = "Normalized expression",
    column_names_gp = gpar(fontsize = 6),
    cluster_rows = FALSE,
    left_annotation = clusters_ha,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    row_names_side = "left",
    show_column_names = TRUE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    width = ncol(scaled_mat[, sample_order_WTC]) * unit(2, "mm"),
    # height = nrow(mat) * unit(5, "mm"),
    col = colorRampPalette(c(
        "black",
        "purple",
        "orange",
        "yellow"
    ))(1000),
)
dev.off()


DE_day_dorso_ventral <- rbind(DEGs_day_dorsal_f, DEGs_day_ventral_f)
DE_day_dorso_ventral$cluster_1 <- clusters_1[rownames(DE_day_dorso_ventral)]
DE_day_dorso_ventral$cluster_2 <- clusters_2[rownames(DE_day_dorso_ventral)]
DE_day_dorso_ventral$cluster_3 <- clusters_3[rownames(DE_day_dorso_ventral)]
DE_day_dorso_ventral %>% View()
write.csv(DE_day_dorso_ventral, "results/tables/DEG_kinetic_day02VSday12_dorsal_and_ventral.csv")

load("/home/jules/Documents/phd/Data/literature/CORTECON/Cortecon_Barebones.RData")
cortecon_counts <- RNA.raw %>% as.matrix()
rm(ex.raw, RNA.raw, h19e, h19eg, lincRNA, RNAs)

rownames(cortecon_counts) <- gene_converter(rownames(cortecon_counts), "ENTREZID", "ENSEMBL")
cortecon_counts <- cortecon_counts[!duplicated(rownames(cortecon_counts)), ]
cortecon_counts <- cortecon_counts[which(!is.na(rownames(cortecon_counts))), ]
colnames(cortecon_counts) <- gsub(" ", "", colnames(cortecon_counts))
colnames(cortecon_counts) <- gsub("\\.", "_", colnames(cortecon_counts))
cortecton_meta <- read.table("/home/jules/Documents/phd/Data/literature/CORTECON/cortecon_meta.csv", sep = ",", header = TRUE)


comm_genes <- intersect(rownames(cortecon_counts), rownames(rawcounts))

merged_counts <- cbind(
    rawcounts[comm_genes, ],
    cortecon_counts[comm_genes, ]
)
merged_meta <- rbind(
    meta,
    cortecton_meta
)
merged_meta$dataset <- ifelse(merged_meta$sample %in% colnames(cortecon_counts), "cortecon", "lab")
merged_meta <- filter(merged_meta, type == "dorsal" & day != "day00")
# merged_meta <- filter(merged_meta, dataset == "cortecon" | (dataset == "lab" & day == "day12"))
merged_meta$grouped_day <- merged_meta$day %>% sapply(function(day) {
    if (day %in% c("day02", "day04", "day06")) {
        return("day02-06")
    } else if (day %in% c("day07", "day08", "day10")) {
        return("day7-10")
    } else if (day %in% c("day26", "day33")) {
        return("day26-33")
    } else if (day %in% c("day49-77")) {
        return("day49-77")
    } else {
        return(day)
    }
})

merged_meta %>% View()
merged_counts <- merged_counts[, merged_meta$sample]
merged_counts <- merged_counts[which(rowSums(merged_counts) >= 50), ]

merged_norm <- varianceStabilizingTransformation(merged_counts)

lab_norm <- varianceStabilizingTransformation(rawcounts[which(rowSums(rawcounts) >= 50), meta$sample])
corte_norm <- varianceStabilizingTransformation(cortecon_counts[which(rowSums(cortecon_counts) >= 50), ])
comm_genes_late <- intersect(rownames(lab_norm), rownames(corte_norm))
merged_norm_late <- cbind(
    lab_norm[comm_genes_late, ],
    corte_norm[comm_genes_late, ]
)
merged_norm_late <- merged_norm_late[, merged_meta$sample]
merged_norm_late <- limma::removeBatchEffect(merged_norm, merged_meta$dataset)
merged_norm <- limma::removeBatchEffect(merged_norm, merged_meta$dataset)


test_kmean <- kmeans(t(merged_norm), centers = 8)
test_kmean$cluster %>%
    as.data.frame() %>%
    View()
pca_res <- ggPCA(
    t(merged_norm),
    ncp = 5,
    # ind.sup = which(colnames(merged_norm_late) %in% filter(merged_meta, dataset == "lab")$sample),
    scale.unit = TRUE,
    graph = FALSE
)

meta_pca <- cbind(pca_res$gg.ind, merged_meta)

png(filename = "results/images/Figure_2A/F2A_cortecon.png", width = 2400, height = 1600, res = 250)
ggplot(data = meta_pca, aes(x = PC1, y = PC2, color = grouped_day, shape = dataset)) +
    geom_point(size = 3) +
    custom_theme()
dev.off()
