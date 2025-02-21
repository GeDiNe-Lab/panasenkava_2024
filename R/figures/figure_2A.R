# Loading packages and functions
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(ggrepel)
library(DESeq2)
library(paletteer)
library(DEGreport)
library(ggsignif)
library(patchwork)
library(grid)
library(Matrix)

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
meta <- filter(rawmeta, sample != "LON71_D12_2", diff == "diff13", line %in% c("LON71", "WTC"), ((manip == "veranika" & day != "day12") | (manip == "lauryane" & day == "day12")))

# filtering out lowly expressed genes
counts <- rawcounts[, meta$sample][which(rowSums(rawcounts[, meta$sample]) >= 25), ]

# making DESeq object with lineage,days and type as covariates
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ line + day + type
)

# Normalization by variance stabilizing transformation
vsd_blind <- vst(dds, blind = TRUE)
vsd <- vst(dds, blind = FALSE)

# PCA with 3000 genes
pca.data <- plotPCA.DESeqTransform(vsd_blind, intgroup = c("type", "day", "line"), returnData = TRUE, ntop = 3000)
percentVar <- round(100 * attr(pca.data, "percentVar"))

pca_plot <- ggplot(pca.data, aes(PC1, PC2, color = line, shape = day)) +
    geom_point(size = 3, stroke = 1.5) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_color_manual(values = c("#868686", "#000000")) +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6)) +
    custom_theme()
ggsave("results/images/Figure_2/Figure2B.png", pca_plot, width = 8, height = 6)

# Get PCA/covariates ANOVA results
PC_covariate_ANOVA <- pca_anova(
    pca_data = pca.data,
    metadata = meta,
    covariates = c("line", "type", "day")
)
#  Saving ANOVA results
write.csv(PC_covariate_ANOVA, "results/tables/Figure_2/Figure2_ANOVA.csv")

# Variance explained by each PCs
varplot <- ggplot(data.frame(perc = percentVar, PC = factor(colnames(pca.data[1:20]), levels = colnames(pca.data[1:20]))), aes(x = PC, y = perc)) +
    geom_bar(stat = "identity") +
    custom_theme(diag_text = TRUE) +
    ylim(0, 100) +
    ggtitle("Variation explained by each PCs with all genes")
ggsave("results/images/Figure_2/Figure2_percentVar.png", pca_plot, width = 8, height = 6)

######################################
######################################
# FIGURE 2 C : PC1-2 correlated genes heatmap

# PC3 is associated with lineage difference, we want genes higlhy correlated with PC1 and PC2 but not PC3
top_PC1 <- cor(pca.data$PC1, t(assay(vsd_blind))) %>% as.vector()
names(top_PC1) <- rownames(vsd_blind)
top_PC1 <- sort(abs(top_PC1), decreasing = TRUE)[1:1000] %>% names()

top_PC2 <- cor(pca.data$PC2, t(assay(vsd_blind))) %>% as.vector()
names(top_PC2) <- rownames(vsd_blind)
top_PC2 <- sort(abs(top_PC2), decreasing = TRUE)[1:1000] %>% names()

top_PC3 <- cor(pca.data$PC3, t(assay(vsd_blind))) %>% as.vector()
names(top_PC3) <- rownames(vsd_blind)
top_PC3 <- sort(abs(top_PC3), decreasing = TRUE)[1:1000] %>% names()

# get Heatmap genes
hm_genes <- setdiff(union(top_PC1, top_PC2), top_PC3)

# making matrix for heatmap, clustering and GO enrichment
vsd_hm <- assay(vsd)[hm_genes, ]

scaled_mat <- t(apply(vsd_hm, 1, scale))
colnames(scaled_mat) <- colnames(vsd_hm)

#  Making heatmap for WTC lineage only
sample_order_WTC <- c(
    filter(meta, type == "ventral" & line == "WTC")$sample[order(filter(meta, type == "ventral" & line == "WTC")$day, decreasing = TRUE)],
    filter(meta, type == "dorsal" & line == "WTC")$sample[order(filter(meta, type == "dorsal" & line == "WTC")$day)]
)

clustering_WTC <- hclust(dist(scaled_mat[, sample_order_WTC]))
clusters_WTC <- cutree(clustering_WTC, k = 4)

# Check out cluster sizes and order
row_split <- factor(
    clusters_WTC[clustering_WTC$order],
    levels = c(1, 2, 3, 4),
    labels = c(
        "Early differentiation\nn=585",
        "Forebrain\nn=428",
        "Ventral forebrain\nn=524",
        "Dorsal forebrain\nn=463"
    )
)

# Getting gene list for Figure 2C
fig_genelist_2 <- row_split %>% as.data.frame()
colnames(fig_genelist_2) <- "Figure2C"
fig_genelist_2$Figure2C <- fig_genelist_2$Figure2C %>% sapply(function(x) {
    str_split(x, "\n")[[1]][1]
})
write.csv(fig_genelist_2, "results/tables/Figure_2/fig_genelist_2.csv")

# Colors for each group
group_colors <- c(
    "Forebrain\nn=428" = "#4d6da5",
    "Dorsal forebrain\nn=463" = "#78588c",
    "Early differentiation\nn=585" = "#b16060",
    "Ventral forebrain\nn=524" = "#5e9a5e"
)

# Block annotation for color blocks (without outline)
clusters_WTC_ha <- rowAnnotation(
    clusters = anno_block(
        gp = gpar(fill = group_colors[levels(row_split)], col = NA), # Remove outline with `col = NA`
        which = "row",
        width = unit(2, "mm") # Thinner blocks
    ),
    labels = anno_block(
        gp = gpar(fill = "white", col = "white"), # Remove outline with `col = NA`
        labels = levels(row_split), # Add group labels
        labels_gp = gpar(fontsize = 10, fontface = "bold"), # Customize label appearance
        labels_rot = -90, # Horizontal labels
        labels_just = "center"
    )
)

WTC_ht_plot <- Heatmap(
    scaled_mat[clustering_WTC$order, sample_order_WTC],
    name = "Normalized expression",
    column_names_gp = gpar(fontsize = 6),
    row_split = row_split,
    row_title = NULL,
    cluster_rows = FALSE,
    right_annotation = clusters_WTC_ha,
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

png_save(
    plot = WTC_ht_plot,
    path = "results/images/Figure_2/Figure2C.png",
    width = 2000,
    height = 1800
)

#  Making heatmap for WTC lineage only
sample_order_LON <- c(
    filter(meta, type == "ventral" & line == "LON71")$sample[order(filter(meta, type == "ventral" & line == "LON71")$day, decreasing = TRUE)],
    filter(meta, type == "dorsal" & line == "LON71")$sample[order(filter(meta, type == "dorsal" & line == "LON71")$day)]
)

clustering_LON <- hclust(dist(scaled_mat[, sample_order_LON]))
clusters_LON <- cutree(clustering_LON, k = 4)
# Check out cluster sizes and order
row_split <- factor(
    clusters_LON[clustering_LON$order],
    levels = c(1, 2, 3, 4),
    labels = c(
        "Early differentiation\nn=553",
        "Forebrain\nn=443",
        "Ventral forebrain\nn=522",
        "Dorsal forebrain\nn=482"
    )
)

# Colors for each group
group_colors <- c(
    "Forebrain\nn=443" = "#4d6da5",
    "Dorsal forebrain\nn=482" = "#78588c",
    "Early differentiation\nn=553" = "#b16060",
    "Ventral forebrain\nn=522" = "#5e9a5e"
)

# Block annotation for color blocks (without outline)
clusters_LON_ha <- rowAnnotation(
    clusters = anno_block(
        gp = gpar(fill = group_colors[levels(row_split)], col = NA), # Remove outline with `col = NA`
        which = "row",
        width = unit(2, "mm") # Thinner blocks
    ),
    labels = anno_block(
        gp = gpar(fill = "white", col = "white"), # Remove outline with `col = NA`
        labels = levels(row_split), # Add group labels
        labels_gp = gpar(fontsize = 10, fontface = "bold"), # Customize label appearance
        labels_rot = -90, # Horizontal labels
        labels_just = "center"
    )
)

LON_ht_plot <- Heatmap(
    scaled_mat[clustering_LON$order, sample_order_LON],
    name = "Normalized expression",
    column_names_gp = gpar(fontsize = 6),
    row_split = row_split,
    row_title = NULL,
    cluster_rows = FALSE,
    right_annotation = clusters_LON_ha,
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

png_save(
    plot = LON_ht_plot,
    path = "results/images/Figure_2/Figure2suppB.png",
    width = 2000,
    height = 1800
)

row_split <- factor(
    clusters_WTC[clustering_WTC$order],
    levels = c(1, 2, 3, 4),
    labels = c(
        "Early differentiation\nn=585",
        "Forebrain\nn=428",
        "Ventral forebrain\nn=524",
        "Dorsal forebrain\nn=463"
    )
)
# GO enrichment for WTC lineage only
GO_1 <- plot_go_term(
    names(clusters_WTC[which(clusters_WTC == 1)]),
    path = "results/images/Figure_2/Figure2supp_WTC_earlydiff.png",
    range = c(1:5),
    imgh = 8,
)
write.csv(GO_1, "results/tables/Figure_2/WTC_earlydiff_GO.csv")
GO_2 <- plot_go_term(
    names(clusters_WTC[which(clusters_WTC == 2)]),
    path = "results/images/Figure_2/Figure2supp_WTC_forebrain.png",
    range = c(1:5),
    imgh = 8,
)
write.csv(GO_2, "results/tables/Figure_2/WTC_forebrain_GO.csv")
GO_3 <- plot_go_term(
    names(clusters_WTC[which(clusters_WTC == 3)]),
    path = "results/images/Figure_2/Figure2supp_WTC_vAN.png",
    range = c(1:5),
    imgh = 8,
)
write.csv(GO_3, "results/tables/Figure_2/WTC_vAN_GO.csv")
GO_4 <- plot_go_term(
    names(clusters_WTC[which(clusters_WTC == 4)]),
    path = "results/images/Figure_2/Figure2supp_WTC_dAN.png",
    range = c(1:5),
    imgh = 8,
)
write.csv(GO_4, "results/tables/Figure_2/WTC_dAN_GO.csv")
# GO enrichment for LON lineage only
GO_1b <- plot_go_term(
    names(clusters_LON[which(clusters_LON == 1)]),
    path = "results/images/Figure_2/Figure2supp_LON_earlydiff.png",
    range = c(1:5),
    imgh = 8,
)
write.csv(GO_1b, "results/tables/Figure_2/LON_earlydiff_GO.csv")
GO_2b <- plot_go_term(
    names(clusters_LON[which(clusters_LON == 2)]),
    path = "results/images/Figure_2/Figure2supp_LON_forebrain.png",
    range = c(1:5),
    imgh = 8,
)
write.csv(GO_2b, "results/tables/Figure_2/LON_forebrain_GO.csv")
GO_3b <- plot_go_term(
    names(clusters_LON[which(clusters_LON == 3)]),
    path = "results/images/Figure_2/Figure2supp_LON_vAN.png",
    range = c(1:5),
    imgh = 8,
)
write.csv(GO_3b, "results/tables/Figure_2/LON_vAN_GO.csv")
GO_4b <- plot_go_term(
    names(clusters_LON[which(clusters_LON == 4)]),
    path = "results/images/Figure_2/Figure2supp_LON_dAN.png",
    range = c(1:5),
    imgh = 8,
)
write.csv(GO_4b, "results/tables/Figure_2/LON_dAN_GO.csv")

# making table of genes used in heatmap with cluster and subcluster annotation
genes_cluster <- data.frame(
    genes = hm_genes %>% gene_converter("ENSEMBL", "SYMBOL"),
    cluster_LON = clusters_LON[hm_genes],
    cluster_WTC = clusters_WTC[hm_genes],
)
rownames(genes_cluster) <- hm_genes

write.csv(genes_cluster, "results/tables/Figure_2/WTC_LON_heatmap_genes.csv")

######################################
######################################
# FIGURE 2 D-E : vAN vs dAN at day02

# DEGs ventral VS dorsal at day02
DEGs_day02 <- DESeqDataSetFromMatrix(
    countData = counts[, filter(meta, day == "day02")$sample],
    colData = filter(meta, day == "day02"),
    design = ~ line + type
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("type", "ventral", "dorsal")) %>%
    as.data.frame() %>%
    na.omit()
DEGs_day02$gene <- gene_converter(rownames(DEGs_day02), "ENSEMBL", "SYMBOL")
DEGs_day02 <- filter(DEGs_day02, !is.na(gene))
DEGs_day02$f1 <- ifelse((abs(DEGs_day02$log2FoldChange) >= 1 & DEGs_day02$padj < 0.01), DEGs_day02$gene, NA)
DEGs_day02$f2 <- ifelse((abs(DEGs_day02$log2FoldChange) >= 2 & DEGs_day02$padj < 0.01), DEGs_day02$gene, NA)

# Building gene list for Figure 2D
fig_genelist_3 <- data.frame(Figure2D = rep("Higly_DE", length(rownames(DEGs_day02)[which(!is.na(DEGs_day02$f2))])))
rownames(fig_genelist_3) <- rownames(DEGs_day02)[which(!is.na(DEGs_day02$f2))]
write.csv(fig_genelist_3, "results/tables/Figure_2/fig_genelist_3.csv")

vplot <- ggplot(DEGs_day02, aes(x = log2FoldChange, y = -log10(padj), label = f2)) +
    geom_point(data = filter(DEGs_day02, !is.na(f2)), aes(x = log2FoldChange, y = -log10(padj)), color = "red", size = 4) +
    geom_point(data = filter(DEGs_day02, is.na(f2)), aes(x = log2FoldChange, y = -log10(padj)), alpha = 0.5, color = "grey", size = 4) +
    geom_text_repel(size = 9, fontface = "bold") +
    xlab("log2(FoldChange)") +
    ylab("-log10(adjusted pvalue)") +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
    theme(
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 20),
    ) +
    custom_theme()
ggsave("results/images/Figure_2/Figure2D.png", plot = vplot, width = 12, height = 10)

######################################
######################################
# FIGURE 2 E : lineplots for selected genes

lp_meta <- filter(rawmeta, sample == "WTC6cipc" | (sample != "LON71_D12_2" & diff == "diff13" & line == "WTC" & type == "ventral" & ((manip == "veranika" & day != "day12") | (manip == "lauryane" & day == "day12"))))
lp_meta[which(lp_meta$sample == "WTC6cipc"), "day"] <- "day0"
# filtering out lowly expressed genes
lp_counts <- rawcounts[, lp_meta$sample][which(rowSums(rawcounts[, lp_meta$sample]) >= 25), ]
# lp_counts <- rawcounts[hm_genes, lp_meta$sample]
lp_vsd <- DESeqDataSetFromMatrix(
    countData = lp_counts,
    colData = lp_meta,
    design = ~day
) %>% vst(blind = FALSE)

genes1 <- c("FOXA2", "PTCH1", "SIX3", "SHH", "GSC", "LRP2", "CHRD", "FREM1")
genes2 <- c("DKK1", "SOX5", "NKX2-1", "SLIT2", "FGF10", "DDC")
genes3 <- c("GRM3", "NOG", "SOX6", "NTNG1", "EDN3", "CTNNA3", "NEDD9", "CLSTN2", "SLC8A1", "PITX2", "TMEFF2", "SIM1", "GADL1", "KCND3", "LRRK2")

symbol_vsd <- assay(lp_vsd)
rownames(symbol_vsd) <- gene_converter(rownames(symbol_vsd), "ENSEMBL", "SYMBOL")
filtered_vsd <- symbol_vsd[Reduce(union, list(genes1, genes2, genes3)), ]
scaled_filtered_vsd <- filtered_vsd - min(filtered_vsd)

genes1_df_long <- scaled_filtered_vsd[genes1, ] %>%
    t() %>%
    cbind(dplyr::select(lp_meta, c("day")), .[lp_meta$sample, ]) %>%
    reshape2::melt(., id.vars = "day", variable.name = "gene", value.name = "expression") %>%
    group_by(., gene, day) %>%
    summarise(., mean_expression = mean(expression), sd_expression = sd(expression))
genes1_df_long$sd_expression[which(is.na(genes1_df_long$sd_expression))] <- 0
genes1_plot <- kinetic_lineplots(genes1_df_long)
ggsave("results/images/Figure_2/Figure2E1.png", genes1_plot, width = 15, height = 10)

genes2_df_long <- scaled_filtered_vsd[genes2, ] %>%
    t() %>%
    cbind(dplyr::select(lp_meta, c("day")), .[lp_meta$sample, ]) %>%
    reshape2::melt(., id.vars = "day", variable.name = "gene", value.name = "expression") %>%
    group_by(., gene, day) %>%
    summarise(., mean_expression = mean(expression), sd_expression = sd(expression))
genes2_df_long$sd_expression[which(is.na(genes2_df_long$sd_expression))] <- 0
genes2_plot <- kinetic_lineplots(genes2_df_long)
ggsave("results/images/Figure_2/Figure2E2.png", genes2_plot, width = 15, height = 10)

genes3_df_long <- scaled_filtered_vsd[genes3, ] %>%
    t() %>%
    cbind(dplyr::select(lp_meta, c("day")), .[lp_meta$sample, ]) %>%
    reshape2::melt(., id.vars = "day", variable.name = "gene", value.name = "expression") %>%
    group_by(., gene, day) %>%
    summarise(., mean_expression = mean(expression), sd_expression = sd(expression))
genes3_df_long$sd_expression[which(is.na(genes3_df_long$sd_expression))] <- 0
genes3_plot <- kinetic_lineplots(genes3_df_long)
ggsave("results/images/Figure_2/Figure2E3.png", genes3_plot, width = 15, height = 10)

######################################
######################################
# FIGURE 2 G : Expression pattern plots

DEGs <- read.csv("/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/Volcano_DEG_dbd_ventral.csv", header = TRUE)$gene %>% gene_converter("SYMBOL", "ENSEMBL")

lp_meta <- filter(rawmeta, (sample != "LON71_D12_2" & diff == "diff13" & line == "WTC" & type == "ventral" & ((manip == "veranika" & day != "day12") | (manip == "lauryane" & day == "day12"))))
rownames(lp_meta) <- lp_meta$sample
lp_meta$day <- as.factor(lp_meta$day)

# filtering out lowly expressed genes
lp_counts <- rawcounts[, lp_meta$sample][which(rowSums(rawcounts[, lp_meta$sample]) >= 25), ]
# making DESeq object with days as covariate
lp_vsd <- DESeqDataSetFromMatrix(
    countData = lp_counts,
    colData = lp_meta,
    design = ~day
) %>% vst(blind = FALSE)

clusters <- degPatterns(
    assay(lp_vsd)[rownames(lp_vsd) %in% DEGs, ],
    meta = lp_meta,
    time = "day",
    reduce = TRUE,
    nClusters = 10,
)
# intervert cluster 1 with 2 and 6 with 1 for convenience
clusters$normalize$cluster <- ifelse(clusters$normalize$cluster == 2, 6, ifelse(clusters$normalize$cluster == 6, 2, clusters$normalize$cluster))
clusters$df$cluster <- ifelse(clusters$df$cluster == 2, 6, ifelse(clusters$df$cluster == 6, 2, clusters$df$cluster))

fig_genelist_4 <- clusters$df %>%
    filter(cluster %in% c(1, 2)) %>%
    dplyr::select(cluster)
colnames(fig_genelist_4) <- "Figure2H"
write.csv(fig_genelist_4, "results/tables/Figure_2/fig_genelist_4.csv")

sign_comp <- list(
    c("day02", "day04"),
    c("day04", "day06"),
    c("day06", "day08"),
    c("day08", "day10"),
    c("day10", "day12")
)

plot_list <- lapply(1:9, function(i) {
    return(MyDegPlotCluster(table = filter(clusters$normalize, cluster == i), time = "day", sign_comp = sign_comp, cluster_i = i))
})
combined_plot <- wrap_plots(plot_list, ncol = 3)
ggsave(paste0("results/images/Figure_2/Figure2G1.png"), combined_plot, width = 35, height = 25)
ggsave(paste0("results/images/Figure_2/Figure2G2.png"), plot_list[[1]], width = 15, height = 10)
ggsave(paste0("results/images/Figure_2/Figure2G3.png"), plot_list[[2]], width = 15, height = 10)

clusters$df$symbol <- clusters$df$genes %>% gene_converter("ENSEMBL", "SYMBOL")
write.csv(clusters$df, file = "/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_2/DEGpattern_WTC.csv", quote = FALSE, row.names = FALSE)

clusters$df
