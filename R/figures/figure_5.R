# Loading packages and functions
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(DESeq2)
library(Seurat)

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
cellids <- read.table("scdata/indices_week345_wholebody.csv", header = TRUE, sep = ",")
genes <- read.table("scdata/genes_week345_wholebody.tsv", header = FALSE)$V1

# set rownames
rownames(sc_counts) <- genes

# format and set colnames
formated_indices <- c(1:nrow(cellids)) %>% sapply(function(i) {
    c(str_split(cellids[i, ]$index, "-")[[1]][1:2], cellids[i, ]$week_stage) %>%
        paste(collapse = "-") %>%
        return()
})
colnames(sc_counts) <- formated_indices

sc_meta <- read.csv("scdata/week345_wholebody_metadata.csv") %>%
    filter(celltype_region_num %in% c(1:8)) %>%
    dplyr::select(c("week_stage", "celltype_region_num", "celltype_region", "barcode"))

sc_meta$celltype_region[which(sc_meta$celltype_region == "Spinal Cord Motor Neuron 1")] <- "SC Motor Neuron 1"
sc_meta$celltype_region[which(sc_meta$celltype_region == "Spinal Cord Motor Neuron 2")] <- "SC Motor Neuron 2"
sc_meta$celltype_region[which(sc_meta$celltype_region == "Spinal Cord Motor Neuron 3")] <- "SC Motor Neuron 3"

sc_meta_weekAll <- sc_meta
seurat_weekAll <- CreateSeuratObject(counts = sc_counts[, sc_meta_weekAll$barcode], meta.data = sc_meta_weekAll, min.cells = 20)
seurat_weekAll <- NormalizeData(seurat_weekAll)

seurat_weekAll <- FindVariableFeatures(seurat_weekAll, selection.method = "vst", nfeatures = 2000)
seurat_weekAll <- ScaleData(seurat_weekAll, features = rownames(seurat_weekAll))
seurat_weekAll <- RunPCA(seurat_weekAll, features = VariableFeatures(object = seurat_weekAll))
seurat_weekAll <- RunUMAP(seurat_weekAll, dims = 1:10)

umap_weekAll <- seurat_weekAll[["umap"]]@cell.embeddings %>% as.data.frame()
sc_meta_weekAll <- cbind(umap_weekAll[sc_meta_weekAll$barcode, ], sc_meta_weekAll)


sc_meta_week3 <- filter(sc_meta, week_stage == "W3-1")
seurat_week3 <- CreateSeuratObject(counts = sc_counts[, sc_meta_week3$barcode], meta.data = sc_meta_week3, min.cells = 20)
seurat_week3 <- NormalizeData(seurat_week3)

seurat_week3 <- FindVariableFeatures(seurat_week3, selection.method = "vst", nfeatures = 2000)
seurat_week3 <- ScaleData(seurat_week3, features = rownames(seurat_week3))
seurat_week3 <- RunPCA(seurat_week3, features = VariableFeatures(object = seurat_week3))
seurat_week3 <- RunUMAP(seurat_week3, dims = 1:10)

umap_week3 <- seurat_week3[["umap"]]@cell.embeddings %>% as.data.frame()
sc_meta_week3 <- cbind(umap_week3[sc_meta_week3$barcode, ], sc_meta_week3)
View(sc_meta_week3)

sc_meta_week4 <- filter(sc_meta, week_stage == "W4-1")
seurat_week4 <- CreateSeuratObject(counts = sc_counts[, sc_meta_week4$barcode], meta.data = sc_meta_week4, min.cells = 20)
seurat_week4 <- NormalizeData(seurat_week4)

seurat_week4 <- FindVariableFeatures(seurat_week4, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_week4)
seurat_week4 <- ScaleData(seurat_week4, features = all.genes)
seurat_week4 <- RunPCA(seurat_week4, features = VariableFeatures(object = seurat_week4))
seurat_week4 <- RunUMAP(seurat_week4, dims = 1:10)

umap_week4 <- seurat_week4[["umap"]]@cell.embeddings %>% as.data.frame()
sc_meta_week4 <- cbind(umap_week4[sc_meta_week4$barcode, ], sc_meta_week4)



sc_meta_week5 <- filter(sc_meta, week_stage == "W5-1")
seurat_week5 <- CreateSeuratObject(counts = sc_counts[, sc_meta_week5$barcode], meta.data = sc_meta_week5, min.cells = 20)
seurat_week5 <- NormalizeData(seurat_week5)

seurat_week5 <- FindVariableFeatures(seurat_week5, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_week5)
seurat_week5 <- ScaleData(seurat_week5, features = all.genes)
seurat_week5 <- RunPCA(seurat_week5, features = VariableFeatures(object = seurat_week5))
seurat_week5 <- RunUMAP(seurat_week5, dims = 1:10)

umap_week5 <- seurat_week5[["umap"]]@cell.embeddings %>% as.data.frame()
sc_meta_week5 <- cbind(umap_week5[sc_meta_week5$barcode, ], sc_meta_week5)

umap_plot <- ggplot(data = sc_meta_weekAll, aes(x = umap_1, y = umap_2, color = celltype_region)) +
    geom_point(size = 2.5) +
    guides(color = guide_legend(
        override.aes = list(size = 5)
    )) +
    labs(color = "Cell identity") + # Modify the legend title # Legend point size
    theme(
        legend.title = element_text(size = 25), # Legend title size
        legend.text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)
    ) +
    custom_theme()
ggsave(paste0("results/images/Figure_5/Zeng_weekAll_celltypes", ".png"), plot = umap_plot, width = 14, height = 10, dpi = 300)

umap_plot <- ggplot(data = sc_meta_weekAll, aes(x = umap_1, y = umap_2, color = week_stage)) +
    geom_point(size = 2.5) +
    guides(color = guide_legend(
        override.aes = list(size = 5)
    )) +
    labs(color = "Week") + # Modify the legend title # Legend point size
    custom_theme() +
    theme(
        legend.text = element_text(size = 15), # Change legend text size
        legend.title = element_text(size = 20), # Change legend title size
        legend.key.size = unit(1, "cm"),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)
    )
ggsave(paste0("results/images/Figure_5/Zeng_weekAll_weekstage", ".png"), plot = umap_plot, width = 14, height = 10, dpi = 300)


# UMAP colored by cellular types :
umap_plot <- ggplot(data = sc_meta_week3, aes(x = umap_1, y = umap_2, color = celltype_region)) +
    geom_point(size = 1.5) +
    custom_theme()
ggsave(paste0("results/images/Figure_5/Zeng_week3_celltypes", ".png"), plot = umap_plot, width = 14, height = 10, dpi = 300)

umap_plot <- ggplot(data = sc_meta_week4, aes(x = umap_1, y = umap_2, color = celltype_region)) +
    geom_point(size = 1.5) +
    scale_color_manual(values = c("#008e00", "#fb29ff", "#ab6f00", "#7d0000", "#00b7c0", "#fc1f1f", "purple", "#000091")) +
    guides(color = guide_legend(override.aes = list(size = 6))) +
    custom_theme()
ggsave(paste0("results/images/Figure_5/Zeng_week4_celltypes", ".png"), plot = umap_plot, width = 14, height = 10, dpi = 300)

umap_plot <- ggplot(data = sc_meta_week5, aes(x = umap_1, y = umap_2, color = celltype_region)) +
    geom_point(size = 1.5) +
    custom_theme()
ggsave(paste0("results/images/Figure_5/Zeng_week5_celltypes", ".png"), plot = umap_plot, width = 14, height = 10, dpi = 300)

View(as.data.frame(rownames(seurat_week3)))
View(as.data.frame(rownames(seurat_week4)))
View(as.data.frame(rownames(seurat_week5)))

blue_clust_df_f <- read.csv("results/tables/Figure_3/SHH_cluster.csv")
blue_clust_df_f %>% View()
marker_genes <- blue_clust_df_f$gene
marker_genes
gc()

genes_ventral <- c("FOXA2", "FREM1", "SHH", "NKX2-1", "PTCH1", "LINC00261", "PLCL1", "CAPN6", "LRRK2", "DDC", "SMIM32", "GADL1", "DRC1", "CLSTN2", "SLIT2")
genes_dorsal <- c("CNTN2", "PAX6", "PAX3", "CNTNAP2", "GLI3", "EMX2", "NELL2", "GDF7", "SYT4", "GRIP2")
genes_kinetic <- c("PTCH1", "FOXA2", "FREM1", "SHH", "LRP2", "SIX3", "FGF10", "NKX2-1", "FGF9", "NKX2-2", "GADL1", "CLSTN2", "DRC1")
genes_bigger <- c("PTCH1", "NKX2-1", "DRC1", "CAPN6", "CNTNAP2", "FREM1", "SIX3", "CLSTN2", "SHH", "GADL1")

sc_meta_week3 %>% head()

marker_genes <- genes_bigger
for (gene in marker_genes) {
    if (gene %in% rownames(seurat_week4)) {
        # Extract expression data
        sc_meta_week4$expression <- seurat_week4@assays$RNA$scale.data[gene, sc_meta_week4$barcode]

        sc_meta_week4 <- sc_meta_week4 %>% arrange(expression)

        # Create UMAP plot
        umap_plot <- ggplot(data = sc_meta_week4, aes(x = umap_1, y = umap_2, color = expression)) +
            geom_point(size = 2.5) +
            ggtitle(gene) +
            scale_color_gradient(low = "lightblue", high = "darkred") +
            labs(color = "Expression") +
            custom_theme() +
            theme(
                legend.text = element_text(size = 15), # Change legend text size
                legend.title = element_text(size = 20), # Change legend title size
                legend.key.size = unit(1, "cm"),
                axis.text.x = element_text(size = 20),
                axis.text.y = element_text(size = 20)
            )

        # Save the plot
        ggsave(paste0("/home/jules/Documents/phd/projects/panasenkava_2024/results/images/Figure_5/Zeng_week4_", gene, ".png"),
            plot = umap_plot, width = 14, height = 10, dpi = 300
        )

        # Garbage collection
        gc()
    } else {
        print(paste0(gene, " not in week4"))
    }
}

test <- AddMetaData(seurat_weekAll, sc_meta_weekAll)

DoHeatmap(seurat_weekAll, group.by = "region_type", features = genes_kinetic) + NoLegend()



png(filename = "results/images/Figure_2A/weekALL_kinetic_genes.png", width = 2400, height = 1600, res = 20)
Heatmap(
    scaled_mat,
    name = "Normalized expression",
    column_names_gp = gpar(fontsize = 6),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    row_names_side = "left",
    show_column_names = FALSE,
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










for (gene in marker_genes) {
    if (gene %in% rownames(seurat_week3)) {
        sc_meta_week3$expression <- seurat_week3@assays$RNA$scale.data[gene, sc_meta_week3$barcode]
        umap_plot <- ggplot(data = sc_meta_week3, aes(x = umap_1, y = umap_2, color = expression)) +
            geom_point(size = 1) +
            ggtitle(gene) +
            scale_color_gradient(low = "lightblue", high = "darkred") +
            custom_theme()
        ggsave(paste0("/home/jules/Documents/phd/projects/panasenkava_2024/results/images/Figure_2A/umaps/Zeng_week3_", gene, ".png"), plot = umap_plot, width = 14, height = 10, dpi = 300)
        gc()
    } else {
        print(paste0(gene, " not in week3"))
    }
}

gc()
marker_genes <- c("SFTA3", "CAPN6", "EPHB1", "ZIC5", "AFF2", "NUAK2")
for (gene in marker_genes) {
    if (gene %in% rownames(seurat_week4)) {
        sc_meta_week4$expression <- seurat_week4@assays$RNA$scale.data[gene, sc_meta_week4$barcode]
        umap_plot <- ggplot(data = sc_meta_week4, aes(x = umap_1, y = umap_2, color = expression)) +
            geom_point(size = 1) +
            ggtitle(gene) +
            scale_color_gradient(low = "lightblue", high = "darkred") +
            custom_theme()
        ggsave(paste0("/home/jules/Documents/phd/projects/panasenkava_2024/results/images/Zeng_week4_", gene, ".png"), plot = umap_plot, width = 14, height = 10, dpi = 300)
        gc()
    } else {
        print(paste0(gene, " not in week4"))
    }
}

gc()
for (gene in marker_genes) {
    if (gene %in% rownames(seurat_week5)) {
        sc_meta_week5$expression <- seurat_week5@assays$RNA$scale.data[gene, sc_meta_week5$barcode]
        umap_plot <- ggplot(data = sc_meta_week5, aes(x = umap_1, y = umap_2, color = expression)) +
            geom_point(size = 1) +
            ggtitle(gene) +
            scale_color_gradient(low = "lightblue", high = "darkred") +
            custom_theme()
        ggsave(paste0("/home/jules/Documents/phd/projects/panasenkava_2024/results/images/Figure_2A/umaps/Zeng_week5_", gene, ".png"), plot = umap_plot, width = 14, height = 10, dpi = 300)
        gc()
    } else {
        print(paste0(gene, " not in week5"))
    }
}










focusgenes <- c("FOXG1", "CNTN6", "ZIC2", "PAX6", "FOXA2", "LINC00261", "JAG1", "QKI", "NUAK2", "TRM9", "EPHB1", "ZIC5", "NKX2-1", "SHH", "CNTNAP2", "CAPN6", "ZFHX4")
focusgenes <- marker_genes
focusgenes
focusduo <- combn(focusgenes, 2, simplify = FALSE)

for (duo in focusduo) {
    if ((!duo[1] %in% rownames(seurat_week4)) | (!duo[2] %in% rownames(seurat_week4))) {
        print("These genes are not in week4:")
        if (!duo[1] %in% rownames(seurat_week4)) {
            print(paste0("    ", duo[1]))
        }
        if (!duo[2] %in% rownames(seurat_week4)) {
            print(paste0("    ", duo[2]))
        }
        next
    }
    rawplots <- FeaturePlot(object = seurat_week4, pt.size = 2, order = TRUE, features = duo, blend = TRUE)
    ftplot <- rawplots[[3]] + rawplots[[4]] + plot_layout(widths = c(6, 1))
    ggsave(paste0("Figure_5/duo_week4/Zeng_week4", paste(duo, collapse = "_"), ".png"), plot = ftplot, width = 14, height = 10, dpi = 300)
}




for (duo in focusduo) {
    duo <- focusduo[[1]]
    if ((!duo[1] %in% rownames(seurat_week3)) | (!duo[2] %in% rownames(seurat_week3))) {
        print("These genes are not in week3:")
        if (!duo[1] %in% rownames(seurat_week3)) {
            print(paste0("    ", duo[1]))
        }
        if (!duo[2] %in% rownames(seurat_week3)) {
            print(paste0("    ", duo[2]))
        }
        next
    }
    rawplots <- FeaturePlot(object = seurat_week4, pt.size = 2, blend.threshold = 0.99, order = TRUE, features = c("FOXA2", "LINC00261"), blend = TRUE)
    ftplot <- rawplots[[3]] + rawplots[[4]] + plot_layout(widths = c(6, 1))
    ftplot
    ggsave(paste0("Figure_5/duo_week3/Zeng_week3", paste(duo, collapse = "_"), ".png"), plot = ftplot, width = 14, height = 10, dpi = 300)
}

for (duo in focusduo) {
    if ((!duo[1] %in% rownames(seurat_week5)) | (!duo[2] %in% rownames(seurat_week5))) {
        print("These genes are not in week5:")
        if (!duo[1] %in% rownames(seurat_week5)) {
            print(paste0("    ", duo[1]))
        }
        if (!duo[2] %in% rownames(seurat_week5)) {
            print(paste0("    ", duo[2]))
        }
        next
    }
    rawplots <- FeaturePlot(object = seurat_week5, pt.size = 2, order = TRUE, features = duo, blend = TRUE)
    ftplot <- rawplots[[3]] + rawplots[[4]] + plot_layout(widths = c(6, 1))
    ggsave(paste0("Figure_5/duo_week5/Zeng_week5", paste(duo, collapse = "_"), ".png"), plot = ftplot, width = 14, height = 10, dpi = 300)
}


# CHECK EXPRESSION cinétique DANS LES UMAPS




library(data.table)
percents <- marker_genes %>%
    lapply(function(gene) {
        if (gene %in% rownames(seurat_weekAll@assays$RNA$counts)) {
            percent <- unique(sc_meta_weekAll$week_stage) %>% sapply(function(ct) {
                test <- seurat_weekAll@assays$RNA$counts[gene, ]
                test <- ifelse(test >= 1, sc_meta_weekAll$week_stage, "none")
                test2 <- table(test)
                ncell <- test2[ct]
                if (is.na(ncell)) {
                    ncell <- 0
                }
                return(unname(ncell / nrow(filter(sc_meta_weekAll, week_stage == ct)) * 100))
            })
        } else {
            percent <- rep(0, 8)
        }
        print(gene)
        return(percent)
    }) %>%
    do.call(rbind, .) %>%
    as.data.frame()
rownames(percents) <- marker_genes
View(percents)
percents_tot <- percents
percents_tot$total <- rowSums(percents) / 3
View(percents_tot)

write.csv(percents_tot, file = "results/images/Figure_5/percents.csv", quote = FALSE)

percents_tot <- read.csv("/home/jules/Documents/phd/projects/panasenkava_2024/Figure_5/SHH_cluster_500_percents.csv", row.names = 1)



percents_tot_pos <- percents_tot[filter(top500, cor > 0)$gene, ]
write.csv(percents_tot_pos, file = "/home/jules/Documents/phd/projects/panasenkava_2024/Figure_5/SHH_cluster_500_percents_pos.csv", quote = FALSE)

percents_tot_neg <- percents_tot[filter(top500, cor < 0)$gene, ]
write.csv(percents_tot_neg, file = "/home/jules/Documents/phd/projects/panasenkava_2024/Figure_5/SHH_cluster_500_percents_neg.csv", quote = FALSE)


View(percents_tot)



sc_meta_week4$celltype_region %>% table()


top500 <- read.csv("/home/jules/Documents/phd/projects/panasenkava_2024/results/tables/Figure_3/SHH_cluster_500.csv")
top500 %>% head()
w4_cellmatch <- sc_meta_week4

marker_genes <- c("SFTA3", "CAPN6", "EPHB1", "AFF2", "CNTNAP2", "ZIC5", "NUAK2", "SPON1", "SFRP1", "SHROOM3", "SALL1")

setdiff(marker_genes, rownames(seurat_week4@assays$RNA$counts))
marker_genes
w4df <- marker_genes %>%
    lapply(function(gene) {
        if (gene %in% rownames(seurat_week4@assays$RNA$counts)) {
            return(seurat_week4@assays$RNA$counts[gene, ] >= 1)
        } else {
            print(paste0(gene, " not in week4"))
            return(rep(FALSE, ncol(seurat_week4@assays$RNA$data)))
        }
    }) %>%
    do.call(cbind, .)
colnames(w4df) <- marker_genes

w4_cellmatch <- cbind(w4_cellmatch, w4df[w4_cellmatch$barcode, ])
w4_cellmatch
which(seurat_week4@assays$RNA$counts["NKX2-1", ] >= 1)
seurat_week4@assays$RNA$counts["NKX2-1", ]
w4df <- union(marker_genes, c("NKX2-1", "PAX6", "FOXA2", "LHX2")) %>%
    lapply(function(gene) {
        if (gene %in% rownames(seurat_week4@assays$RNA$counts)) {
            return(seurat_week4@assays$RNA$counts[gene, ] > 1)
        } else {
            print(paste0(gene, " not in week4"))
            return(rep(FALSE, ncol(seurat_week4@assays$RNA$data)))
        }
    }) %>%
    do.call(cbind, .)
colnames(w4df) <- union(marker_genes, c("NKX2-1", "PAX6", "FOXA2", "LHX2"))

w4_cellmatch <- cbind(w4_cellmatch, w4df[w4_cellmatch$barcode, ])

pairwise_df <- expand.grid(marker_genes, c("NKX2-1", "PAX6", "FOXA2", "LHX2"))
pairwise_df$Var1 <- as.vector(pairwise_df$Var1)
pairwise_df$Var2 <- as.vector(pairwise_df$Var2)

pairwise_df %>% head()
w4_cellmatch[, pairwise_df$Var1[1]] %>% table()
w4_cellmatch[, pairwise_df$Var2[1]] %>% table()
which(w4_cellmatch[, "NKX2-1"] == TRUE)
w4_cellmatch$expression
# execute this !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
for (i in c(1:nrow(pairwise_df))) {
    gene1 <- pairwise_df$Var1[i]
    gene2 <- pairwise_df$Var2[i]
    w4_cellmatch$gene1 <- w4_cellmatch[, pairwise_df$Var1[i]]
    w4_cellmatch$gene2 <- w4_cellmatch[, pairwise_df$Var2[i]]
    w4_cellmatch$duo <- sapply(c(1:nrow(w4_cellmatch)), function(i) {
        if (w4_cellmatch$gene1[i] == TRUE & w4_cellmatch$gene2[i] == TRUE) {
            return("1-both")
        } else if (w4_cellmatch$gene1[i] == TRUE) {
            return(paste0("2-", gene1))
        } else if (w4_cellmatch$gene2[i] == TRUE) {
            return(paste0("3-", gene2))
        } else {
            return("4-none")
        }
    })
    w4_cellmatch$expression <- w4_cellmatch$duo
    w4_cellmatch <- w4_cellmatch[order(w4_cellmatch$expression, decreasing = TRUE), ]

    # w4_cellmatch[, c("umap1", "umap2", "expression")] %>% head()
    colors <- c("#ffdd00", "grey", "#ff0000", "green")
    names(colors) <- c("1-both", "4-none", paste0("2-", pairwise_df$Var1[i]), paste0("3-", pairwise_df$Var2[i]))
    plot <- ggplot(data = w4_cellmatch[, c("umap_1", "umap_2", "expression", "gene1", "gene2")], aes(x = umap_1, y = umap_2, color = expression)) +
        geom_point(size = 2.5) +
        ggtitle(paste0("Week4: gene1 = ", gene1, " and gene2 = ", gene2)) +
        scale_color_manual(values = colors) +
        custom_theme() +
        theme(
            legend.text = element_text(size = 15), # Change legend text size
            legend.title = element_text(size = 20), # Change legend title size
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20)
        )
    ggsave(paste0("results/images/Figure_5/final/", gene1, "_", gene2, ".png"), plot = plot, width = 14, height = 10, dpi = 300)
}


for (i in c(1:nrow(pairwise_df))) {
    gene1 <- pairwise_df$Var1[i]
    gene2 <- pairwise_df$Var2[i]
    w4_cellmatch$gene1 <- w4_cellmatch[, pairwise_df$Var1[i]]
    w4_cellmatch$gene2 <- w4_cellmatch[, pairwise_df$Var2[i]]
    w4_cellmatch$duo <- sapply(c(1:nrow(w4_cellmatch)), function(i) {
        if (w4_cellmatch$gene1[i] == TRUE & w4_cellmatch$gene2[i] == TRUE) {
            return("both")
        } else if (w4_cellmatch$gene1[i] == TRUE) {
            return(gene1)
        } else if (w4_cellmatch$gene2[i] == TRUE) {
            return(gene2)
        } else {
            return("none")
        }
    })
    w4_cellmatch$expression <- w4_cellmatch$duo
}













sc_meta_week3$GLI1 <- seurat_week3@assays$RNA$scale.data["GLI1", sc_meta_week3$barcode]
sc_meta_week3$GLI1 <- ifelse(sc_meta_week3$GLI1 > median(sc_meta_week3$GLI1), "GLI1", " ")
sc_meta_week3$GLI2 <- seurat_week3@assays$RNA$scale.data["GLI2", sc_meta_week3$barcode]
sc_meta_week3$GLI2 <- ifelse(sc_meta_week3$GLI2 > median(sc_meta_week3$GLI2), "GLI2", " ")
sc_meta_week3$GLI3 <- seurat_week3@assays$RNA$scale.data["GLI3", sc_meta_week3$barcode]
sc_meta_week3$GLI3 <- ifelse(sc_meta_week3$GLI3 > median(sc_meta_week3$GLI3), "GLI3", " ")

sc_meta_week4$GLI1 <- seurat_week4@assays$RNA$scale.data["GLI1", sc_meta_week4$barcode]
sc_meta_week4$GLI1 <- ifelse(sc_meta_week4$GLI1 > median(sc_meta_week4$GLI1), "GLI1", " ")
sc_meta_week4$GLI2 <- seurat_week4@assays$RNA$scale.data["GLI2", sc_meta_week4$barcode]
sc_meta_week4$GLI2 <- ifelse(sc_meta_week4$GLI2 > median(sc_meta_week4$GLI2), "GLI2", " ")
sc_meta_week4$GLI3 <- seurat_week4@assays$RNA$scale.data["GLI3", sc_meta_week4$barcode]
sc_meta_week4$GLI3 <- ifelse(sc_meta_week4$GLI3 > median(sc_meta_week4$GLI3), "GLI3", " ")

sc_meta_week5$GLI1 <- seurat_week5@assays$RNA$scale.data["GLI1", sc_meta_week5$barcode]
sc_meta_week5$GLI1 <- ifelse(sc_meta_week5$GLI1 > median(sc_meta_week5$GLI1), "GLI1", " ")
sc_meta_week5$GLI2 <- seurat_week5@assays$RNA$scale.data["GLI2", sc_meta_week5$barcode]
sc_meta_week5$GLI2 <- ifelse(sc_meta_week5$GLI2 > median(sc_meta_week5$GLI2), "GLI2", " ")
sc_meta_week5$GLI3 <- seurat_week5@assays$RNA$scale.data["GLI3", sc_meta_week5$barcode]
sc_meta_week5$GLI3 <- ifelse(sc_meta_week5$GLI3 > median(sc_meta_week5$GLI3), "GLI3", " ")

for (gene in marker_genes) {
    if (gene %in% rownames(seurat_week3)) {
        sc_meta_week3$expression <- seurat_week3@assays$RNA$scale.data[gene, sc_meta_week3$barcode]
        sc_meta_week3$expression <- ifelse(sc_meta_week3$expression > median(sc_meta_week3$expression), gene, " ")
        sc_meta_week3$GLI1_expression <- sapply(c(1:nrow(sc_meta_week3)), function(i) {
            if (sc_meta_week3$GLI1[i] == "GLI1" & sc_meta_week3$expression[i] == gene) {
                return(paste0("GLI1_", gene))
            } else if (sc_meta_week3$GLI1[i] == "GLI1") {
                return("GLI1")
            } else if (sc_meta_week3$expression[i] == gene) {
                return(gene)
            } else {
                return(" ")
            }
        })
        sc_meta_week3$GLI2_expression <- sapply(c(1:nrow(sc_meta_week3)), function(i) {
            if (sc_meta_week3$GLI2[i] == "GLI2" & sc_meta_week3$expression[i] == gene) {
                return(paste0("GLI2_", gene))
            } else if (sc_meta_week3$GLI2[i] == "GLI2") {
                return("GLI2")
            } else if (sc_meta_week3$expression[i] == gene) {
                return(gene)
            } else {
                return(" ")
            }
        })
        sc_meta_week3$GLI3_expression <- sapply(c(1:nrow(sc_meta_week3)), function(i) {
            if (sc_meta_week3$GLI3[i] == "GLI3" & sc_meta_week3$expression[i] == gene) {
                return(paste0("GLI3_", gene))
            } else if (sc_meta_week3$GLI3[i] == "GLI3") {
                return("GLI3")
            } else if (sc_meta_week3$expression[i] == gene) {
                return(gene)
            } else {
                return(" ")
            }
        })

        custom_colors <- c("#e9e9e9", "lightgreen", "#ffae51", "red")
        names(custom_colors) <- c(" ", "GLI1", gene, paste0("GLI1_", gene))
        umap_plot <- ggplot(data = sc_meta_week3, aes(x = umap_1, y = umap_2, color = GLI1_expression)) +
            geom_point(size = 1.5) +
            scale_color_manual(values = custom_colors) +
            custom_theme()
        ggsave(paste0("results/images/Figure_5/week3_GLI1/Zeng_week3_GLI1_", gene, ".png"), plot = umap_plot, width = 14, height = 10, dpi = 300)

        custom_colors <- c("#e9e9e9", "lightgreen", "#ffae51", "red")
        names(custom_colors) <- c(" ", "GLI2", gene, paste0("GLI2_", gene))
        umap_plot <- ggplot(data = sc_meta_week3, aes(x = umap_1, y = umap_2, color = GLI2_expression)) +
            geom_point(size = 1.5) +
            scale_color_manual(values = custom_colors) +
            custom_theme()
        ggsave(paste0("results/images/Figure_5/week3_GLI2/Zeng_week3_GLI2", gene, ".png"), plot = umap_plot, width = 14, height = 10, dpi = 300)

        custom_colors <- c("#e9e9e9", "lightgreen", "#ffae51", "red")
        names(custom_colors) <- c(" ", "GLI3", gene, paste0("GLI3_", gene))
        umap_plot <- ggplot(data = sc_meta_week3, aes(x = umap_1, y = umap_2, color = GLI3_expression)) +
            geom_point(size = 1.5) +
            scale_color_manual(values = custom_colors) +
            custom_theme()
        ggsave(paste0("results/images/Figure_5/week3_GLI3/Zeng_week3_GLI3", gene, ".png"), plot = umap_plot, width = 14, height = 10, dpi = 300)
        gc()
    } else {
        print(paste0(gene, " not in week3"))
    }
}

for (gene in marker_genes) {
    if (gene %in% rownames(seurat_week4)) {
        sc_meta_week4$expression <- seurat_week4@assays$RNA$scale.data[gene, sc_meta_week4$barcode]
        sc_meta_week4$expression <- ifelse(sc_meta_week4$expression > median(sc_meta_week4$expression), gene, " ")
        sc_meta_week4$GLI1_expression <- sapply(c(1:nrow(sc_meta_week4)), function(i) {
            if (sc_meta_week4$GLI1[i] == "GLI1" & sc_meta_week4$expression[i] == gene) {
                return(paste0("GLI1_", gene))
            } else if (sc_meta_week4$GLI1[i] == "GLI1") {
                return("GLI1")
            } else if (sc_meta_week4$expression[i] == gene) {
                return(gene)
            } else {
                return(" ")
            }
        })
        sc_meta_week4$GLI2_expression <- sapply(c(1:nrow(sc_meta_week4)), function(i) {
            if (sc_meta_week4$GLI2[i] == "GLI2" & sc_meta_week4$expression[i] == gene) {
                return(paste0("GLI2_", gene))
            } else if (sc_meta_week4$GLI2[i] == "GLI2") {
                return("GLI2")
            } else if (sc_meta_week4$expression[i] == gene) {
                return(gene)
            } else {
                return(" ")
            }
        })
        sc_meta_week4$GLI3_expression <- sapply(c(1:nrow(sc_meta_week4)), function(i) {
            if (sc_meta_week4$GLI3[i] == "GLI3" & sc_meta_week4$expression[i] == gene) {
                return(paste0("GLI3_", gene))
            } else if (sc_meta_week4$GLI3[i] == "GLI3") {
                return("GLI3")
            } else if (sc_meta_week4$expression[i] == gene) {
                return(gene)
            } else {
                return(" ")
            }
        })

        custom_colors <- c("#e9e9e9", "lightgreen", "#ffae51", "red")
        names(custom_colors) <- c(" ", "GLI1", gene, paste0("GLI1_", gene))
        umap_plot <- ggplot(data = sc_meta_week4, aes(x = umap_1, y = umap_2, color = GLI1_expression)) +
            geom_point(size = 1.5) +
            scale_color_manual(values = custom_colors) +
            custom_theme()
        ggsave(paste0("results/images/Figure_5/week4_GLI1/Zeng_week4_GLI1_", gene, ".png"), plot = umap_plot, width = 14, height = 10, dpi = 300)

        custom_colors <- c("#e9e9e9", "lightgreen", "#ffae51", "red")
        names(custom_colors) <- c(" ", "GLI2", gene, paste0("GLI2_", gene))
        umap_plot <- ggplot(data = sc_meta_week4, aes(x = umap_1, y = umap_2, color = GLI2_expression)) +
            geom_point(size = 1.5) +
            scale_color_manual(values = custom_colors) +
            custom_theme()
        ggsave(paste0("results/images/Figure_5/week4_GLI2/Zeng_week4_GLI2", gene, ".png"), plot = umap_plot, width = 14, height = 10, dpi = 300)

        custom_colors <- c("#e9e9e9", "lightgreen", "#ffae51", "red")
        names(custom_colors) <- c(" ", "GLI3", gene, paste0("GLI3_", gene))
        umap_plot <- ggplot(data = sc_meta_week4, aes(x = umap_1, y = umap_2, color = GLI3_expression)) +
            geom_point(size = 1.5) +
            scale_color_manual(values = custom_colors) +
            custom_theme()
        ggsave(paste0("results/images/Figure_5/week4_GLI3/Zeng_week4_GLI3", gene, ".png"), plot = umap_plot, width = 14, height = 10, dpi = 300)
        gc()
    } else {
        print(paste0(gene, " not in week4"))
    }
}

for (gene in marker_genes) {
    if (gene %in% rownames(seurat_week5)) {
        sc_meta_week5$expression <- seurat_week5@assays$RNA$scale.data[gene, sc_meta_week5$barcode]
        sc_meta_week5$expression <- ifelse(sc_meta_week5$expression > median(sc_meta_week5$expression), gene, " ")
        sc_meta_week5$GLI1_expression <- sapply(c(1:nrow(sc_meta_week5)), function(i) {
            if (sc_meta_week5$GLI1[i] == "GLI1" & sc_meta_week5$expression[i] == gene) {
                return(paste0("GLI1_", gene))
            } else if (sc_meta_week5$GLI1[i] == "GLI1") {
                return("GLI1")
            } else if (sc_meta_week5$expression[i] == gene) {
                return(gene)
            } else {
                return(" ")
            }
        })
        sc_meta_week5$GLI2_expression <- sapply(c(1:nrow(sc_meta_week5)), function(i) {
            if (sc_meta_week5$GLI2[i] == "GLI2" & sc_meta_week5$expression[i] == gene) {
                return(paste0("GLI2_", gene))
            } else if (sc_meta_week5$GLI2[i] == "GLI2") {
                return("GLI2")
            } else if (sc_meta_week5$expression[i] == gene) {
                return(gene)
            } else {
                return(" ")
            }
        })
        sc_meta_week5$GLI3_expression <- sapply(c(1:nrow(sc_meta_week5)), function(i) {
            if (sc_meta_week5$GLI3[i] == "GLI3" & sc_meta_week5$expression[i] == gene) {
                return(paste0("GLI3_", gene))
            } else if (sc_meta_week5$GLI3[i] == "GLI3") {
                return("GLI3")
            } else if (sc_meta_week5$expression[i] == gene) {
                return(gene)
            } else {
                return(" ")
            }
        })

        custom_colors <- c("#e9e9e9", "lightgreen", "#ffae51", "red")
        names(custom_colors) <- c(" ", "GLI1", gene, paste0("GLI1_", gene))
        umap_plot <- ggplot(data = sc_meta_week5, aes(x = umap_1, y = umap_2, color = GLI1_expression)) +
            geom_point(size = 1.5) +
            scale_color_manual(values = custom_colors) +
            custom_theme()
        ggsave(paste0("results/images/Figure_5/week5_GLI1/Zeng_week5_GLI1_", gene, ".png"), plot = umap_plot, width = 14, height = 10, dpi = 300)

        custom_colors <- c("#e9e9e9", "lightgreen", "#ffae51", "red")
        names(custom_colors) <- c(" ", "GLI2", gene, paste0("GLI2_", gene))
        umap_plot <- ggplot(data = sc_meta_week5, aes(x = umap_1, y = umap_2, color = GLI2_expression)) +
            geom_point(size = 1.5) +
            scale_color_manual(values = custom_colors) +
            custom_theme()
        ggsave(paste0("results/images/Figure_5/week5_GLI2/Zeng_week5_GLI2", gene, ".png"), plot = umap_plot, width = 14, height = 10, dpi = 300)

        custom_colors <- c("#e9e9e9", "lightgreen", "#ffae51", "red")
        names(custom_colors) <- c(" ", "GLI3", gene, paste0("GLI3_", gene))
        umap_plot <- ggplot(data = sc_meta_week5, aes(x = umap_1, y = umap_2, color = GLI3_expression)) +
            geom_point(size = 1.5) +
            scale_color_manual(values = custom_colors) +
            custom_theme()
        ggsave(paste0("results/images/Figure_5/week5_GLI3/Zeng_week5_GLI3", gene, ".png"), plot = umap_plot, width = 14, height = 10, dpi = 300)
        gc()
    } else {
        print(paste0(gene, " not in week5"))
    }
}
