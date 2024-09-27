# Loading packages and functions
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(ggpubr)
library(DESeq2)
library(gridExtra)

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

Het_Con_DE <- read.csv("results/tables/Figure_4/DE_vAN_het_vs_control.csv")
Hom_Con_DE <- read.csv("results/tables/Figure_4/DE_vAN_homo_vs_control.csv")
Hom_Het_DE <- read.csv("results/tables/Figure_4/DE_vAN_homo_vs_het.csv")

load("results/tables/Figure_3/cyclo_genes_df.RData")
cyclopamin <- cyclo_genes_df
# Loading the VennDiagram package
library(VennDiagram)

# CRISPR DEGs
Het_Con <- filter(Het_Con_DE, abs(log2FoldChange) >= 1, padj < 0.05, !is.na(gene))$gene
Hom_Con <- filter(Hom_Con_DE, abs(log2FoldChange) >= 1, padj < 0.05, !is.na(gene))$gene
Hom_Het <- filter(Hom_Het_DE, abs(log2FoldChange) >= 1, padj < 0.05, !is.na(gene))$gene

# Cyclopamin DEGs
High_No <- filter(cyclopamin, HvsN_thr == TRUE, !is.na(genes))$genes
High_Low <- filter(cyclopamin, HvsL_thr == TRUE, !is.na(genes))$genes
Low_No <- filter(cyclopamin, LvsN_thr == TRUE, !is.na(genes))$genes

# 4way venn diag
# Het_Con vs Hom_Con vs High_No
venn.diagram(
    x = list(Hom_Con, Het_Con, High_No, Low_No),
    disable.logging = TRUE,
    category.names = c("Hom_Con", "Het_Con", "High_No", "Low_No"),
    main = "Hom_Con   Het_Con    High_No",
    filename = "results/images/Figure_4/venn_CRISPR_cyclo_general.png"
)

# Generate the Venn diagram and capture the output
venn_output <- VennDiagram::get.venn.partitions(
    list(
        "Hom_Con" = Hom_Con,
        "Het_Con" = Het_Con,
        "High_No" = High_No,
        "Low_No" = Low_No
    )
) %>% dplyr::select(..set.., ..values..)
colnames(venn_output) <- c("set", "genes")
venn_df <- data.frame(genes = Reduce(union, list(Hom_Con, Het_Con, High_No, Low_No)))

for (i in seq_len(nrow(venn_output))) {
    set <- venn_output$set[i]
    set_genes <- venn_output$genes[[i]]

    # Create a logical column (TRUE/FALSE) indicating if the gene is in the set
    venn_df[[set]] <- venn_df$genes %in% set_genes
}
venn_df %>% View()
View(dplyr::select(venn_output, ..set.., ..count..))
# Extract the intersections


# Write the dataframe to a CSV file
write.csv(venn_df, "results/tables/Figure_4/venn_CRISPR_cyclo_general.csv", row.names = FALSE)

# 2way venn diag

# Het_Con vs High_No
venn.diagram(
    x = list(Het_Con, High_No),
    disable.logging = TRUE,
    category.names = c("Het_Con", "High_No"),
    main = "Het_Con    High_No",
    filename = "results/images/Figure_4/venn_Het_Con_vs_High_No.png"
)
Het_Con_vs_High_No <- data.frame(
    genes = union(Het_Con, High_No)
)
Het_Con_vs_High_No$state <- Het_Con_vs_High_No$genes %>% sapply(function(gene) {
    if (gene %in% Het_Con & gene %in% High_No) {
        return("common")
    } else if (gene %in% Het_Con) {
        return("Het_Con")
    } else {
        return("High_No")
    }
})
write.csv(Het_Con_vs_High_No, "results/tables/Figure_4/venn_Het_Con_vs_High_No.csv")


# Het_Con vs High_Low
venn.diagram(
    x = list(Het_Con, High_Low),
    category.names = c("Het_Con", "High_Low"),
    disable.logging = TRUE,
    main = "Het_Con    High_Low",
    filename = "results/images/Figure_4/venn_Het_Con_vs_High_Low.png"
)
Het_Con_vs_High_Low <- data.frame(
    genes = union(Het_Con, High_Low)
)
Het_Con_vs_High_Low$state <- Het_Con_vs_High_Low$genes %>% sapply(function(gene) {
    if (gene %in% Het_Con & gene %in% High_Low) {
        return("common")
    } else if (gene %in% Het_Con) {
        return("Het_Con")
    } else {
        return("High_Low")
    }
})
write.csv(Het_Con_vs_High_Low, "results/tables/Figure_4/venn_Het_Con_vs_High_Low.csv")


# Het_Con vs Low_No
venn.diagram(
    x = list(Het_Con, Low_No),
    category.names = c("Het_Con", "Low_No"),
    disable.logging = TRUE,
    main = "Het_Con    Low_No",
    filename = "results/images/Figure_4/venn_Het_Con_vs_Low_No.png"
)
Het_Con_vs_Low_No <- data.frame(
    genes = union(Het_Con, Low_No)
)
Het_Con_vs_Low_No$state <- Het_Con_vs_Low_No$genes %>% sapply(function(gene) {
    if (gene %in% Het_Con & gene %in% Low_No) {
        return("common")
    } else if (gene %in% Het_Con) {
        return("Het_Con")
    } else {
        return("Low_No")
    }
})
write.csv(Het_Con_vs_Low_No, "results/tables/Figure_4/venn_Het_Con_vs_Low_No.csv")


# Hom_Con vs High_No
venn.diagram(
    x = list(Hom_Con, High_No),
    category.names = c("Hom_Con", "High_No"),
    disable.logging = TRUE,
    main = "Hom_Con    High_No",
    filename = "results/images/Figure_4/venn_Hom_Con_vs_High_No.png"
)
Hom_Con_vs_High_No <- data.frame(
    genes = union(Hom_Con, High_No)
)
Hom_Con_vs_High_No$state <- Hom_Con_vs_High_No$genes %>% sapply(function(gene) {
    if (gene %in% Hom_Con & gene %in% High_No) {
        return("common")
    } else if (gene %in% Hom_Con) {
        return("Hom_Con")
    } else {
        return("High_No")
    }
})
write.csv(Hom_Con_vs_High_No, "results/tables/Figure_4/venn_Hom_Con_vs_High_No.csv")


# Hom_Con vs High_Low
venn.diagram(
    x = list(Hom_Con, High_Low),
    category.names = c("Hom_Con", "High_Low"),
    disable.logging = TRUE,
    main = "Hom_Con    High_Low",
    filename = "results/images/Figure_4/venn_Hom_Con_vs_High_Low.png"
)
Hom_Con_vs_High_Low <- data.frame(
    genes = union(Hom_Con, High_Low)
)
Hom_Con_vs_High_Low$state <- Hom_Con_vs_High_Low$genes %>% sapply(function(gene) {
    if (gene %in% Hom_Con & gene %in% High_Low) {
        return("common")
    } else if (gene %in% Hom_Con) {
        return("Hom_Con")
    } else {
        return("High_Low")
    }
})
write.csv(Hom_Con_vs_High_Low, "results/tables/Figure_4/venn_Hom_Con_vs_High_Low.csv")


# Hom_Con vs Low_No
venn.diagram(
    x = list(Hom_Con, Low_No),
    category.names = c("Hom_Con", "Low_No"),
    disable.logging = TRUE,
    main = "Hom_Con    Low_No",
    filename = "results/images/Figure_4/venn_Hom_Con_vs_Low_No.png"
)
Hom_Con_vs_Low_No <- data.frame(
    genes = union(Hom_Con, Low_No)
)
Hom_Con_vs_Low_No$state <- Hom_Con_vs_Low_No$genes %>% sapply(function(gene) {
    if (gene %in% Hom_Con & gene %in% Low_No) {
        return("common")
    } else if (gene %in% Hom_Con) {
        return("Hom_Con")
    } else {
        return("Low_No")
    }
})
write.csv(Hom_Con_vs_Low_No, "results/tables/Figure_4/venn_Hom_Con_vs_Low_No.csv")


# Hom_Het vs High_No
venn.diagram(
    x = list(Hom_Het, High_No),
    category.names = c("Hom_Het", "High_No"),
    disable.logging = TRUE,
    main = "Hom_Het    High_No",
    filename = "results/images/Figure_4/venn_Hom_Het_vs_High_No.png"
)
Hom_Het_vs_High_No <- data.frame(
    genes = union(Hom_Het, High_No)
)
Hom_Het_vs_High_No$state <- Hom_Het_vs_High_No$genes %>% sapply(function(gene) {
    if (gene %in% Hom_Het & gene %in% High_No) {
        return("common")
    } else if (gene %in% Hom_Het) {
        return("Hom_Het")
    } else {
        return("High_No")
    }
})
write.csv(Hom_Het_vs_High_No, "results/tables/Figure_4/venn_Hom_Het_vs_High_No.csv")


# Hom_Het vs High_Low
venn.diagram(
    x = list(Hom_Het, High_Low),
    category.names = c("Hom_Het", "High_Low"),
    disable.logging = TRUE,
    main = "Hom_Het    High_Low",
    filename = "results/images/Figure_4/venn_Hom_Het_vs_High_Low.png"
)
Hom_Het_vs_High_Low <- data.frame(
    genes = union(Hom_Het, High_Low)
)
Hom_Het_vs_High_Low$state <- Hom_Het_vs_High_Low$genes %>% sapply(function(gene) {
    if (gene %in% Hom_Het & gene %in% High_Low) {
        return("common")
    } else if (gene %in% Hom_Het) {
        return("Hom_Het")
    } else {
        return("High_Low")
    }
})
write.csv(Hom_Het_vs_High_Low, "results/tables/Figure_4/venn_Hom_Het_vs_High_Low.csv")

# Hom_Het vs Low_No
venn.diagram(
    x = list(Hom_Het, Low_No),
    category.names = c("Hom_Het", "Low_No"),
    disable.logging = TRUE,
    main = "Hom_Het    Low_No",
    filename = "results/images/Figure_4/venn_Hom_Het_vs_Low_No.png"
)
Hom_Het_vs_Low_No <- data.frame(
    genes = union(Hom_Het, Low_No)
)
Hom_Het_vs_Low_No$state <- Hom_Het_vs_Low_No$genes %>% sapply(function(gene) {
    if (gene %in% Hom_Het & gene %in% Low_No) {
        return("common")
    } else if (gene %in% Hom_Het) {
        return("Hom_Het")
    } else {
        return("Low_No")
    }
})
write.csv(Hom_Het_vs_Low_No, "results/tables/Figure_4/venn_Hom_Het_vs_Low_No.csv")

meta_CRISPR <- filter(rawmeta, diff == "diff12", type == "ventral")
meta_cyclo <- filter(rawmeta, diff == "diff9", type != "dorsal", sample != "L9C1_2")

SHH_pathway_genes <- read.table("data/SHH_pathway_genes.tab", sep = "\t", header = TRUE)

markers <- SHH_pathway_genes$Symbol %>% gene_converter("SYMBOL", "ENSEMBL")
# Make sure to keep only markers in the matrix
markers <- intersect(markers, rownames(rawcounts))

# filtering out lowly expressed genes and keeping only seleted samples
retained_row_CRISPR <- rawcounts[rownames(rawcounts) %in% markers, meta_CRISPR$sample]
counts_CRISPR <- rawcounts[, meta_CRISPR$sample][rowSums(rawcounts[, meta_CRISPR$sample]) >= 25, ]

retained_row_cyclo <- rawcounts[rownames(rawcounts) %in% markers, meta_cyclo$sample]
counts_cyclo <- rawcounts[, meta_cyclo$sample][rowSums(rawcounts[, meta_cyclo$sample]) >= 25, ]

# putting back potentially filtered out markers (posterior markers for example as they should not be expressed)
if (length(which(!rownames(retained_row_CRISPR) %in% rownames(counts_CRISPR))) == 1) {
    counts_CRISPR <- rbind(counts_CRISPR, retained_row_CRISPR[which(!rownames(retained_row_CRISPR) %in% rownames(counts_CRISPR)), ])
    rownames(counts_CRISPR)[nrow(counts_CRISPR)] <- rownames(retained_row_CRISPR)[which(!rownames(retained_row_CRISPR) %in% rownames(counts_CRISPR))]
} else if (length(which(!rownames(retained_row_CRISPR) %in% rownames(counts_CRISPR))) > 1) {
    counts_CRISPR <- rbind(counts_CRISPR, retained_row_CRISPR[which(!rownames(retained_row_CRISPR) %in% rownames(counts_CRISPR)), ])
}

if (length(which(!rownames(retained_row_cyclo) %in% rownames(counts_cyclo))) == 1) {
    counts_cyclo <- rbind(counts_cyclo, retained_row_cyclo[which(!rownames(retained_row_cyclo) %in% rownames(counts_cyclo)), ])
    rownames(counts_cyclo)[nrow(counts_cyclo)] <- rownames(retained_row_cyclo)[which(!rownames(retained_row_cyclo) %in% rownames(counts_cyclo))]
} else if (length(which(!rownames(retained_row_cyclo) %in% rownames(counts_cyclo))) > 1) {
    counts_cyclo <- rbind(counts_cyclo, retained_row_cyclo[which(!rownames(retained_row_cyclo) %in% rownames(counts_cyclo)), ])
}


dds_cyclo <- DESeqDataSetFromMatrix(
    countData = counts_cyclo,
    colData = meta_cyclo,
    design = ~cyclo_dose_qual
)
vsd_cyclo <- vst(dds_cyclo, blind = FALSE)
dds_CRISPR <- DESeqDataSetFromMatrix(
    countData = counts_CRISPR,
    colData = meta_CRISPR,
    design = ~CRISPR
)
vsd_CRISPR <- vst(dds_CRISPR, blind = FALSE)

# Looking at SHH co-expressed genes expression in ventral CRISPR samples
scale_cyclo <- assay(vsd_cyclo) - min(assay(vsd_cyclo))
rownames(scale_cyclo) <- gene_converter(rownames(scale_cyclo), "ENSEMBL", "SYMBOL")
scale_CRISPR <- assay(vsd_CRISPR) - min(assay(vsd_CRISPR))
rownames(scale_CRISPR) <- gene_converter(rownames(scale_CRISPR), "ENSEMBL", "SYMBOL")

meta_bp_CRISPR <- meta_CRISPR
meta_bp_CRISPR$CRISPR_type <- paste(meta_bp_CRISPR$CRISPR, meta_bp_CRISPR$type, sep = "_")

meta_bp_cyclo <- meta_cyclo

for (gene in SHH_pathway_genes$Symbol) {
    meta_bp_CRISPR$gene <- scale_CRISPR[gene, ]
    red_df_CRISPR <- dplyr::select(meta_bp_CRISPR, c("CRISPR", "gene"))
    stat_df_CRISPR <- data.frame(control = filter(red_df_CRISPR, CRISPR == "control")$gene, hetero = filter(red_df_CRISPR, CRISPR == "hetero")$gene, homo = filter(red_df_CRISPR, CRISPR == "homo")$gene)

    grp_df_CRISPR <- meta_bp_CRISPR %>%
        group_by(CRISPR) %>%
        summarize(
            gene_mean = mean(gene),
            gene_sd = sd(gene)
        )

    grp_df_CRISPR$CRISPR <- factor(c("+/+", "+/-", "-/-"), levels = c("+/+", "+/-", "-/-"))
    c_het_p <- wilcox.test(stat_df_CRISPR$control, stat_df_CRISPR$hetero)$p.value %>% round(4)
    c_hom_p <- wilcox.test(stat_df_CRISPR$control, stat_df_CRISPR$homo)$p.value %>% round(4)
    het_hom_p <- wilcox.test(stat_df_CRISPR$hetero, stat_df_CRISPR$homo)$p.value %>% round(4)
    signif <- c(c_het_p, c_hom_p, het_hom_p) %>% sapply(function(pval) {
        if (is.nan(pval)) {
            return("***")
        } else if (pval < 0.001) {
            return("***")
        } else if (pval < 0.01) {
            return("**")
        } else if (pval < 0.05) {
            return("*")
        } else {
            return(" ")
        }
    })
    pval_data <- data.frame(
        group1 = c("+/+", "+/+", "+/-"),
        group2 = c("+/-", "-/-", "-/-"),
        p_value = paste(c(c_het_p, c_hom_p, het_hom_p), signif, sep = " "),
        y_position = c(max(scale_CRISPR) + max(scale_CRISPR) * 0.1, max(scale_CRISPR) + max(scale_CRISPR) * 0.2, max(scale_CRISPR) + max(scale_CRISPR) * 0.3) # Adjust these based on your plot's scale
    )
    bp_CRISPR <- ggbarplot(grp_df_CRISPR,
        x = "CRISPR", y = "gene_mean",
        fill = "CRISPR", palette = c("#80AD3C", "#b9e27b", "#dfe981"),
        add = "gene_sd",
        width = 0.9,
        ggtheme = custom_theme(diag_text = TRUE, hide_legend = TRUE)
    ) +
        ylim(-0.2, max(scale_CRISPR) + max(scale_CRISPR) * 0.3) +
        geom_errorbar(data = grp_df_CRISPR, aes(ymin = gene_mean - gene_sd, ymax = gene_mean + gene_sd), width = 0.2) +
        stat_pvalue_manual(
            pval_data,
            y.position = "y_position",
            label = "p_value"
        )


    meta_bp_cyclo$gene <- scale_cyclo[gene, ]
    red_df_cyclo <- dplyr::select(meta_bp_cyclo, c("cyclo_dose_qual", "gene"))

    grp_df_cyclo <- meta_bp_cyclo %>%
        group_by(cyclo_dose_qual) %>%
        summarize(
            gene_mean = mean(gene),
            gene_sd = sd(gene)
        )
    grp_df_cyclo <- grp_df_cyclo[c(3, 2, 1), ] # making sure everything is in the right order
    grp_df_cyclo$cyclo_dose_qual <- factor(c("none", "low", "high"), levels = c("none", "low", "high"))
    no_low_p <- wilcox.test(filter(red_df_cyclo, cyclo_dose_qual == "none")$gene, filter(red_df_cyclo, cyclo_dose_qual == "low")$gene)$p.value %>% round(4)
    no_high_p <- wilcox.test(filter(red_df_cyclo, cyclo_dose_qual == "none")$gene, filter(red_df_cyclo, cyclo_dose_qual == "high")$gene)$p.value %>% round(4)
    low_high_p <- wilcox.test(filter(red_df_cyclo, cyclo_dose_qual == "low")$gene, filter(red_df_cyclo, cyclo_dose_qual == "high")$gene)$p.value %>% round(4)
    signif <- c(no_low_p, no_high_p, low_high_p) %>% sapply(function(pval) {
        if (is.nan(pval)) {
            return("")
        } else if (pval < 0.001) {
            return("***")
        } else if (pval < 0.01) {
            return("**")
        } else if (pval < 0.05) {
            return("*")
        } else {
            return("")
        }
    })
    pval_data <- data.frame(
        group1 = c("none", "none", "low"),
        group2 = c("low", "high", "high"),
        p_value = paste(c(no_low_p, no_high_p, low_high_p), signif, sep = " "),
        y_position = c(max(scale_cyclo) + max(scale_cyclo) * 0.1, max(scale_cyclo) + max(scale_cyclo) * 0.2, max(scale_cyclo) + max(scale_cyclo) * 0.3) # Adjust these based on your plot's scale
    )
    bp_cyclo <- ggbarplot(grp_df_cyclo,
        x = "cyclo_dose_qual", y = "gene_mean",
        fill = "cyclo_dose_qual", palette = c("#80AD3C", "#f6a563", "#f86110"),
        add = "gene_sd",
        width = 0.9,
        ggtheme = custom_theme(diag_text = TRUE, hide_legend = TRUE)
    ) +
        ylim(-0.2, max(scale_cyclo) + max(scale_cyclo) * 0.3) +
        geom_errorbar(data = grp_df_cyclo, aes(ymin = gene_mean - gene_sd, ymax = gene_mean + gene_sd), width = 0.2) +
        stat_pvalue_manual(
            pval_data,
            y.position = "y_position",
            label = "p_value"
        )
    title <- textGrob(paste0("Scaled normalized expression of ", gene), gp = gpar(fontsize = 16, fontface = "bold"))

    # Arrange the two plots and add the title
    finalplot <- grid.arrange(
        title, # Add the title
        arrangeGrob(bp_CRISPR, bp_cyclo, ncol = 2), # Add the plots side by side
        nrow = 2, # The layout will have two rows, one for title and one for the plots
        heights = c(0.1, 1) # Relative heights (title takes 10% of the height)
    )
    ggsave(filename = paste0("results/images/Figure_4/vAN_SHH_pathway_genes_barplot/barplot_", gene, ".png"), plot = finalplot, units = "px", width = 1800, height = 1400, dpi = 250)
}
SHH_pathway_genes$Symbol
