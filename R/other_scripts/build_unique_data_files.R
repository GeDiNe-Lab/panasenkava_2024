# Loading packages and functions
library(Matrix)
library(tidyverse)
library(org.Hs.eg.db)

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
F1_rawcounts <- readcounts("/home/jules/Documents/phd/Data/Article_veranika/bulk/counts.csv", sep = ",", header = TRUE)
F1_meta <- read.table("/home/jules/Documents/phd/Data/Article_veranika/bulk/metadata.csv", sep = ",", header = TRUE)

F2_rawcounts <- readcounts("/home/jules/Documents/phd/Data/lab_RNAseq/diff13/diff13_counts.csv", sep = ",", header = TRUE)
F2_meta <- read.table("/home/jules/Documents/phd/Data/lab_RNAseq/diff13/diff13_meta.csv", sep = ",", header = TRUE)

View(F1_meta)
View(F2_meta)

# Formating metadata for diff9 and diff12
F1_meta$CRISPR[which(F1_meta$line == "LON")] <- "no"
F1_meta$sexe <- rep("H", nrow(F1_meta))
F1_meta$day <- rep("day12", nrow(F1_meta))
F1_meta$manip <- rep("veranika", nrow(F1_meta))
F1_meta$line[which(F1_meta$line == "LON")] <- "LON71"
F1_meta$sequencing <- ifelse(F1_meta$line == "LON71", "batch1", "batch2")
F1_meta <- F1_meta[2:ncol(F1_meta)]

# formating metadata for diff13
F2_meta$sequencing <- rep("batch3", nrow(F2_meta))
F2_meta$state <- rep("neuroectoderm", nrow(F2_meta))
F2_meta$cyclo_dose_quant <- rep(0, nrow(F2_meta))
F2_meta$cyclo_dose_qual <- rep("none", nrow(F2_meta))
F2_meta$CRISPR <- rep("no", nrow(F2_meta))

# Merging metadata
meta <- rbind(F1_meta, F2_meta)

# Formating merged metadata
meta$diff <- rep("diff", nrow(meta))
meta[which(meta$sequencing == "batch1"), "diff"] <- "diff9"
meta[which(meta$sequencing == "batch2"), "diff"] <- "diff12"
meta[which(meta$sequencing == "batch3"), "diff"] <- "diff13"
View(meta)

# when creating individual rawcounts matrix we got rid of genes with duplicated names
# especially since they all had 0 counts

# get genes not in common
missing_genes <- setdiff(rownames(F1_rawcounts), rownames(F2_rawcounts)) # gene in F1 not in F2
setdiff(rownames(F2_rawcounts), rownames(F1_rawcounts)) # gene in F2 not in F1 (there is none)

# Adding missing genes to F2_rawcounts (they all have 0 counts)
adding_row <- matrix(0, nrow = length(missing_genes), ncol = ncol(F2_rawcounts))
colnames(adding_row) <- colnames(F2_rawcounts)
rownames(adding_row) <- missing_genes
F2_rawcounts <- rbind(F2_rawcounts, adding_row)

# making sure genes are the same and in same order in both matrix
union(setdiff(rownames(F1_rawcounts), rownames(F2_rawcounts)), setdiff(rownames(F2_rawcounts), rownames(F1_rawcounts)))
F1_rawcounts <- F1_rawcounts[sort(rownames(F1_rawcounts)), ]
F2_rawcounts <- F2_rawcounts[sort(rownames(F2_rawcounts)), ]

# Merging rawcounts
rawcounts <- cbind(F1_rawcounts, F2_rawcounts)

# ordering columns according to metadata
rawcounts <- rawcounts[, meta$sample]

# checking SHH expression
rawcounts["ENSG00000164690", ] # all good

# saving data
write.csv(rawcounts, "/home/jules/Documents/phd/projects/panasenkava_2024/data/rawcounts.csv")
write.csv(meta, "/home/jules/Documents/phd/projects/panasenkava_2024/data/meta.csv", row.names = FALSE)
