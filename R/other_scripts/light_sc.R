# Loading packages and functions
library(Matrix)
library(tidyverse)

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

sc_meta <- read.csv("scdata/week345_wholebody_metadata.csv") %>%
    dplyr::select(c("week_stage", "celltype_region_num", "celltype_region", "barcode"))

filtered_ids <- which(sc_meta$celltype_region_num %in% c(1:8))

sc_counts <- readMM("scdata/week345_wholebody.mtx") %>% t()
sc_counts_f <- sc_counts[, filtered_ids]

filtered_genes <- which(rowSums(sc_counts_f) > 0)
filtered_genes

read.table("scdata/indices_week345_wholebody.csv", header = TRUE, sep = ",")[filtered_ids, ] %>%
    write.csv("scdata/light_scdata/indices_week345_wholebody.csv", row.names = FALSE)
read.table("scdata/genes_week345_wholebody.tsv", header = TRUE, sep = "\t")[filtered_genes, ] %>%
    write.csv("scdata/light_scdata/genes_week345_wholebody.tsv", row.names = FALSE)

writeMM(t(sc_counts_f[filtered_genes, ]), "scdata/light_scdata/week345_wholebody.mtx")
