# Loading packages and functions
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(DESeq2)

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

F3_rawcounts <- readcounts("/home/jules/Documents/phd/Data/lab_RNAseq/manip4/manip4_counts.csv")
F3_meta <- read.table("/home/jules/Documents/phd/Data/lab_RNAseq/manip4/manip4_metadata.csv", sep = ",", header = T)

F4_rawcounts <- readcounts("/home/jules/Documents/phd/Data/lab_RNAseq/IPSdiff12/IPSdiff12_counts.csv")
F4_rawmeta <- read.table("/home/jules/Documents/phd/Data/lab_RNAseq/IPSdiff12/IPSdiff12_metadata.csv", sep = ",", header = T)

View(F1_meta)
View(F2_meta)
View(F3_meta)
View(F4_rawmeta)
