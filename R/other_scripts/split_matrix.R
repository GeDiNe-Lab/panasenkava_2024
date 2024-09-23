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
# Loading data
rawcounts <- readcounts("data/rawcounts.csv", sep = ",", header = TRUE)
rawmeta <- read.table("data/meta.csv", sep = ",", header = TRUE)
rawmeta %>% View()
cin <- filter(rawmeta, diff == "diff13")
cinLON71 <- filter(rawmeta, diff == "diff13", line == "LON71")
cinWTC <- filter(rawmeta, diff == "diff13", line == "WTC")
cinLON80 <- filter(rawmeta, diff == "diff13", line == "LON80")
cycloLON <- filter(rawmeta, diff == "diff9")
CRISPR_control <- filter(rawmeta, diff == "diff12", CRISPR == "control")
CRISPR_hetero <- filter(rawmeta, diff == "diff12", CRISPR == "hetero")
CRISPR_homo <- filter(rawmeta, diff == "diff12", CRISPR == "homo")
CRISPR <- filter(rawmeta, diff == "diff12")


write.csv(rawcounts[, cin$sample], "data/kinetic.csv")
write.csv(rawcounts[, cinLON71$sample], "data/kinetic_LON71.csv")
write.csv(rawcounts[, cinWTC$sample], "data/kinetic_WTC.csv")
write.csv(rawcounts[, cinLON80$sample], "data/kinetic_LON80.csv")
write.csv(rawcounts[, cycloLON$sample], "data/cyclo_gradient.csv")
write.csv(rawcounts[, CRISPR$sample], "data/CRISPR.csv")
write.csv(rawcounts[, CRISPR_control$sample], "data/CRISPR_control.csv")
write.csv(rawcounts[, CRISPR_hetero$sample], "data/CRISPR_hetero.csv")
write.csv(rawcounts[, CRISPR_homo$sample], "data/CRISPR_homo.csv")
