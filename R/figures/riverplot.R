# Loading packages and functions
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(DESeq2)
library(Seurat)
library(reshape2)
library(networkD3)

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

figure_genes <- read.csv("results/tables/Figure_genelist.csv", header = TRUE)
# Example gene dataframe
df <- figure_genes
df <- dplyr::select(df, -c(genes, X))
colnames(df) <- c("F1E", "F2C", "F2D", "DEGpattern", "WGCNA", "single-cell", "F5")
df$DEGpattern <- as.character(df$DEGpattern)


layer1 <- c("single-cell") %>%
    lapply(function(col) {
        lk <- c("DEGpattern") %>% lapply(function(other_col) {
            return(table(paste(paste(paste0(col, " "), df[, col], sep = "-"), paste(other_col, df[, other_col], sep = "-"), sep = "_")))
        })
        return(lk)
    }) %>%
    unlist()
links1 <- names(layer1) %>%
    lapply(function(link) {
        return(data.frame(source = str_split(link, "_")[[1]][1], target = str_split(link, "_")[[1]][2], value = layer1[[link]]))
    }) %>%
    bind_rows()

layer2 <- c("DEGpattern") %>%
    lapply(function(col) {
        lk <- c("WGCNA") %>% lapply(function(other_col) {
            return(table(paste(paste(col, df[, col], sep = "-"), paste(other_col, df[, other_col], sep = "-"), sep = "_")))
        })
        return(lk)
    }) %>%
    unlist()
links2 <- names(layer2) %>%
    lapply(function(link) {
        return(data.frame(source = str_split(link, "_")[[1]][1], target = str_split(link, "_")[[1]][2], value = layer2[[link]]))
    }) %>%
    bind_rows()

layer3 <- c("WGCNA") %>%
    lapply(function(col) {
        lk <- c("single-cell") %>% lapply(function(other_col) {
            return(table(paste(paste(col, df[, col], sep = "-"), paste(other_col, df[, other_col], sep = "-"), sep = "_")))
        })
        return(lk)
    }) %>%
    unlist()
links3 <- names(layer3) %>%
    lapply(function(link) {
        return(data.frame(source = str_split(link, "_")[[1]][1], target = str_split(link, "_")[[1]][2], value = layer3[[link]]))
    }) %>%
    bind_rows()

links <- bind_rows(links1, links2, links3)
# links <- filter(links, !source %in% c("DEGpattern-NA", "WGCNA-NA"), !target %in% c("DEGpattern-NA", "WGCNA-NA"))
# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
    name = c(
        as.character(links$source),
        as.character(links$target)
    ) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name) - 1
links$IDtarget <- match(links$target, nodes$name) - 1

# Make the Network
p <- sankeyNetwork(
    Links = links, Nodes = nodes,
    Source = "IDsource", Target = "IDtarget",
    Value = "value", NodeID = "name",
    sinksRight = FALSE
)
p
