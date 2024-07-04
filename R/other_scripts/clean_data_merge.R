# Loading packages and functions
library(ggplot2)
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(tibble)

source("/home/jules/Documents/phd/custom_fct.R")

# Loading data
d9_counts <- readcounts("/home/jules/Documents/phd/Data/Article_veranika/bulk/raw_counts_diff9.csv")
d9_meta <- read.table(file = "/home/jules/Documents/phd/Data/Article_veranika/bulk/metadata_diff9.csv", sep = ",", header = T)
d12_counts <- readcounts("/home/jules/Documents/phd/Data/Article_veranika/bulk/raw_counts_diff12.csv")
d12_meta <- read.table(file = "/home/jules/Documents/phd/Data/Article_veranika/bulk/metadata_diff12.csv", sep = ",", header = T)



d12_counts["ENSG00000164690", filter(d12_meta, CRISPR == "control", type %in% c("dorsal", "ventral"))$SAMPLE_NAME]


rownames(d12_counts) <- sub("\\..*", "", rownames(d12_counts))
rownames(d9_counts) <- sub("\\..*", "", rownames(d9_counts))

doublets <- which(table(rownames(d9_counts)) >= 2) %>% names()
min_doub <- doublets %>% sapply(function(doublet) {
    ind <- which(rownames(d9_counts) %in% doublet)
    min_ind <- which(rowSums(d9_counts[ind, ]) == min(rowSums(d9_counts[ind, ])))

    if (length(min_ind) > 1) {
        min_ind <- min_ind[1]
    }
    return(ind[min_ind])
})
d9_counts <- d9_counts[-min_doub, ]

doublets <- which(table(rownames(d12_counts)) >= 2) %>% names()
doublets # noe genes in here

rows <- setdiff(rownames(d9_counts), rownames(d12_counts))
setdiff(rownames(d12_counts), rownames(d9_counts)) # no genes in here

additional_rows <- rows %>% lapply(function(ens) {
    return(rep(0, ncol(d12_counts)))
})
additional_rows <- do.call(rbind, additional_rows)
rownames(additional_rows) <- rows
d12_counts <- rbind(d12_counts, additional_rows)

gene_order <- intersect(rownames(d9_counts), rownames(d12_counts))
counts <- cbind(d9_counts[gene_order, ], d12_counts[gene_order, ])

meta <- data.frame(
    sample = c(d9_meta$new_name, d12_meta$SAMPLE_NAME),
    line = c(rep("LON", nrow(d9_meta)), rep("WTC", nrow(d12_meta))),
    type = c(d9_meta$type, d12_meta$type),
    CRISPR = c(rep("control", nrow(d9_meta)), d12_meta$CRISPR),
    cyclo_dose_quant = c(d9_meta$cyclo_dose, d12_meta$cyclo_dose),
    cyclo_dose_qual = sapply(c(d9_meta$cyclo_dose, d12_meta$cyclo_dose), function(d) {
        if (d >= 0.5) {
            return("high")
        } else if (d >= 0.125) {
            return("low")
        } else {
            return("none")
        }
    })
)
View(meta)

# Saving data
write.csv(counts, "/home/jules/Documents/phd/Data/Article_veranika/bulk/counts.csv")
write.csv(meta, "/home/jules/Documents/phd/Data/Article_veranika/bulk/metadata.csv")
