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

fig_genelist_1 <- read.table("results/tables/Figure_1/fig_genelist_1.csv", sep = ",", header = TRUE)
fig_genelist_2 <- read.table("results/tables/Figure_2/fig_genelist_2.csv", sep = ",", header = TRUE)
fig_genelist_3 <- read.table("results/tables/Figure_2/fig_genelist_3.csv", sep = ",", header = TRUE)
fig_genelist_4 <- read.table("results/tables/Figure_2/fig_genelist_4.csv", sep = ",", header = TRUE)
fig_genelist_5 <- read.table("results/tables/Figure_4/fig_genelist_5.csv", sep = ",", header = TRUE)
fig_genelist_1
fig_genelist_2
fig_genelist_3
fig_genelist_4
fig_genelist_5
# List of dataframes
df_list <- list(fig_genelist_1, fig_genelist_2, fig_genelist_3, fig_genelist_4, fig_genelist_5)

# Merge using full_join iteratively
merged_df <- Reduce(function(x, y) full_join(x, y, by = "X"), df_list)
merged_df$genes <- merged_df$X %>% gene_converter("ENSEMBL", "SYMBOL")
# Print result
View(merged_df)

write.csv(merged_df, "results/tables/Figure_genelist.csv", row.names = FALSE)
