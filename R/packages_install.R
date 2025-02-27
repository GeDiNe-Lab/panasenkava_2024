sup_install <- function() {
    BiocManager::install(c("ComplexHeatmap", "DEGreport", "DESeq2", "org.Hs.eg.db"))
    remotes::install_github("wilkelab/cowplot")
    install.packages("colorspace", repos = "http://R-Forge.R-project.org")
    remotes::install_github("clauswilke/colorblindr")
}

sup_install()
