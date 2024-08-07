#' Perform Principal Component Analysis (PCA) using FactoMineR package
#'
#' This function performs Principal Component Analysis (PCA) using the FactoMineR package.
#' It calculates the PCA scores for individuals and variables, as well as the eigenvalues
#' and percentages of variance explained by each principal component.
#'
#' @param X The data matrix or data frame containing the variables to be analyzed.
#' @param scale.unit Logical value indicating whether the variables should be scaled to have unit variance.
#' @param ncp The number of dimensions to keep in the PCA analysis.
#' @param ind.sup Optional matrix or data frame containing supplementary individuals.
#' @param quanti.sup Optional matrix or data frame containing quantitative supplementary variables.
#' @param quali.sup Optional matrix or data frame containing qualitative supplementary variables.
#' @param row.w Optional vector of weights for the individuals.
#' @param col.w Optional vector of weights for the variables.
#' @param graph Logical value indicating whether to display the PCA graph.
#' @param axes Numeric vector indicating the dimensions to display in the PCA graph.
#'
#' @return A list containing the following components:
#'   \item{ind}{A data frame with the PCA scores for individuals.}
#'   \item{var}{A data frame with the PCA scores for variables.}
#'   \item{eig}{A data frame with the eigenvalues (variances) of the principal components.}
#'   \item{gg.ind}{A data frame with the PCA scores for individuals, including supplementary individuals.}
#'   \item{gg.var}{A data frame with the PCA scores for variables, including supplementary variables.}
#'   \item{gg.prct}{A list with the percentages of variance explained by each principal component.}
#'
#' @examples
#' # Load the data
#' data(iris)
#'
#' # Perform PCA
#' result <- ggPCA(iris[, 1:4])
#'
#' # Access the PCA scores for individuals
#' result$ind
#'
#' # Access the PCA scores for variables
#' result$var
#'
#' # Access the eigenvalues (variances) of the principal components
#' result$eig
#'
#' # Access the percentages of variance explained by each principal component
#' result$gg.prct
#'
#' @import FactoMineR
#' @importFrom dplyr round
#' @importFrom purrr sapply
#' @importFrom stats as.data.frame
#' @importFrom stats rbind
#' @importFrom stats scale
#' @importFrom stats sum
#' @importFrom utils data

ggPCA <- function(X, scale.unit = TRUE, ncp = 5, ind.sup = NULL,
                  quanti.sup = NULL, quali.sup = NULL, row.w = NULL,
                  col.w = NULL, graph = TRUE, axes = c(1, 2)) {
  pca_result <- FactoMineR::PCA(
    X = X, scale.unit = scale.unit, ncp = ncp, ind.sup = ind.sup,
    quanti.sup = quanti.sup, quali.sup = quali.sup, row.w = row.w,
    col.w = col.w, graph = graph, axes = axes
  )

  ind_pca_scores <- pca_result$ind$coord[, 1:ncp]
  ind_pca_scores <- rbind(ind_pca_scores, pca_result$ind.sup$coord[, 1:ncp])
  ind_pca_scores <- as.data.frame(ind_pca_scores)
  colnames(ind_pca_scores) <- c(1:ncp) %>% sapply(function(i) {
    return(paste("PC", i, sep = ""))
  })

  var_pca_scores <- pca_result$var$coord[, 1:ncp]
  var_pca_scores <- rbind(var_pca_scores, pca_result$quali.sup$coord[, 1:ncp])
  var_pca_scores <- rbind(var_pca_scores, pca_result$quanti.sup$coord[, 1:ncp])
  var_pca_scores <- as.data.frame(var_pca_scores)
  colnames(var_pca_scores) <- c(1:ncp) %>% sapply(function(i) {
    return(paste("PC", i, sep = ""))
  })

  # Access the eigenvalues (variances) and percentages of variance explained
  eigenvalues <- pca_result$eig %>% as.data.frame()
  percent_variance <- eigenvalues$eigenvalue / sum(eigenvalues$eigenvalue) * 100
  percent_variance <- percent_variance %>%
    round(2) %>%
    as.list()
  names(percent_variance) <- c(1:length(percent_variance)) %>% sapply(function(i) {
    return(paste("PC", i, sep = ""))
  })
  pca_result$gg.ind <- ind_pca_scores
  pca_result$gg.var <- var_pca_scores
  pca_result$gg.prct <- percent_variance

  return(pca_result)
}

ggAddPCA <- function(meta, pca_result, idcol = NULL) {
  if (is.null(idcol)) {
    stop("Error must inquire an sample name column of the meta table")
  }
  if (is.null(pca_result$gg.ind)) {
    stop("Error on custom ggPCA result format")
  }
  if (sum(colnames(pca_result$gg.ind) %in% colnames(meta)) == length(colnames(pca_result$gg.ind))) {
    pca_result$gg.ind <- pca_result$gg.ind[meta[, idcol], ]
    print("PCA result were here and might have been updated")
    meta[, which(colnames(meta) %in% colnames(pca_result$gg.ind))] <- pca_result$gg.ind
    return(meta)
  } else {
    pca_result$gg.ind <- pca_result$gg.ind[meta[, idcol], ]
    return(cbind(meta, pca_result$gg.ind))
  }
}

#' Custom Theme Function
#'
#' This function creates a custom theme for plots.
#'
#' @param diag_text A logical value indicating whether to display diagonal text on the x-axis labels. Default is \code{FALSE}.
#' @param hide_x_lab A logical value indicating whether to hide the x-axis labels. Default is \code{FALSE}.
#' @param hide_legend A logical value indicating whether to hide the legend. Default is \code{FALSE}.
#'
#' @return A \code{theme} object with customized settings for plot elements.
#'
#' @details This function allows you to create a custom theme for your plots by specifying various settings for plot elements such as title, background, grid lines, axis lines, axis labels, legend text, and legend title.
#'
#' @examples
#' custom_theme(diag_text = TRUE, hide_x_lab = FALSE, hide_legend = FALSE)
#'
#' @export
custom_theme <- function(diag_text = FALSE, hide_legend = FALSE, hide_x_lab = FALSE) {
  if (hide_legend == FALSE) {
    hide_legend <- "none"
  } else {
    hide_legend <- "right"
  }
  if (diag_text == TRUE) {
    text_angle <- 45
  } else {
    text_angle <- 0
  }
  if (hide_x_lab == TRUE) {
    axis_text_x <- element_blank()
  } else {
    axis_text_x <- element_text(size = 12, angle = text_angle, hjust = 1)
  }
  return(theme(
    plot.title = element_text(size = 20),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(color = "#dcdcdc"),
    axis.line = element_line(size = 0.2, colour = "Black"),
    axis.title = element_text(size = 20),
    axis.text.x = axis_text_x,
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.position = hide_legend
  ))
}

#   axis.text.x = element_text(angle = text_angle,hjust = 1,size = 12)

load_library <- function() {
  library(org.Hs.eg.db)
  library(tidyverse)
  library(FactoMineR)
  library(factoextra)
  library(DESeq2)
  library(dplyr)
  library(clusterProfiler)
  library(EnhancedVolcano)
  library(ggrepel)
  library(MASS)
  library(gridExtra)
  library(ggsci)
  library(pheatmap)
  library(data.table)
  library(parallel)
  library(Matrix)
  library(distances)
  library(fastcluster)
  library(parallelDist)
  library(rmarkdown)
  library(umap)
}


normalize_ctmat <- function(ctmat) {
  # Convert the count expression matrix into a cpm matrix
  #
  # Args:
  #   ctmat : count expression matrix
  #
  # Return:
  #   ctmat : count per million expression matrix
  print("Calc libsize...")
  lib_sizes <- ctmat %>% colSums()

  print("Getting Intervals...")
  # Retrieve columns as intervals in non zero values vector
  intervals <- c(1:(length(ctmat@p) - 1)) %>% lapply(function(i) {
    c(ctmat@p[i] + 1, ctmat@p[i + 1])
  })

  print("Normalize Values...")
  # Normalized the values
  ctmat@x <- c(1:length(lib_sizes)) %>%
    lapply(function(i) {
      ctmat@x[unlist(intervals[i])[1]:unlist(intervals[i])[2]] / lib_sizes[i] * 1e6
    }) %>%
    unlist()

  return(ctmat)
}

clean_ctmat <- function(ctmat, gene_thr = 0.1, sample_thr = 2) {
  # Filter out Genes and Cells from ctmat that do not pass given treshold
  #
  # Args:
  #   ctmat : expression matrix
  #   gene_thr : % of cells in which a gene need to be expressed to be kept
  #   sample_thr : standard deviation of n_genes beyond which cells are filtered out
  #
  # Return:
  #   ctmat : cleaned expression matrix

  # Filter out Genes
  genes <- rownames(ctmat)
  genes <- genes[
    rowSums(ctmat > 0) %>% sapply(function(row_count) row_count > ncol(ctmat) * gene_thr) == T
  ]

  ctmat <- ctmat[rownames(ctmat) %in% genes, ]
  gc()

  # Filter out Cells
  valid_cells <- tibble(cell = colnames(ctmat), n_gene = colSums(ctmat > 0))

  n_genes <- list(m = mean(valid_cells$n_gene), sd = sd(valid_cells$n_gene))

  valid_cells <- valid_cells %>%
    filter(
      n_gene > (n_genes$m - (n_genes$sd * sample_thr)),
      n_gene < (n_genes$m + (n_genes$sd * sample_thr))
    )

  ctmat <- ctmat[, valid_cells$cell]

  return(ctmat)
}

load_go_annotation <- function(file_path = "/space/scratch/jgarreau/data/pro_GO.csv") {
  # Load Gene ontology annotation file
  #
  # Args:
  #   file_path : annotation file path
  #
  # Return:
  #   GO_unique_filtered : GO annotation tibble

  GO <- read.delim(file = file_path, sep = ",", stringsAsFactors = TRUE)
  GO_unique <- data.frame(table(GO$GO.ID))
  colnames(GO_unique) <- c("GO", "count")
  GO_unique_filtered <- filter(GO_unique, count >= 20) %>% as.tibble()
  return(GO_unique_filtered)
}


load_manip4 <- function(counts_path = "/home/jules/Documents/phd/Data/VERANIKA/manip4_imports-counts_30082023.csv",
                        meta_path = "/home/jules/Documents/phd/Data/VERANIKA/vera_manip4_metadata.csv", L9C1_2 = F) {
  # Load manip 4 data
  manip4_ds <- read.table(counts_path, header = T, sep = ",", check.names = F) %>% as.matrix()
  meta <- read.table(meta_path, header = T, sep = ",")

  # get rid of L9C1-2
  if (isFALSE(L9C1_2)) {
    print("get rid of L9C1-2...")
    meta <- meta %>% dplyr::filter(new_name != "L9C1_2")
  }
  # select correct samples
  print("get correct samples...")
  manip4_ds <- manip4_ds[, which(colnames(manip4_ds) %in% c("ids", meta$new_name))]

  # remove ".xx" in ensembleids and set ids as rownames
  print("format matrix...")
  rownames(manip4_ds) <- manip4_ds[, "ids"] %>% sapply(function(id) {
    return(unlist(strsplit(id, split = "\\."))[1])
  })

  # remove duplicate ensemble ids
  manip4_ds <- manip4_ds[rownames(manip4_ds) %>% unique(), ]

  # remove ids columns
  manip4_ds <- manip4_ds[, meta$new_name]

  # convert counts string to counts integer
  print("convert counts from strings to integers...")
  manip4_ds <- manip4_ds %>% apply(MARGIN = c(1, 2), FUN = function(input_string) {
    cleaned_string <- gsub("[[:space:]]+", "", input_string)
    return(cleaned_string)
  })
  dim_matrix <- dim(manip4_ds)
  rownames_matrix <- rownames(manip4_ds)
  colnames_matrix <- colnames(manip4_ds)
  manip4_ds <- as.numeric(manip4_ds) %>% matrix(nrow = dim_matrix[1], ncol = dim_matrix[2])
  rownames(manip4_ds) <- rownames_matrix
  colnames(manip4_ds) <- colnames_matrix

  return(list(meta = meta, counts = manip4_ds))
}

gene_converter <- function(gene_vec, from, to) {
  return(mapIds(org.Hs.eg.db,
    key = gene_vec,
    keytype = from, column = to
  ) %>% unname())
}

vstmm <- function(object, blind = TRUE, fitType = "parametric") {
  if (is.null(colnames(object))) {
    colnames(object) <- seq_len(ncol(object))
  }
  if (is.matrix(object)) {
    matrixIn <- TRUE
    object <- DESeqDataSetFromMatrix(
      object, DataFrame(row.names = colnames(object)),
      ~1
    )
  } else {
    matrixIn <- FALSE
  }
  if (is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
    object <- estimateSizeFactors(object)
  }
  if (blind) {
    design(object) <- ~1
  }
  if (blind | is.null(attr(dispersionFunction(object), "fitType"))) {
    object <- estimateDispersionsGeneEst(object, quiet = TRUE)
    object <- estimateDispersionsFit(object,
      quiet = TRUE,
      fitType
    )
  }
  vsd <- getVarianceStabilizedData(object)
  if (matrixIn) {
    return(vsd)
  }
  se <- SummarizedExperiment(
    assays = vsd, colData = colData(object),
    rowRanges = rowRanges(object), metadata = metadata(object)
  )
  DESeqTransform(se)
}


auroc <- function(ranking, links) {
  # Compute the AUROC for the given ranked genes_pair and list of true positive links
  #
  # Args:
  #   ranking: pair of genes ranked by correlation (descending) for network 1
  #   links: list of true positive pair of genes (top 1% of correlation) for network 2
  #
  # Return:
  #   (vector):
  #     n_true : number of true positives
  #     n_full : number of ranked genes pair
  #     auroc : Fraction of the area under the ROC curve
  #     p_value : p_value of the mann-whitney test
  #
  rank_len <- length(ranking)
  links_len <- length(links)
  rank <- c(1:rank_len)
  in_set <- c(ranking %chin% links)
  setStatus <- tibble(item = ranking, rank = rank, in_set = in_set)
  rm(ranking)
  rm(links)
  gc()
  result <- wilcox.test(rank ~ in_set, data = setStatus)
  ns <- setStatus %>%
    group_by(in_set) %>%
    summarize(n = n()) %>%
    arrange(in_set)
  rm(setStatus)
  gc()

  # auc
  auc <- result$statistic / (as.double(ns$n[1]) * as.double(ns$n[2]))
  return(auc)
  # return(c(n_true = links_len, n_full = rank_len, auroc = auc, pvalue = result$p.value))
}

readcounts <- function(file, sep = ",", header = T, row.names = 1) {
  return(as.matrix(read.table(file, sep = sep, header = header, row.names = row.names)))
}

compute_perc_mat <- function(coex_mat) {
  # Compute the percentile matrix for the given coexpression matrix
  #
  # Args:
  #   coex_mat : the coexpression matrix
  #
  # Return:
  #   perc_mat : the percentile matrix

  coex_mat_copy <- coex_mat
  # Set lower triangle to NA so we don't rank pairs twice
  coex_mat_copy[lower.tri(coex_mat_copy, diag = T)] <- NA
  coexVec <- coex_mat_copy %>% as.vector()
  perc_mat <- coexVec %>%
    percent_rank() %>%
    matrix(ncol = ncol(coex_mat))
  rownames(perc_mat) <- rownames(coex_mat)
  colnames(perc_mat) <- colnames(coex_mat)
  # Lower triangle is then set to 0, diagonal to 1
  diag(perc_mat) <- 1
  perc_mat[lower.tri(perc_mat, diag = F)] <- 0
  perc_mat <- perc_mat %>% as("dgCMatrix")
  gc()

  return(perc_mat)
}

compute_lnks_mat <- function(perc_mat, thresh = 0.99) {
  # Compute the links matrix for the given percentile matrix
  #
  # Args:
  #   perc_mat : the percentile matrix
  #   thresh : percentile threshold
  #
  # Return:
  #   lnks_mat : the links matrix
  lnks_mat <- perc_mat

  # Value higher than the treshold (top 0.1% best correlation rank) are links
  lnks_mat@x <- replace(lnks_mat@x, lnks_mat@x >= thresh, 1)

  lnks_mat@x <- replace(lnks_mat@x, lnks_mat@x < thresh, 0)

  # Going back and forth so the 0 values or getting "pushed out" of the values vector
  lnks_mat <- lnks_mat %>%
    as.matrix() %>%
    as("dgCMatrix")
  gc()
  return(lnks_mat)
}


#' get_outlier_genes Function
#'
#' This function identifies outlier genes based on principal component analysis (PCA) of gene expression counts.
#'
#' @param counts A matrix or data frame containing gene expression counts.
#' @param threshold The threshold value for identifying outlier genes. Default is 10.
#'
#' @return A list containing the following elements:
#'   - plot: A ggplot object showing the PCA plot with outlier genes highlighted.
#'   - genes: A character vector containing the names of the identified outlier genes.
#'
#' @import ggplot2
#' @import ggfortify
#' @importFrom scales custom_theme
#'
#' @examples
#' counts <- read.csv("gene_counts.csv")
#' outlier_genes <- get_outlier_genes(counts, threshold = 15)
#'
#' @export
get_outlier_genes <- function(counts, threshold = 10) {
  gene_pca <- ggPCA(counts, ncp = 5, graph = F, scale.unit = T)
  outlier1 <- (gene_pca$gg.ind$PC1 < threshold * -sd(gene_pca$gg.ind$PC1) | gene_pca$gg.ind$PC1 > threshold * sd(gene_pca$gg.ind$PC1))
  outlier2 <- (gene_pca$gg.ind$PC2 < threshold * -sd(gene_pca$gg.ind$PC2) | gene_pca$gg.ind$PC2 > threshold * sd(gene_pca$gg.ind$PC2))
  outlier <- (outlier1 | outlier2)

  geneplot <- ggplot(data = data.frame(
    PC1 = gene_pca$gg.ind$PC1,
    PC2 = gene_pca$gg.ind$PC2,
    outlier = outlier
  ), aes(x = PC1, y = PC2, color = outlier)) +
    geom_point(size = 1) +
    custom_theme()

  return(list(plot = geneplot, genes = rownames(counts)[outlier]))
}
