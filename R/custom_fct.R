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

plotPCA.DESeqTransform <- function(object, intgroup = "condition",
                                   ntop = 500, returnData = FALSE, pcsToUse = 1:20, batch = NULL) {
  message(paste0("using ntop=", ntop, " top features by variance"))

  # calculate the variance for each gene
  rv <- rowVars(assay(object))

  # select the ntop genes by variance
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]

  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select, ]))

  pca_data_fe <- factoextra::get_pca_var(pca)

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)

  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }

  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])

  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  } else {
    colData(object)[[intgroup]]
  }

  # assembly the data
  pcs <- paste0("PC", pcsToUse)
  d <- pcsToUse %>% lapply(function(PC) {
    return(pca$x[, PC])
  })
  d <- do.call(cbind, d)
  d <- cbind(d, data.frame(group = group, intgroup.df, name = colnames(object)))
  colnames(d)[pcsToUse] <- pcs

  attr(d, "percentVar") <- percentVar[pcsToUse]
  attr(d, "factoextra") <- pca_data_fe
  attr(d, "pca_var") <- prcomp(assay(object)[select, ])
  return(d)
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
  if (hide_legend == TRUE) {
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


gene_converter <- function(gene_vec, from, to) {
  return(mapIds(org.Hs.eg.db,
    key = gene_vec,
    keytype = from, column = to
  ) %>% unname())
}


readcounts <- function(file, sep = ",", header = T, row.names = 1) {
  return(as.matrix(read.table(file, sep = sep, header = header, row.names = row.names)))
}
