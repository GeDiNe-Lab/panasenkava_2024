kinetic_lineplots <- function(data) {
  plot <- ggplot(data, aes(x = day, y = mean_expression, color = gene, group = gene, label = gene)) +
    geom_point(size = 2) +
    geom_line(size = 2) +
    geom_errorbar(
      aes(
        ymin = mean_expression - sd_expression,
        ymax = mean_expression + sd_expression
      ),
      width = 0, # Width of the horizontal bar on the error bar
      size = 1 # Thickness of the error bars
    ) +
    geom_text_repel(
      data = filter(data, day == "day12"),
      aes(x = 7.75),
      nudge_x = 0,
      nudge_y = 0.1,
      direction = "y",
      size = 10,
      segment.color = NA,
      max.overlaps = Inf
    ) +
    scale_x_discrete(breaks = c("day0", "day02", "day04", "day06", "day08", "day10", "day12"), expand = c(0, 0)) +
    ylim(0, 5) +
    ylab("Scaled normalized expression") +
    custom_theme() +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 30),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 20),
      plot.margin = margin(t = 10, r = 30, b = 10, l = 10)
    )
  return(plot)
}
png_save <- function(plot, path, width, height, res = 250) {
  #' Save a plot as a PNG file
  #'
  #' This function saves a given plot to a specified file path in PNG format.
  #'
  #' @param plot The plot object to be saved.
  #' @param path A character string specifying the file path where the PNG file will be saved.
  #' @param width The width of the PNG file in pixels.
  #' @param height The height of the PNG file in pixels.
  #' @param res The resolution of the PNG file in dots per inch (DPI). Default is 250.
  #'
  #' @return None. The function is called for its side effect of saving the plot to a file.
  #'
  #' @examples
  #' \dontrun{
  #'   library(ggplot2)
  #'   p <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point()
  #'   png_save(p, "plot.png", width = 800, height = 600)
  #' }
  png_save <- function(plot, path, width, height, res = 250) {
    png(filename = path, width = width, height = height, res = res)
    draw(plot)
    dev.off()
  }
  png(filename = path, width = width, height = height, res = res)
  draw(plot)
  dev.off()
}

pca_anova <- function(pca_data, metadata, covariates) {
  # Building dataframe with first 5 PC and covariates
  PC_covariate <- cbind(pca_data[, 1:5], metadata %>%
    dplyr::select(all_of(covariates)) %>%
    apply(2, function(x) {
      return(as.numeric(factor(x)) - 1)
    }) %>%
    as.matrix())

  # Computing PC-covariate correlation and ANOVA
  PC_covariate_cor <- cor(PC_covariate[, 1:5], PC_covariate[, 6:ncol(PC_covariate)]) %>% abs()
  PC_covariate_ANOVA <- c(6:ncol(PC_covariate)) %>% lapply(function(i) {
    apply(PC_covariate[, 1:5], 2, function(x) {
      aov(x ~ PC_covariate[, i])
    }) %>% sapply(function(x) {
      summary(x)[[1]]$`Pr(>F)`[1]
    })
  })
  PC_covariate_ANOVA <- Reduce(cbind, PC_covariate_ANOVA)
  colnames(PC_covariate_ANOVA) <- colnames(PC_covariate)[6:ncol(PC_covariate)]

  return(PC_covariate_ANOVA)
}


plot_go_term <- function(genelist, path, range = c(1:20), cut = 40, textsize = 20, imgw = 17, imgh = 10) {
  GO_enrichment <- clusterProfiler::enrichGO(genelist,
    OrgDb = "org.Hs.eg.db",
    keyType = "ENSEMBL",
    ont = "BP"
  )
  GO_results <- GO_enrichment@result
  GO_results$GeneRatio <- sapply(GO_enrichment@result$GeneRatio, function(x) {
    eval(parse(text = x))
  }) %>% unname()
  GO_results$rank <- rank(-GO_results$GeneRatio, ties.method = "first")

  GO_results_f <- GO_results[order(GO_results$GeneRatio, decreasing = TRUE)[range], ]

  GO_results_f$Description <- str_wrap(GO_results_f$Description, width = cut) %>% str_to_upper()
  GO_results_f$Description <- factor(GO_results_f$Description, levels = rev(GO_results_f$Description))

  goplot <- ggplot(GO_results_f, aes(x = GeneRatio, y = reorder(Description, GeneRatio), fill = p.adjust)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Description),
      hjust = 1.01, # Move text inside the bar, adjust as needed
      color = "#000000", # Make the text white for better visibility
      size = 13
    ) + # Adjust size to fit the text inside the bar
    custom_theme() +
    scale_fill_gradient(name = "p-value", low = "#e06663", high = "#327eba") +
    theme(
      axis.title.x = element_text(size = 30), # Adjusts the x-axis title size
      axis.text.x = element_text(size = textsize),
      axis.text.y = element_blank(), # Remove y-axis text
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      legend.text = element_text(size = 20), # Adjusts the legend text size
      legend.title = element_text(size = 30), # Adjusts the legend title size
      legend.key.size = unit(2, "lines")
    )
  ggsave(paste0(path, ".png"), goplot, width = imgw, height = imgh)
  return(GO_results)
}

MyDegPlotCluster <- function(table, time, sign_comp, cluster_i, color = NULL,
                             min_genes = 10,
                             process = FALSE,
                             cluster_column = "cluster",
                             prefix_title = "Group: ") {
  stopifnot(class(table)[1] == "data.frame")

  if (cluster_column %in% colnames(table)) {
    table[["cluster"]] <- table[[cluster_column]]
  }
  if (process) {
    table <- .process(table, time, color)
    print("huu ?")
  }

  if ("cluster" %in% colnames(table)) {
    counts <- table(distinct(table, genes, cluster)[["cluster"]])
    counts <- counts[counts >= min_genes]
    if (length(counts) == 0) {
      stop("No clusters with min_genes > ", min_genes)
    }
    table <- inner_join(table,
      data.frame(
        cluster = as.integer(names(counts)),
        title = paste(
          prefix_title,
          names(counts),
          "- genes:",
          counts
        ),
        stringsAsFactors = FALSE
      ),
      by = "cluster"
    )
  }

  if (is.null(color)) {
    color <- "dummy"
    table[[color]] <- ""
    lines <- FALSE
  }
  table[["line_group"]] <- paste(
    table[["genes"]],
    table[[color]]
  )
  splan <- length(unique(table[[time]])) - 1L
  table$title <- table$title %>% as.factor()
  old <- table$title %>%
    levels() %>%
    stringi::stri_extract_first(., regex = "\\d+") %>%
    as.integer()
  index <- sort(old, index.return = TRUE)[[2]]
  table$title <- factor(table$title, levels = levels(table$title)[index])

  plot <- ggplot(table, aes_string(
    x = time, y = "value"
  )) +
    geom_boxplot(
      alpha = 0,
      outlier.size = 0,
      outlier.shape = NA,
      color = "#80AD3C", # Fixed color for boxplot outline
      fill = "#80AD3C"
    ) +
    geom_jitter(
      alpha = 0.4, size = 1,
      color = "#80AD3C",
      fill = "#80AD3C"
    ) +
    stat_smooth(
      aes_string(
        x = time, y = "value",
        group = color, color = color
      ),
      se = FALSE,
      method = "lm", formula = y ~ poly(x, splan),
      color = "black"
    ) +
    geom_line(aes_string(group = "line_group"), alpha = 0.1, color = "#80AD3C") +
    geom_signif(
      comparisons = sign_comp, # Groups being compared
      map_signif_level = TRUE, # Automatically converts p-values to significance stars
      y_position = rep(2, length(sign_comp)),
      color = "black",
      size = 1,
      step_increase = 0.05,
      textsize = 10
    ) +
    ggtitle(paste0("Cluster ", cluster_i, ", ", length(unique(table$genes)), " genes")) +
    ylab("Z-score of gene abundance") +
    xlab("") +
    ylim(-2.1, 3.2) +
    custom_theme(hide_legend = TRUE) +
    theme(
      plot.title = element_text(size = 35),
      axis.title.y = element_text(size = 30),
      axis.text.x = element_text(size = 30)
    )
  plot
}



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
                                   ntop = 500, returnData = FALSE, pcsToUse = 1:15, batch = NULL) {
  message(paste0("using ntop=", ntop, " top features by variance"))

  # calculate the variance for each gene

  rv <- MatrixGenerics::rowVars(assay(object))

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
  print(pcsToUse)

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


gene_converter <- function(gene_vec, from, to, species = "human") {
  if (species == "human") {
    db <- org.Hs.eg.db
  } else if (species == "mouse") {
    db <- org.Mm.eg.db
  }
  SFTA3_i <- which(gene_vec == "ENSG00000229415")
  sym <- mapIds(db,
    key = gene_vec,
    keytype = from, column = to
  ) %>% unname()
  if (!is.null(SFTA3_i)) {
    sym[SFTA3_i] <- "SFTA3"
  }
  return(sym)
}


readcounts <- function(file, sep = ",", header = T, row.names = 1) {
  return(as.matrix(read.table(file, sep = sep, header = header, row.names = row.names)))
}
