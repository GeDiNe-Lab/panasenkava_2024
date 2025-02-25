function(
    ma, metadata, minc = 15, summarize = "merge", time = "time",
    col = NULL, consensusCluster = FALSE, reduce = FALSE, cutoff = 0.7,
    scale = TRUE, pattern = NULL, groupDifference = NULL, eachStep = FALSE,
    plot = TRUE, fixy = NULL, nClusters = NULL, skipDendrogram = TRUE) {
    benchmarking <- NULL
    metadata <- as.data.frame(metadata)
    ma <- ma[, row.names(metadata)]
    rownames(ma) <- make.names(rownames(ma))
    if (is.null(col)) {
        col <- "colored"
        metadata[, col] <- rep("one_group", nrow(metadata))
    }
    if (!summarize %in% names(metadata)) {
        metadata[, summarize] <- as.factor(paste0(
            metadata[, col],
            metadata[, time]
        ))
    }
    metadata[, summarize] <- droplevels(metadata[, summarize])
    stopifnot(class(metadata)[1] == "data.frame")
    stopifnot(class(ma)[1] == "matrix" | class(ma)[1] == "data.frame")
    stopifnot(summarize %in% names(metadata))
    stopifnot(time %in% names(metadata))
    ma <- as.matrix(ma)
    if (!is.null(fixy)) {
        stopifnot(length(fixy) == 2)
    }
    if (nrow(ma) > 3000 & is.null(pattern)) {
        message(
            "A large number of genes was given-- please, ",
            "make sure this is not an error. Normally, ", "only DE genes will be useful for this function."
        )
    }
    message("Working with ", nrow(ma), " genes.")
    counts_group <- .summarize_scale(
        ma, metadata[[summarize]],
        FALSE
    )
    if (!is.null(groupDifference)) {
        counts_group <- .remove_low_difference(
            counts_group,
            groupDifference, eachStep
        )
    }
    if (scale) {
        norm_sign <- t(apply(counts_group, 1, .scale))
    } else {
        norm_sign <- counts_group
    }
    colnames(norm_sign) <- colnames(counts_group)
    metadata_groups <- metadata %>% dplyr::distinct(!!sym(summarize),
        .keep_all = TRUE
    )
    rownames(metadata_groups) <- metadata_groups[, summarize]
    norm_sign <- norm_sign[, row.names(metadata_groups), drop = FALSE]
    if (nrow(ma) == 1) {
        p <- .plot_cluster(
            norm_sign, as.character(rownames(norm_sign)),
            metadata_groups[, time], metadata_groups[, col],
            rownames(norm_sign), fixy
        )
        cluster_genes <- c()
        groups <- as.factor("1")
        names(groups) <- c(row.names(ma))
    } else if (!consensusCluster & is.null(pattern)) {
        cluster_genes <- .make_clusters(counts_group)
        groups <- .select_genes(cluster_genes, norm_sign, minc,
            reduce = reduce, cutoff = cutoff, nClusters = nClusters
        )
        benchmarking <- .benckmark_cutoff(
            cluster_genes, norm_sign,
            minc
        )
    } else if (consensusCluster & is.null(pattern)) {
        cluster_genes <- .make_concensus_cluster(counts_group)
        groups <- .select_concensus_genes(cluster_genes)
    } else if (is.character(pattern)) {
        stopifnot(pattern %in% rownames(counts_group))
        if (length(pattern) > 1) {
            reference <- rowMedians(counts_group[pattern, ])
        } else {
            reference <- counts_group[pattern, ]
        }
        cluster_genes <- .find_pattern(counts_group, reference)
        groups <- .select_pattern(cluster_genes)
    } else if (is.numeric(pattern)) {
        stopifnot(length(pattern) == ncol(counts_group))
        cluster_genes <- .find_pattern(counts_group, pattern)
        groups <- .select_pattern(cluster_genes)
    }
    temp <- names(groups)
    dend_plot <- NA
    if (length(unique(groups)) > 0 & is.null(nClusters) & !skipDendrogram) {
        dend <- cluster_genes
        h <- dend$dc
        clust <- cutree(as.hclust(dend), h = h)
        clust.cutree <- dendextend:::cutree(dend, h = h, order_clusters_as_data = FALSE)
        dend <- as.dendrogram(dend, h = h)
        idx <- order(names(clust.cutree))
        clust.cutree <- clust.cutree[idx]
        df.merge <- merge(clust, clust.cutree, by = "row.names")
        df.merge.sorted <- df.merge[order(df.merge$y), ]
        lbls <- unique(df.merge.sorted$x)
        dend_plot <- dendextend::color_branches(dend,
            h = h,
            groupLabels = TRUE, warn = FALSE
        ) %>%
            dendextend::set(
                "labels",
                ""
            ) %>%
            suppressWarnings()
        groups <- match(groups, lbls)
        if (plot) {
            plot(dend_plot,
                xlab = "", ylab = "", main = "",
                sub = "", axes = FALSE, cex = 2
            )
        }
    }
    if (length(unique(groups)) > 0 & is.numeric(nClusters) &
        !skipDendrogram) {
        dend <- cluster_genes
        clust <- cutree(as.hclust(dend), k = nClusters)
        clust.cutree <- dendextend:::cutree(dend,
            k = nClusters,
            order_clusters_as_data = FALSE
        )
        dend <- as.dendrogram(dend, k = nClusters)
        idx <- order(names(clust.cutree))
        clust.cutree <- clust.cutree[idx]
        df.merge <- merge(clust, clust.cutree, by = "row.names")
        df.merge.sorted <- df.merge[order(df.merge$y), ]
        lbls <- unique(df.merge.sorted$x)
        dend_plot <- dendextend::color_branches(dend,
            k = nClusters,
            groupLabels = TRUE, warn = FALSE
        ) %>%
            dendextend::set(
                "labels",
                ""
            ) %>%
            suppressWarnings()
        groups <- match(groups, lbls)
        if (plot) {
            plot(dend_plot,
                xlab = "", ylab = "", main = "",
                sub = "", axes = FALSE, cex = 2
            )
        }
    }
    names(groups) <- temp
    df <- data.frame(
        genes = names(groups), cluster = groups,
        stringsAsFactors = FALSE
    )
    raw <- counts_group %>%
        as.data.frame() %>%
        rownames_to_column("genes") %>%
        gather(!!sym(summarize), "value", -genes) %>%
        inner_join(metadata_groups %>%
            mutate_if(is.factor, as.character)) %>%
        inner_join(df,
            by = "genes"
        )
    summarise <- raw %>%
        group_by(
            !!sym(summarize), !!sym("cluster"),
            !!sym(time), !!sym(col)
        ) %>%
        summarise(
            abundance = median(value),
            n_genes = n()
        ) %>%
        ungroup()
    normalized <- norm_sign %>%
        as.data.frame() %>%
        rownames_to_column("genes") %>%
        gather(!!sym(summarize), "value", -genes) %>%
        inner_join(metadata_groups %>%
            mutate_if(is.factor, as.character)) %>%
        inner_join(df,
            by = "genes"
        )
    if (!is.null(benchmarking)) {
        normalized <- normalized %>% left_join(benchmarking[["genes"]],
            by = "genes"
        )
    }
    normalized[[time]] <- factor(normalized[[time]], levels = levels(metadata[[time]]))
    plot_benchmarking <- .plot_benchmarking(
        normalized, benchmarking,
        time, col
    )
    plot_benchmarking_curve <- .plot_benchmarking_curve(benchmarking)
    plotting_data <- list(
        norm = normalized, time = time, col = col,
        min_genes = minc
    )
    if (length(unique(groups)) > 0) {
        p <- degPlotCluster(normalized, time, col, min_genes = minc)
        if (!is.null(fixy)) {
            p <- p + ylim(fixy[1], fixy[2])
        }
        if (plot) {
            print(p)
        }
    }
    invisible(list(
        df = df, pass = unique(groups), plot = p,
        dend = dend_plot, hr = cluster_genes, normalized = normalized,
        summarise = summarise, raw = raw, counts = ma[raw[["genes"]], ], benchmarking = plot_benchmarking, benchmarking_curve = plot_benchmarking_curve
    ))
}
