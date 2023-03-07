fold_change <- function(data, group, genes, ref_genes) {
  pairs <- utils::combn(unique(data[[group]]), 2L)
  gene_cols <- setdiff(names(data), c(group, ref_genes))
  n_pairs <- ncol(pairs)
  n_data <- data[, .N, by = group]$N[1L]
  data_pairs <- data.table::data.table(
    group_a = rep(pairs[1L, ], each = n_data),
    group_b = rep(pairs[2L, ], each = n_data)
  )
  data_pairs[ , pair := paste0(group_b, " vs ", group_a)]
  idx <- seq_len(n_data) - n_data
  for (i in seq_len(n_pairs)) {
    idx <- idx + n_data
    idx_x <- data[[group]] == pairs[1L, i]
    idx_y <- data[[group]] == pairs[2L, i]
    genes_x <- data[idx_x, .SD, .SDcols = gene_cols]
    genes_y <- data[idx_y, .SD, .SDcols = gene_cols]
    ref_x <- 0.0
    ref_y <- 0.0
    if (!is.null(ref_genes)) {
      ref_x <- data[idx_x, rowMeans(.SD), .SDcols = ref_genes]
      ref_y <- data[idx_y, rowMeans(.SD), .SDcols = ref_genes]
    }
    data_pairs[idx, (gene_cols) := genes_x - ref_x - genes_y + ref_y]
  }
  data_long <- melt(
    data_pairs,
    measure.vars = gene_cols,
    variable.name = "gene",
    value.name = "logct"
  )
  data_summ <- data_long[,
    {
      mean <- exp(mean(logct))
      sd <- sd(exp(logct))
      list(mean = mean, lo = mean - sd, hi = mean + sd)
    },
    by = c("pair", "gene")
  ]
  # avoid NSE warnings in R CMD check
  gene <- pair <- lo <- hi <- NULL
  ggplot2::ggplot(data_summ, ggplot2::aes(x = gene, y = mean, color = pair)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(0.75)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lo, ymax = hi),
      position = ggplot2::position_dodge(0.75)
    ) +
    ggplot2::xlab("") +
    ggplot2::ylab("Fold Change (mean +- SD)") +
    ggplot2::geom_hline(yintercept = 1.0, alpha = 0.33)
}
