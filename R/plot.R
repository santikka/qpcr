fold_change_plots <- function(data, ref_genes) {
  has_ref <- !is.null(ref_genes)
  pairs <- utils::combn(unique(data$.group), 2L)
  genes <- setdiff(names(data), c(".group", ref_genes))
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
    if (has_ref) {
      ref_x <- data[idx_x, rowMeans(.SD), .SDcols = ref_genes]
      ref_y <- data[idx_y, rowMeans(.SD), .SDcols = ref_genes]
    }
    data_pairs[idx, (gene_cols) := genes_x - ref_x - genes_y + ref_y]
  }
  data_long <- data.table::melt(
    data_pairs,
    measure.vars = gene_cols,
    variable.name = "gene",
    value.name = "logct"
  )
  data_fold <- data_long[,
    {
      mean <- mean(logct)
      sd <- sd(logct)
      list(mean = exp(mean), lo = exp(mean - sd), hi = exp(mean + sd))
    },
    by = c("pair", "gene")
  ]
  data_log_fold <- data_long[,
    {
      mean <- ifelse_(has_ref, mean(logct), mean(logct / log(2.0)))
      sd <- ifelse_(has_ref, sd(logct), sd(logct / log(2.0)))
      list(mean = mean, lo = mean - sd, hi = mean + sd)
    },
    by = c("pair", "gene")
  ]
  # avoid NSE warnings in R CMD check
  gene <- pair <- lo <- hi <- NULL
  plot_fold <-
    ggplot2::ggplot(data_fold, ggplot2::aes(x = gene, y = mean, color = pair)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(0.75)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lo, ymax = hi),
      position = ggplot2::position_dodge(0.75)
    ) +
    ggplot2::xlab("") +
    ggplot2::ylab(
      expression(
        "Fold Change" * phantom(0) * 2^{-Delta * Delta * "Ct"} * phantom(0) *
          bgroup("(", 2^{"mean " %+-% " SD"}, ")")
      )
    ) +
    ggplot2::geom_hline(yintercept = 1.0, alpha = 0.33) +
    ggplot2::scale_y_continuous(
      limits = c(0.0, NA),
      expand = ggplot2::expansion(mult = c(0, 0.1))
    ) +
    ggplot2::scale_color_discrete(name = "") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    )
  plot_log_fold <-
    ggplot2::ggplot(
      data_log_fold,
      ggplot2::aes(x = gene, y = mean, color = pair)
    ) +
    ggplot2::geom_point(position = ggplot2::position_dodge(0.75)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lo, ymax = hi),
      position = ggplot2::position_dodge(0.75)
    ) +
    ggplot2::xlab("") +
    ggplot2::ylab(
      expression(
        "Log"[2] * " Fold Change" * phantom(0)
          -Delta * Delta * "Ct" * phantom(0) * "(mean " %+-% " SD)"
      )
    ) +
    ggplot2::geom_hline(yintercept = 0.0, alpha = 0.33) +
    ggplot2::scale_color_discrete(name = "") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    )

  list(fold = plot_fold, log_fold = plot_log_fold)
}
