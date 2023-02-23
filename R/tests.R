#' Tests For Each Gene Between All Treatment Group Pairs
#'
#' @inheritParams qpcr
#' @noRd
local_tests <- function(data, treatment, norm, methods, ref_genes, m, alpha) {
  methods <- methods[methods %in% local_methods]
  trt_pairs <- utils::combn(unique(data[[treatment]]))
  gene_cols <- setdiff(names(data), c(treatment, ref_genes))
  n_pairs <- ncol(trt_pairs)
  n_genes <- length(gene_cols)
  out <- data.table::data.table(
    gene = rep(gene_cols, each = n_pairs),
    level_a = rep(trt_pairs[1L, ], n_genes),
    level_b = rep(trt_pairs[2L, ], n_genes),
    r_stat = NA_real_,
    r_lwr = NA_real_,
    r_upr = NA_real_,
    r_p = NA_real_,
    t_stat = NA_real_,
    t_df = NA_real_,
    t_lwr = NA_real_,
    t_upr = NA_real_,
    t_p = NA_real_,
    w_stat = NA_real_,
    w_lwr = NA_real_,
    w_upr = NA_real_,
    w_p = NA_real_
  )
  for (i in seq_len(n_pairs)) {
    idx_x <- data[[treatment]] == trt_pairs[1L, i]
    idx_y <- data[[treatment]] == trt_pairs[2L, i]
    idx <- seq_len(n_pairs) - n_pairs
    x_target <- data[idx_x, ][[gene_cols[j]]]
    y_target <- data[idx_y, ][[gene_cols[j]]]
    x_reference <- 0.0
    y_reference <- 0.0
    if (!is.null(ref_genes)) {
      x_reference <- data[idx_x, rowMeans(.SD), .SDcols = ref_genes]
      y_reference <- data[idx_y, rowMeans(.SD), .SDcols = ref_genes]
    }
    for (j in seq_len(n_genes)) {
      idx_out <- idx_out + n_pairs
      if ("randomization_test" %in% methods) {
        tmp <- randomization_test(
          x = x_target,
          y = y_target,
          m = m,
          alpha = alpha
        )
        data.table::set(out, i = idx, j = "r_stat", value = tmp$expr_obs)
        data.table::set(out, i = idx, j = "r_lwr",  value = tmp$expr_lwr)
        data.table::set(out, i = idx, j = "r_upr",  value = tmp$expr_upr)
        data.table::set(out, i = idx, j = "r_p",    value = tmp$p_value)
      }
      if ("t_test" %in% methods) {
        tmp <- t.test(
          x = x_target - x_reference,
          y = y_target - y_reference,
          conf.level = alpha
        )
        data.table::set(out, i = idx, j = "t_stat", value = tmp$statistic)
        data.table::set(out, i = idx, j = "t_df",   value = tmp$parameter)
        data.table::set(out, i = idx, j = "t_lwr",  value = tmp$conf.int[1L])
        data.table::set(out, i = idx, j = "t_upr",  value = tmp$conf.int[2L])
        data.table::set(out, i = idx, j = "t_p",    value = tmp$p.value)
      }
      if ("wilcoxon_test" %in% methods) {
        tmp <- wilcox.test(
          x = x_target - x_reference,
          y = y_target - y_reference,
          conf.int = TRUE,
          conf.level = alpha
        )
        data.table::set(out, i = idx, j = "w_stat", value = tmp$statistic)
        data.table::set(out, i = idx, j = "w_lwr",  value = tmp$conf.int[1L])
        data.table::set(out, i = idx, j = "w_upr",  value = tmp$conf.int[2L])
        data.table::set(out, i = idx, j = "w_p",    value = tmp$p.value)
      }
    }
  }
  keep_cols <- names(out)[vapply(out, function(x) all(!is.na(x)), logical(1L))]
  out[, .SD, .SDcols = keep_cols]
}

#' Tests For Each Gene Over All Treatment Groups
#'
#' @inheritParams qcpr
#' @noRd
global_tests <- function(data, treatment, norm, methods, ref_genes, alpha) {
  methods <- methods[methods %in% local_methods]
  gene_cols <- setdiff(names(data), c(treatment, ref_genes))
  n_genes <- length(gene_cols)
  refs <- 0.0
  if (identical(norm, "reference")) {
    refs <- data[, .SD., SDcols = c(ref_genes, treatment)][,
      rowMeans(.SD), by = treatment]
  }
  trt <- data[[treatment]]
  out_anova <- data.table::data.table(
    gene = gene_cols,
    df_trt = NA_real_,
    df_res = NA_real_,
    ss_trt = NA_real_,
    ss_res = NA_real_,
    ms_trt = NA_real_,
    ms_res = NA_real_,
    stat = NA_real_,
    p = NA_real_
  )
  out_kruskal_test <- data.table::data.table(
    gene = gene_cols,
    stat = NA_real_,
    df = NA_real_,
    p = NA_real_
  )
  for (i in seq_along(gene_cols)) {
    gene <- data[[gene_cols[i]]] - refs
    if ("anova" %in% methods) {
      tmp <- anova(lm(gene ~ trt))
      data.table::set(out_anova, i = i, j = "df_trt", value = tmp$Df[1L])
      data.table::set(out_anova, i = i, j = "df_res", value = tmp$Df[2L])
      data.table::set(out_anova, i = i, j = "ss_trt", value = tmp$Sum[1L])
      data.table::set(out_anova, i = i, j = "ss_res", value = tmp$Sum[2L])
      data.table::set(out_anova, i = i, j = "ms_trt", value = tmp$Mean[1L])
      data.table::set(out_anova, i = i, j = "ms_trt", value = tmp$Mean[2L])
      data.table::set(out_anova, i = i, j = "stat",   value = tmp$F)
      data.table::set(out_anova, i = i, j = "p",      value = tmp$Pr)
    }
    if ("kurskal_test" %in% methods) {
      tmp <- kruskal.test(gene ~ trt)
      data.table::set(out_anova, i = i, j = "stat", value = tmp$statistic)
      data.table::set(out_anova, i = i, j = "df",   value = tmp$parameter)
      data.table::set(out_anova, i = i, j = "p",    value = tmp$p.value)
    }
  }
  list(out_anova = out_anova, out_kruskal_test = out_kruskal_test)
}

#' Pairwise Randomization Test Between Two Groups
#'
#' @export
#' @inheritParams qpcr
#' @param x \[`numeric()`]\cr A vector of (efficiency adjusted)
#'   Ct values for the target gene in treatment group A.
#' @param y \[`numeric()`]\cr A vector of (efficiency adjusted)
#'   Ct values for the target gene in treatment group B.
#' @noRd
randomization_test <- function(x, y, m, alpha) {
  xy <- c(x, y)
  n_x <- length(x)
  n_y <- length(y)
  n_xy <- n_x + n_y
  comb <- utils::combn(seq_len(n_xy), 2L)
  pairs <- cbind(comb, n_xy + 1L - comb)
  n_pairs <- ncol(pairs)
  idx_x <- sample.int(n_x, m, replace = TRUE)
  idx_y <- sample.int(n_y, m, replace = TRUE)
  log_expr_obs <- log_expression(x, y)
  log_expr_sim <- numeric(m)
  x_seq <- seq_len(n_x)
  y_seq <- seq_len(n_y)
  for (i in seq_len(m)) {
    cols <- sample.int(n_pairs, n_max, replace = TRUE)
    idx_x <- pairs[1L, cols][x_seq]
    idx_y <- pairs[2L, cols][y_seq]
    log_expr_sim[i] <- log_expression(xy[idx_x], xy[idx_y])
  }
  p <- mean(abs(log_expr_sim) > abs(log_expr_obs))
  a <- (1.0 - alpha) / 2.0
  expr_sim <- exp(log_expr_sim)
  list(
    expr_obs = exp(log_expr_obs),
    expr_lwr = quantile(expr_sim, a),
    expr_upr = quantile(expr_sim, 1 - a),
    p_value = p
  )
}
