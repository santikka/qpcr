#' Tests For Each Gene Between All Group Pairs
#'
#' @inheritParams qpcr
#' @noRd
local_tests <- function(data, group, norm, methods, ref_genes, m, alpha) {
  methods <- methods[methods %in% local_methods]
  pairs <- utils::combn(unique(data[[group]]), 2L)
  gene_cols <- setdiff(names(data), c(group, ref_genes))
  n_pairs <- ncol(pairs)
  n_genes <- length(gene_cols)
  out <- data.table::data.table(
    gene = rep(gene_cols, each = n_pairs),
    group_a = rep(pairs[1L, ], n_genes),
    group_b = rep(pairs[2L, ], n_genes),
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
    idx_x <- data[[group]] == pairs[1L, i]
    idx_y <- data[[group]] == pairs[2L, i]
    idx <- i - n_pairs
    ref_x <- 0.0
    ref_y <- 0.0
    if (!is.null(ref_genes)) {
      ref_x <- data[idx_x, rowMeans(.SD), .SDcols = ref_genes]
      ref_y <- data[idx_y, rowMeans(.SD), .SDcols = ref_genes]
    }
    for (j in seq_len(n_genes)) {
      idx <- idx + n_pairs
      target_x <- data[idx_x, ][[gene_cols[j]]]
      target_y <- data[idx_y, ][[gene_cols[j]]]
      if ("randomization_test" %in% methods) {
        tmp <- randomization_test(
          x = target_x - ref_x,
          y = target_y - ref_y,
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
          x = target_x - ref_x,
          y = target_y - ref_y,
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
          x = target_x - ref_x,
          y = target_y - ref_y,
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

#' Tests For Each Gene Over All Groups
#'
#' @inheritParams qcpr
#' @noRd
global_tests <- function(data, group, norm, methods, ref_genes, alpha, ...) {
  methods <- methods[methods %in% global_methods]
  gene_cols <- setdiff(names(data), c(group, ref_genes))
  n_genes <- length(gene_cols)
  refs <- 0.0
  if (identical(norm, "reference")) {
    refs <- data[, .SD, .SDcols = c(ref_genes, group)][,
      list(ref = rowMeans(.SD)), by = group
    ]$ref
  }
  grp <- data[[group]]
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
  out_kruskal <- data.table::data.table(
    gene = gene_cols,
    stat = NA_real_,
    df = NA_real_,
    p = NA_real_
  )
  for (i in seq_along(gene_cols)) {
    gene <- data[[gene_cols[i]]] - refs
    if ("anova" %in% methods) {
      tmp <- anova(lm(gene ~ grp))
      data.table::set(out_anova, i = i, j = "df_trt", value = tmp$Df[1L])
      data.table::set(out_anova, i = i, j = "df_res", value = tmp$Df[2L])
      data.table::set(out_anova, i = i, j = "ss_trt", value = tmp$Sum[1L])
      data.table::set(out_anova, i = i, j = "ss_res", value = tmp$Sum[2L])
      data.table::set(out_anova, i = i, j = "ms_trt", value = tmp$Mean[1L])
      data.table::set(out_anova, i = i, j = "ms_res", value = tmp$Mean[2L])
      data.table::set(out_anova, i = i, j = "stat",   value = tmp$F[1L])
      data.table::set(out_anova, i = i, j = "p",      value = tmp$Pr[1L])
    }
    if ("kruskal_test" %in% methods) {
      tmp <- kruskal.test(gene ~ grp)
      data.table::set(out_kruskal, i = i, j = "stat", value = tmp$statistic)
      data.table::set(out_kruskal, i = i, j = "df",   value = tmp$parameter)
      data.table::set(out_kruskal, i = i, j = "p",    value = tmp$p.value)
    }
  }
  list(out_anova = out_anova, out_kruskal = out_kruskal)
}

#' Pairwise Randomization Test Between Two Groups
#'
#' @export
#' @inheritParams qpcr
#' @param x \[`numeric()`]\cr A vector of (efficiency adjusted)
#'   Ct values for the target gene in group A.
#' @param y \[`numeric()`]\cr A vector of (efficiency adjusted)
#'   Ct values for the target gene in group B.
#' @noRd
randomization_test <- function(x, y, m, alpha) {
  xy <- c(x, y)
  n_x <- length(x)
  n_y <- length(y)
  n_xy <- n_x + n_y
  n_max <- max(n_x, n_y)
  comb <- utils::combn(seq_len(n_xy), 2L)
  pairs <- cbind(comb, n_xy + 1L - comb)
  n_pairs <- ncol(pairs)
  idx_x <- sample.int(n_x, m, replace = TRUE)
  idx_y <- sample.int(n_y, m, replace = TRUE)
  log_expr_obs <- mean(x) - mean(y)
  log_expr_sim <- numeric(m)
  x_seq <- seq_len(n_x)
  y_seq <- seq_len(n_y)
  for (i in seq_len(m)) {
    cols <- sample.int(n_pairs, n_max, replace = TRUE)
    idx_x <- pairs[1L, cols][x_seq]
    idx_y <- pairs[2L, cols][y_seq]
    log_expr_sim[i] <- mean(xy[idx_x]) - mean(xy[idx_y])
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
