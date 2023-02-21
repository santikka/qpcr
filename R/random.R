# Randomization Test Without Reference Gene -------------------------------

#' Pairwise Randomization Test For Each Gene Between All Groups
#'
#' @inheritParams qpcr
#' @param m \[`integer(1)`]\cr Number of replications to use.
randomization_test <- function(data, treatment, m) {
  trt_pairs <- utils::combn(unique(data[[treatment]]))
  gene_cols <- setdiff(names(data), treatment)
  n_pairs <- ncol(trt_pairs)
  n_genes <- length(gene_cols)
  out <- data.frame(
    gene = rep(gene_cols, each = n_pairs),
    level_a = trt_pairs[1L, ],
    level_b = trt_pairs[2L, ],
    expression = NA_real_,
    pvalue = NA_real_
  )
  for (i in seq_len(n_pairs * n_genes)) {
    idx_x <- data[[treatment]] == out$level_a[i]
    idx_y <- data[[treatment]] == out$level_b[i]
    tmp <- randomize(
      x = data[idx_x, ][[out$gene[i]]],
      y = data[idx,y, ][[out$gene[i]]],
      m = m
    )
    out$expression[i] <- tmp$expr_obs
    out$pvalue[i] <- tmp$p
  }
}

#' Pairwise Randomization Test Between Two Groups
#'
#' @param x \[`numeric()`]\cr A vector of Ct values for treatment group A.
#' @param y \[`numeric()`]\cr A vector of Ct values for treatment group B.
#' @param m \[`integer(1)`]\cr Number of replications to use.
#' @noRd
randomize <- function(x, y, m) {
  xy <- c(x, y)
  n_x <- length(x)
  n_y <- length(y)
  n_xy <- n_x + n_y
  comb <- utils::combn(seq_len(n_xy), 2L)
  pairs <- cbind(comb, n_xy + 1L - comb)
  n_pairs <- ncol(pairs)
  idx_x <- sample.int(n_x, m, replace = TRUE)
  idx_y <- sample.int(n_y, m, replace = TRUE)
  log_expr_obs <- mean(x) - mean(y)
  #expr_sim <- 2.0^(mean(x[idx_x]) - mean(y[idx_y]))
  #d_obs <- 2.0^(x[idx_x] - y[idx_y])
  log_expr_sim <- numeric(m)
  for (i in seq_len(m)) {
    cols <- sample.int(n_pairs, n_max, replace = TRUE)
    idx_x <- pairs[1L, cols][seq_len(n_x)]
    idx_y <- pairs[2L, cols][seq_len(n_y)]
    log_expr_sim[i] <- mean(xy[idx_x]) - mean(xy[idx_y])
  }
  p <- mean(abs(log_expr_sim) > abs(log_expr_obs))
  list(
    expr_obs = exp(log_expr_obs),
    #expr_sim = exp(log_expr_sim),
    #d_obs = d_obs,
    #d_hyp = d_hyp,
    p = p
  )
}

# Randomization Test With Reference Gene ----------------------------------

#' Pairwise Randomization Test For Each Non-Reference Gene Between All Groups
#'
#' @inheritParams qpcr
#' @param m \[`integer(1)`]\cr Number of replications to use.
randomization_test_ref <- function(data, treatment, ref_genes, eff, m) {
  trt_pairs <- utils::combn(unique(data[[treatment]]))
  gene_cols <- setdiff(names(data), treatment)
  n_pairs <- ncol(trt_pairs)
  n_genes <- length(gene_cols)
  out <- data.frame(
    gene = rep(gene_cols, each = n_pairs),
    level_a = trt_pairs[1L, ],
    level_b = trt_pairs[2L, ],
    expression = NA_real_,
    pvalue = NA_real_
  )
  for (i in seq_len(n_pairs * n_genes)) {
    idx_x <- data[[treatment]] == out$level_a[i]
    idx_y <- data[[treatment]] == out$level_b[i]
    tmp <- randomize(
      x = data[idx_x, ][[out$gene[i]]],
      y = data[idx,y, ][[out$gene[i]]],
      m = m
    )
    out$expression[i] <- tmp$expr_obs
    out$pvalue[i] <- tmp$p
  }
}

#' Pairwise Randomization Test Between Two Groups Using a Reference Gene
#'
#' @param x_target \[`numeric()`]\cr
#'   A vector of Ct values of the target gene for treatment group A.
#' @param y_target \[`numeric()`]\cr
#'   A vector of Ct values of the target gene for treatment group B.
#' @param x_ref \[`numeric()`]\cr
#'   A vector of Ct values of the reference gene for treatment group A.
#' @param y_ref \[`numeric()`]\cr
#'   A vector of Ct values of the reference gene for treatment group B.
#' @param eff_x \[`numeric(1)`]\cr Efficiency value for treatment group A.
#' @param eff_y \[`numeric(1)`]\cr Efficiency value for treatment group B.
#' @param m \[`integer(1)`]\cr Number of replications to use.
#' @noRd
randomize_ref <- function(x_target, y_target, x_ref, y_ref, eff_x, eff_y, m) {

}
