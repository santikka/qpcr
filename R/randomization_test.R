#' Pairwise Randomization Test
#'
#' @param x \[`numeric()`] A vector of Ct values for treatment group A.
#' @param y \[`numeric()`] A vector of Ct values for treatment group B.
#' @param m \[`integer(1)`] Number of replications to use.
#' @noRd
randomization_test <- function(x, y, m) {
  xy <- c(x, y)
  n_x <- length(x)
  n_y <- length(y)
  n_xy <- n_x + n_y
  n_max <- max(n_a, n_b)
  pairs <- expand.grid(seq_len(n_xy), seq_len(n_xy))
  pairs <- pairs[pairs[, 1L] != pairs[, 2L], ]
  n_pairs <- nrow(pairs)
  idx_x <- sample.int(n_x, m, replace = TRUE)
  idx_y <- sample.int(n_y, m, replace = TRUE)
  eff_obs <- 2.0^(mean(x) - mean(y))
  eff_sim <- 2.0^(mean(x[idx_x]) - mean(y[idx_y]))
  d_obs <- 2.0^(x[idx_x] - y[idx_y])
  d_hyp <- numeric(m)
  for (i in seq_len(m)) {
    rows <- sample.int(n_rows, n_max, replace = TRUE)
    idx_x <- pairs[rows, 1L][seq_len(n_x)]
    idx_y <- pairs[rows, 2L][seq_len(n_y)]
    d_hyp[i] <- 2.0^(mean(xy[idx_x]) - mean(xy[idx_y]))
  }
  p <- mean(abs(log(d_hyp)) > abs(log(eff_obs)))
  list(
    eff_obs = eff_obs,
    eff_sim = eff_sim,
    d_obs = d_obs,
    d_hyp = d_hyp,
    p = p
  )
}
