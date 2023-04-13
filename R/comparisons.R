
#' Construct Group/Population Comparison Pairs
#'
#' @inheritParams qpcr
#' @noRd
parse_comparisons <- function(data, group, ref_genes comparisons) {
  n_group <- length(group)
  genes <- setdiff(names(data), union(".group", ref_genes))
  levels_grid <- expand.grid(gene = genes, group = levels(data$.group))
  levels <- paste(levels_grid$gene, levels_grid$group, sep = ":")
  all_pairs <- utils::combn(levels, 2L)
  matched_pairs <- logical(ncol(all_pairs))
  comparisons <- ifelse_(
    is.null(comparisons),
    paste0(
      paste(rep(":", 1L + n_group), collapse = ""),
      "/",
      paste(rep(":", 1L + n_group), collapse = "")
    ),
    comparisons
  )
  for (i in seq_along(comparisons)) {
    m <- regexec(
      pattern = "([^/]+)/([^/]+)",
      text = comparisons[i],
      perl = TRUE
    )
    match <- regmatches(comparisons[i], m)[[1L]][-1L]
    stopifnot_(
      !identical(m[[1L]], -1),
      "Invalid comparison"
    )
    lhs <- parse_comparison(match[1L], levels, n_group)
    rhs <- parse_comparison(match[2L], levels, n_group)
    pairs <- expand.grid(lhs, rhs)
    for (j in seq_len(nrow(pairs))) {
      idx <- c(
        which(all_pairs[1, ] == pairs[j, 1] & all_pairs[2, ] == pairs[j, 2]),
        which(all_pairs[2, ] == pairs[j, 1] & all_pairs[1, ] == pairs[j, 2])
      )
      matched_pairs[idx] <- TRUE
    }
  }
  all_pairs[, matched_pairs]
}

#' Parse A Single Comparison
#'
#' @param x
#' @param genes
#' @param groups
#' @param populations
#' @noRd
parse_comparison <- function(x, levels, n_group) {
  pattern <- paste(rep("([^:]*)", n_group + 1L), collapse = ":")
  m <- regexec(pattern = pattern, text = x, perl = TRUE)
  stopifnot_(
    !identical(m[[1L]], -1),
    "Invalid comparison"
  )
  match <- regmatches(x, m)[[1L]][-1L]
  match[nzchar(match)] <- "[^:]+"
  pattern <- paste0(match)
  m <- regexec(pattern = pattern, text = levels, perl = TRUE)
  unlist(regmatches(levels, m))
}
