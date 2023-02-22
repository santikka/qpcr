normalize_normagene <- function(x, treatment) {
  out <- x[, log(.SD), by = treatment][,
    rowMeans(colMeans(.SD, na.rm = TRUE)[col(.SD)] - .SD, na.rm = TRUE) + .SD,
    by = treatment]
}

normalize_reference <- function(x, treatment, ref_genes) {
  x[,]
}

normalize_none <- function(x, treatment) {
  x[, log(.SD), by = treatment]
}
