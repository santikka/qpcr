normagene <- function(x, treatment, log = FALSE) {
  out <- x[, log(.SD), by = treatment][,
    rowMeans(colMeans(.SD, na.rm = TRUE)[col(.SD)] - .SD, na.rm = TRUE) + .SD,
    by = treatment]
  if (log) {
    out
  } else {
    out[, exp(.SD), by = treatment]
  }
}
