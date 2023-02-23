normagene <- function(x, treatment) {
  out <- x[, log(.SD), by = treatment][,
    exp(
      rowMeans(colMeans(.SD, na.rm = TRUE)[col(.SD)] - .SD, na.rm = TRUE) +
        .SD
    ),
    by = treatment]
}
