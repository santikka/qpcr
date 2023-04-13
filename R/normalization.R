normagene <- function(x) {
  x[, log(.SD), by = ".group"][,
    exp(
      rowMeans(colMeans(.SD, na.rm = TRUE)[col(.SD)] - .SD, na.rm = TRUE) + .SD
    ),
    by = ".group"
  ]
}
