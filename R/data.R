parse_data <- function(data, treatment) {
  data <- data.table::as.data.table(data.table::copy(data))
  data.table::set(
    x = data,
    j = treatment,
    value = as.factor(data[[treatment]])
  )
  cols <- names(data)
  gene_cols <- setdiff(cols, treatment)
  col_types <- vapply(data[, .SD, .SDcols = gene_cols], typeof, character(1L))
  valid_cols <- col_types == "double"
  stopifnot_(
    all(valid_cols),
    c(
      "Column{?s} {.var {gene_cols[!valid_cols]}} of {.arg data}
      {?is/are} invalid:",
      `x` = "Column type{?s} {.cls {col_types[!valid_cols]}}
             {?is/are} not supported."
    )
  )
  finite_cols <- vapply(
    data,
    function(x) all(is.finite(x)),
    logical(1L)
  )
  stopifnot_(
    all(finite_cols),
    "Non-finite or missing values were found in variable{?s}
    {.var {cols[!finite_cols]}} of {.arg data}."
  )
  data
}
