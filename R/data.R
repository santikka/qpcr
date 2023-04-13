#' Parse Input Data For Analysis
#'
#' @inheritParams qpcr
#' @noRd
parse_data <- function(data, group, norm, ref_genes, eff, eff_adjust) {
  stopifnot_(
    is.data.frame(data),
    "Argument {.arg data} must be a {.cls data.frame} object."
  )
  stopifnot_(
    checkmate::test_character(x = group),
    "Argument {.arg group} must be a {.cls character} vector."
  )
  stopifnot_(
    is.null(ref_genes) || checkmate::test_character(x = ref_genes),
    "Argument {.arg ref_genes} must be a {.cls character} vector."
  )
  cols <- names(data)
  invalid <- group[!group %in% cols]
  stopifnot_(
    identical(length(invalid), 0L),
    c(
      "Argument {.arg group} must contain only column names of {.arg data}:",
      `x` = "Columns {.var {invalid}} {?does/do} not exist in {.arg data}."
    )
  )
  data <- data.table::as.data.table(data.table::copy(data))
  data.table::set(
    x = data,
    j = ".group",
    value = interaction(data[, c(group)], sep = ":")
  )
  data <- data[, .SD, .SDcols = c(setdiff(data_names, group), ".group")]
  cols <- names(data)
  invalid <- ref_genes[!reg_genes %in% cols]
  stopifnot_(
    is.null(ref_genes) || identical(length(invalid), 0L),
    c(
      "Argument {.arg ref_genes} contains unknown columns:",
      `x` = "Column{?s} {.val {invalid}} {?does/do} not exist in {.arg data}."
    )
  )
  invalid <- intersect(group, ref_genes)
  stopifnot_(
    is.null(ref_genes) || identical(length(invalid), 0L),
    "Grouping column{?s} {.val {invalid}} {?is/are}
     also defined in {.arg ref_genes}."
  )
  genes <- setdiff(cols, union(".group", ref_genes))
  gene_cols <- c(genes, ref_genes)
  col_types <- vapply(
    data[, .SD, .SDcols = gene_cols],
    typeof,
    character(1L)
  )
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
  if (identical(norm, "normagene")) {
    data <- normagene(data)
  }
  parse_efficiency(data, norm, ref_genes, eff, eff_adjust)
}

#' Parse Efficiency Values For Analysis
#'
#' @inheritParams qpcr
#' @noRd
parse_efficiency <- function(data,  norm, ref_genes, eff, eff_adjust) {
  if (!identical(norm, "reference")) {
    return(data[, .SD * log(2.0), by = ".group"])
  }
  stopifnot_(
    !is.null(eff) && is.list(eff),
    "Argument {.arg eff} must be a {.cls list} or a {.cls data.frame} object."
  )
  eff_names <- names(eff)
  stopifnot_(
    !is.null(eff_names),
    "Argument {.arg eff} must have names."
  )
  gene_cols <- setdiff(names(data), ".group")
  eff <- eff[gene_cols]
  stopifnot_(
    all(gene_cols %in% eff_names),
    c(
      "Efficiency values must be supplied for all genes:",
      `x` = "Column{?s} {.val {gene_cols[!gene_cols %in% eff_names]}}
            {?do/does} not have a corresponding efficiency value in
            argument {.arg eff}."
    )
  )
  eff <- unlist(eff)
  invalid <- isTRUE(eff < 0.0) | !is.finite(eff)
  stopifnot_(
    identical(length(invalid), 0L),
    c(
      "Primer efficiency values must be greater than 0.",
      `x` = "Non-positive or non-finite efficiency value{?s} {?was/were} found
             for gene{?s} {.val {gene_cols[invalid]}}"
    )
  )
  eff <- switch(
    eff_adjust,
    `none` = eff ,
    `limit` = pmin(eff, 1.0),
    `scale` = eff / max(eff)
  )
  data[, .SD * log(eff)[col(.SD)], by = ".group"]
}
