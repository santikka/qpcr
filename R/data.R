#' Parse Input Data For Analysis
#'
#' @inheritParams qpcr
#' @noRd
parse_data <- function(data, treatment, norm, genes, ref_genes) {
  stopifnot_(
    is.data.frame(data),
    "Argument {.arg data} must be a {.cls data.frame} object."
  )
  stopifnot_(
    checkmate::test_string(x = treatment),
    "Argument {.arg treatment} must be a single character string."
  )
  stopifnot_(
    !is.null(data[[treatment]]),
    c(
      "Argument {.arg treatment} must be a column name of {.arg data}:",
      `x` = "Column {.var {treatment}} does not exist in {.arg data}."
    )
  )
  data <- data.table::as.data.table(data.table::copy(data))
  data.table::set(
    x = data,
    j = treatment,
    value = as.factor(data[[treatment]])
  )
  cols <- names(data)
  stopifnot_(
    is.null(genes) || all(genes %in% cols),
    c(
      "Argument {.arg genes} contains unknown columns:",
      `x` = "Column{?s} {.val {genes[!genes %in% cols]}}
             were not found in {.arg data}."
    )
  )
  stopifnot_(
    is.null(genes) || !treatment %in% genes,
    "Treatment assignment column {.val {treatment}}
     is also defined in {.arg genes}."
  )
  stopifnot_(
    is.null(ref_genes) || checkmate::test_character(x = ref_genes),
    "Argument {.arg ref_genes} must be a {.cls character} vector."
  )
  stopifnot_(
    is.null(ref_genes) || all(ref_genes %in% cols),
    c(
      "Argument {.arg ref_genes} contains unknown columns:",
      `x` = "Column{?s} {.val {ref_genes[!ref_genes %in% cols]}}
             {?was/were} not found in {.arg data}."
    )
  )
  stopifnot_(
    is.null(ref_genes) || !treatment %in% genes,
    "Treatment assignment column {.val {treatment}}
     is also defined in {.arg ref_genes}."
  )
  stopifnot_(
    is.null(genes) || is.null(ref_genes) || !any(ref_genes %in% genes),
    c(
      "Arguments {.arg genes} and {.arg ref_genes} must not overlap:",
      `x` = "Column{?s} {.val {ref_genes}} {?was/were}
             found in both arguments."
    )
  )
  genes <- ifelse_(
    is.null(genes),
    setdiff(cols, union(treatment, ref_genes)),
    genes
  )
  data <- data[, .SD, .SDcols = c(genes, ref_genes, treatment)]
  gene_cols <- c(genes, ref_genes)
  col_types <- vapply(
    data[, .SD, .SDcols = gene_cols],
    typeof,
    character(1L),
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
}

#' Parse Efficiency Values For Analysis
#'
#' @inheritParams qpcr
#' @noRd
parse_efficiency <- function(data, treatment, norm, ref_genes, eff) {
  if (!identical(norm, "reference")) {
    return()
  }
  stopifnot_(
    !is.null(eff) && is.list(eff),
    "Argument {.arg eff} must be a {.cls list} or a {.cls data.frame} object."
  )
  eff_names <- names(eff)
  stopifnot(
    !is.null(names(eff_names)),
    "Argument {.arg eff} must have names."
  )
  gene_cols <- setdiff(names(data), treatment)
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
}
