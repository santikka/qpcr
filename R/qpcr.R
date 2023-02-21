#' qPCR Analysis
#'
#' @param data \[`data.frame` or `data.table`]\cr Data containing the
#'   Ct values and the treatment assignment. Should have one column per gene.
#' @param treatment \[`character(1)`]\cr Column name of `data` naming the
#'   column that defines the treatment group assignment.
#' @param norm \[`character(1)`]\cr Normalization method to use for the
#'   analysis. Option `"normagene"` (the default) applies the NORMA-gene
#'   normalization method of Heckmann et al. to the data. Option `"reference"`
#'   uses reference gene normalization instead, and the reference genes should
#'   be defined in the argument `ref_genes` along with the efficiencies in
#'   the argument `eff`.  Option `"none"` uses the data as is without applying
#'   any normalization.
#' @param methods \[`character()`]\cR A vector describing the names of the
#'   statistical methods to apply to the data. The default option `"rand"`
#'   applies a randomization test to the data for each gene and treatment
#'   level pair. Other possible options to include are `"t"`, `"anova"`,
#'   `"wilcoxon"`, `"kw"`, and `"mw"` for the t-test, analysis of variance,
#'   Wilcoxon signed rank test, Kruskal-Wallis H test, and Mann-Whitney U test,
#'   respectively.
#' @param ref_genes \[`character()`]\cr A vector of column names of `data`
#'   which should be treated as reference genes when `nor == "reference"`.
#'   This is argument is ignored for the `"normagene"` and `"none"`
#'   normalization options.
#' @param eff \[`data.frame`]\cr Efficiencies as a `data.frame` or
#'   a named `list` for each gene when when `norm == "reference"`.
#'   This is argument is ignored for the `"normagene"` and `"none"`
#'   normalization options.
#' @examples
#' # TODO examples
#'
qpcr <- function(data, treatment, norm = c("normagene", "reference", "none"),
                 methods = c("rand", "t", "anova", "wilcoxon", "kw", "mw"),
                 ref_genes = NULL, eff = NULL) {
  stopifnot_(
    !missing(data),
    "Argument {.arg data} is missing."
  )
  stopifnot_(
    is.data.frame(data),
    "Argument {.arg data} must be a {.cls data.frame} object."
  )
  stopifnot_(
    !missing(treatment),
    "Argument {.arg treatment} is missing."
  )
  stopifnot_(
    checkmate::test_string(x = treatment),
    "Argument {.arg treatment} must be a single character string."
  )
  stopifnot_(
    !is.null(data[[treatment]]),
    c(
      "Argument {.arg treatment} must be a column name of {.arg data}.",
      `x` = "Column {.var {treatment}} does not exist in {.arg data}."
    )
  )
  norm <- onlyif(is.character(norm), tolower(norm))
  norm <- try(
    match.arg(norm, c("normagene", "reference", "none")),
    silent = TRUE
  )
  stopifnot_(
    !inherits(norm, "try-error"),
    "Argument {.arg norm} must be either
    {.val normagene}, {.val reference}, or {.val none}."
  )
  methods <- onlyif(is.character(methods), tolower(methods))
  methods <- try(
    match.arg(methods, c("rand", "t", "anova", "wilcoxon", "kw", "mw"), TRUE),
    silent = TRUE
  )
  stopifnot_(
    !inherits(methods, "try-error"),
    "Argument {.arg methods} must contain one or more of the following:
    {.val rand}, {.val t}, {.val anova}, {.val wilcoxon}, {.val kw}, {.val mw}"
  )
  data <- parse_data(data, treatment)
  if (identical(norm, "normagene")) {
    data <- normagene(data, treatment)
  }
  if (identical(norm, "reference")) {
    stopifnot_(
      !is.null(ref_genes) && checkmate::test_character(x = ref_genes),
      "Argument {.arg ref_genes} must be a {.cls characte}r vector."
    )
    stopifnot_(
      !is.null(eff) && is.list(eff),
      "Argument {.arg eff} must be a {.cls list} or a {.cls data.frame}."
    )
    qpcr_reference(data, treatment, methods, ref_genes, eff)
  } else {
    qpcr_default(data, treatment, methods)
  }
}

#' qPCR Analysis Using Reference Genes

#' qPCR Analysis Without Reference Genes
#'
#' @inheritParams qpcr
#' @noRd
qpcr_default <- function(data, treatment, methods) {
  # TODO
}
