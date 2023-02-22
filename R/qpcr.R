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
#' @param genes \[`character()`]\cr An optional vector of column names of
#'   `data` naming the genes that should be included in the analysis.
#'   If `NULL`, it is assumed that all columns expect those named in
#'   `treatment` and `ref_genes` arguments should be analyzed.
#' @param ref_genes \[`character()`]\cr A vector of column names of `data`
#'   which should be treated as reference genes when `nor == "reference"`.
#'   This is argument is ignored for the `"normagene"` and `"none"`
#'   normalization options.
#' @param eff \[`data.frame` or `list`]\cr Efficiencies when using the
#'   reference gene normalization. Should be either a `data.frame` with a
#'   column for each gene and a single row providing the efficiency values.
#'   Alternatively, a named `list` or a `vector` giving the efficiency values
#'   can be provided.
#' @param m \[`integer(1)`]\cr Number of replications to use for the
#'   randomization test when using the method. The default is `10000`
#' @examples
#' # TODO examples
#'
qpcr <- function(data, treatment, norm = c("normagene", "reference", "none"),
                 methods = c("rand", "t", "anova", "wilcoxon", "kw", "mw"),
                 genes = NULL, ref_genes = NULL, eff = NULL, m = 10000L) {
  stopifnot_(
    !missing(data),
    "Argument {.arg data} is missing."
  )
  stopifnot_(
    !missing(treatment),
    "Argument {.arg treatment} is missing."
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
  stopifnot_(
    checkmate::test_int(x = m),
    "Argument {.arg m} must be a single {.cls integer} value."
  )
  data <- parse_data(data, treatment, norm, genes)
  ref_genes <- onlyif(identical(norm, "reference"), ref_genes)
  eff <- onlyif(identical(norm, "reference"), unlist(eff))
  parse_efficiency(data, treatment, norm, ref_genes, eff)
  qpcr_(data, treatment, norm, methods, ref_genes, eff, m)
}

#' qPCR Analysis
#'
#' @inheritParams qpcr
#' @noRd
qpcr_ <- function(data, treatment, norm, methods, ref_genes, eff) {
  if ("rand" %in% methods) {
    do.call(
      what = paste0("randomization_test_", norm),
      args = list(
        data = data,
        treatment = treatment,
        ref_genes = ref_genes,
        eff = eff,
        m = m
      )
    )
  }
  data <- do.call(
    what = paste0("normalize_", norm),
    args = list(
      data = data,
      treatment = treatment,
      ref_genes = ref_genes,
      eff = eff
    )
  )
}
