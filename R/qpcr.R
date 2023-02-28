#' qPCR Analysis
#'
#' @export
#' @param data \[`data.frame` or `data.table`]\cr Data containing the
#'   Ct values and the group definitions. Should have one column per gene.
#' @param group \[`character(1)`]\cr Column name of `data` that defines the
#'   groups.
#' @param norm \[`character(1)`]\cr Normalization method to use for the
#'   analysis. Option `"normagene"` (the default) applies the NORMA-gene
#'   normalization method of Heckmann et al. to the data. Option `"reference"`
#'   uses reference gene normalization instead, and the reference genes should
#'   be defined in the argument `ref_genes` along with the efficiencies in
#'   the argument `eff`.  Option `"none"` uses the data as is without applying
#'   any normalization.
#' @param methods \[`character()`]\cr A vector describing the names of the
#'   statistical methods to apply to the data. Pairwise comparisons between
#'   groups for each gene can be carried out via options
#'   `"randomization_test"`, `"t_test"`, and `"wilcoxon_test"`, which apply a
#'   randomization test, student t-test and the wilcoxon rank-sum test,
#'   respectively. See the package vignette for more information on the
#'   randomization test. See [stats::t.test], and [stats::wilcox.test] for more
#'   information on the other tests. Comparisons between multiple
#'   groups can be carried out via options `anova` and `kruskal_test`,
#'   which apply the analysis of variance and the Kruskal-Wallis
#'   rank sum test to the data, respectively. See [stats::kruskal.test] for
#'   more information on the test. By default, only the randomization test
#'   is applied. If `NULL`, applies all methods.
#' @param genes \[`character()`]\cr An optional vector of column names of
#'   `data` naming the genes that should be included in the analysis.
#'   If `NULL`, it is assumed that all columns expect those named in
#'   `group` and `ref_genes` arguments should be analyzed.
#' @param ref_genes \[`character()`]\cr A vector of column names of `data`
#'   which should be treated as reference genes when `norm == "reference"`.
#'   This is argument is ignored for the `"normagene"` and `"none"`
#'   normalization options.
#' @param eff \[`data.frame` or `list`]\cr Efficiencies when using the
#'   reference gene normalization. Should be either a `data.frame` with a
#'   column for each gene and a single row providing the efficiency values.
#'   Alternatively, a named `list` or a `vector` giving the efficiency values
#'   can be provided.
#' @param m \[`integer(1)`]\cr Number of replications to use for the
#'   randomization test. The default is `10000` replications.
#' @param alpha \[`numeric(1)`]\cr Confidence level to use for the
#'   confidence intervals (CI). The default is `0.95` for 95% CIs.
#' @examples
#' # TODO examples
#'
qpcr <- function(data, group, norm = c("normagene", "reference", "none"),
                 methods = c(
                   "randomization_test", "t_test", "wilcoxon_test",
                   "anova", "kruskal_test"
                  ),
                 genes = NULL, ref_genes = NULL, eff = NULL,
                 m = 10000L, alpha = 0.95) {
  stopifnot_(
    !missing(data),
    "Argument {.arg data} is missing."
  )
  stopifnot_(
    !missing(group),
    "Argument {.arg group} is missing."
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
  all_methods <- c(local_methods, global_methods)
  methods <- onlyif(is.character(methods), tolower(methods))
  methods <- ifelse_(
    is.null(methods) || length(methods) == 0L,
    all_methods,
    methods
  )
  methods <- try(
    match.arg(methods, all_methods, several.ok = TRUE),
    silent = TRUE
  )
  stopifnot_(
    !inherits(methods, "try-error"),
    "Argument {.arg methods} must contain one or more of the following:
    {.val {all_methods}}"
  )
  stopifnot_(
    checkmate::test_int(x = m, lower = 1L),
    "Argument {.arg m} must be a single positive {.cls integer} value."
  )
  stopifnot_(
    checkmate::test_double(x = alpha, lower = 0.0, upper = 1.0),
    "Argument {.arg alpha} must be a single
     {.cls numeric} value between 0 and 1."
  )
  ref_genes <- onlyif(identical(norm, "reference"), ref_genes)
  eff <- onlyif(identical(norm, "reference"), unlist(eff))
  data <- parse_data(data, group, norm, genes, ref_genes, eff)
  qpcr_(data, group, norm, methods, ref_genes, m, alpha)
}

#' qPCR Analysis
#'
#' @inheritParams qpcr
#' @noRd
qpcr_ <- function(data, group, norm, methods, ref_genes, m, alpha) {
  local <- NULL
  global <- NULL
  args <- as.list(match.call()[-1L])
  if (any(local_methods  %in% methods)) {
    local <- do.call("local_tests", args)
  }
  if (any(global_methods %in% methods)) {
    global <- do.call("global_tests", args)
  }
  list(local = local, global = global)
}
