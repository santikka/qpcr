#' qPCR Analysis
#'
#' @export
#' @param data \[`data.frame` or `data.table`]\cr Data containing the
#'   Ct values, group specification, and optionally population specification.
#'   Should have one column per gene.
#' @param group \[`character(1)`]\cr Column names of `data` whose combination
#'   defines the unique groups. Columns that are not listed in this argument
#'   are assumed to contain Ct values of genes to be analyzed (or to be used
#'   as reference genes).
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
#'   randomization test, student t-test and the Wilcoxon rank-sum test,
#'   respectively. See the package vignette for more information on the
#'   randomization test. See [stats::t.test], and [stats::wilcox.test] for more
#'   information on the other tests. Comparisons between multiple
#'   groups can be carried out via options `anova` and `kruskal_test`,
#'   which apply the analysis of variance and the Kruskal-Wallis
#'   rank sum test to the data, respectively. See [stats::kruskal.test] for
#'   more information on the test. By default, only the randomization test
#'   is applied. If `NULL`, applies all methods.
#' @param comparisons \[`character()`]\cr An optional vector describing the
#'   comparisons to carry out between groups. For example, suppose that
#'   `groups` consists of two columns: `treatment` and `population` (in this
#'   order). Then each element should follow the syntax:
#'   `"geneA:treatmentA:populationA/geneB:groupB:populationB"` which defines
#'   a comparison of `geneA` for `treatmentA` in `populationA` against
#'   `geneB` for `treatmentB` in `populationB`. Individual terms can be
#'   omitted, which corresponds to a wildcard of all possible comparisons
#'   instead, for example: `geneA::/geneA::` would compare `geneA` between all
#'   combinations of the treatments and populations. Similarly,
#'   `"geneA::/geneA:control:` would create all comparisons against a
#'   control group (assuming `"control"` is a level of `treatment` in `data`),
#'   between all populations for `geneA`. If `comparisons` is `NULL`, all
#'   possible pairs are compared (i.e., the same as `"::/::"`). There should
#'   be the same number of semicolons as there are variables in `group`,
#'   e.g., say we only had the `treatment` column in the previous example, then
#'   the elements of `comparisons` should omit the last colon, e.g.,
#'   `geneA:/geneA:control`.
#' @param ref_genes \[`character()`]\cr A vector of column names of `data`
#'   which should be treated as reference genes when `norm == "reference"`.
#'   This is argument is ignored for the `"normagene"` and `"none"`
#'   normalization options.
#' @param eff \[`data.frame` or `list`]\cr Efficiencies when using the
#'   reference gene normalization. Should be either a `data.frame` with a
#'   column for each gene and a single row providing the efficiency values.
#'   Alternatively, a named `list` or a `vector` giving the efficiency values
#'   can be provided.
#' @param eff_adjust \[`character(1)`]\cr How should efficiencies above 1 be
#'   adjusted? The default `"none"` uses the efficiencies as is. Option
#'   `"limit"` caps efficiencies to a maximum of 1. Option `"scale"` adjust all
#'   efficiency values such that the highest is 1.
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
                 comparisons = NULL, ref_genes = NULL, eff = NULL,
                 eff_adjust = c("none", "limit", "scale"),
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
  eff_adjust <- try(
    match.arg(eff_adjust, c("none", "limit", "scale")),
    silent = TRUE
  )
  stopifnot_(
    !inherits(eff_adjust, "try-error"),
    "Argument {.arg eff_adjust} must be either
    {.val none}, {.val limit}, or {.val scale}."
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
  data <- parse_data(data, group, norm, ref_genes, eff, eff_adjust)
  comparisons <- parse_comparisons(data, group, ref_genes, comparisons)
  qpcr_(data, group, norm, methods, comparisons, ref_genes, m, alpha)
}

#' qPCR Analysis
#'
#' @inheritParams qpcr
#' @noRd
qpcr_ <- function(data, group, norm,
                  methods, comparisons, ref_genes, m, alpha) {
  local <- NULL
  global <- NULL
  args <- as.list(match.call()[-1L])
  if (any(local_methods  %in% methods)) {
    local <- do.call("local_tests", args)
  }
  #if (any(global_methods %in% methods)) {
  #  global <- do.call("global_tests", args)
  #}
  plots <- NULL
  #plots <- fold_change_plots(data, ref_genes)
  list(local = local, global = global, plots = plots)
}
