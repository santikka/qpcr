#' qPCR Analysis
#'
#' @param data \[`data.frame`]\cr Data containing the Ct values and the
#'   treatment assignment. Should have one column per gene.
#' @param treatment \[`character(1)`]\cr Column name of `data` naming the
#'   column that defines the treatment group assignment. Should be `NULL` if
#'   there is only one treatment.
#' @param method \[`character(1)`]\cr Method to use for the analysis. Option
#'   `"normagene"` (the default) applies the NORMA-gene normalization method of
#'   Heckmann et al. to the data. Option `"reference"` uses reference genes
#'   instead, which should be defined in the argument `"ref_genes"`.
#' @param ref_genes \[`character()`]\cr A vector of column names of `data`
#'   which should be treated as reference genes when `method == "reference"`.
#'   This is argument is ignored for the `"normagene"` method.
qpcr <- function(data, treatment = NULL, method = c("normagene", "reference"),
                 ref_genes = NULL) {

}
