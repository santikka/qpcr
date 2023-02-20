#' The `qpcr` package.
#'
#' @description A package for qPCR analysis
#'
#' # Authors
#' | Santtu Tikka (author) | <santtuth@@gmail.com>  |
#' | --------------------- | ---------------------- |
#' | Aigi Margus (author)  | <aigi.margus@@jyu.fi>  |
#'
#' @docType package
#' @name qpcr-package
#' @importFrom data.table :=
#' @importFrom loo loo
NULL

# Data table awareness
.datatable.aware <- TRUE

# For data.table and tidyselect
utils::globalVariables(c(".", ".I", ".N", ".SD", "where"))
