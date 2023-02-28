# Data table awareness
.datatable.aware <- TRUE

# For data.table and tidyselect
utils::globalVariables(c(".", ".I", ".N", ".SD", "where"))

#' Stop Function Execution Unless Condition Is True
#'
#' @inheritParams stop_
#' @param cond \[`logical(1)`] Condition to evaluate.
#' @noRd
stopifnot_ <- function(cond, message, ..., call = rlang::caller_env()) {
  if (!cond) {
    cli::cli_abort(message, ..., .envir = parent.frame(), call = call)
  }
}

#' Shorthand for `if (test) yes else no`
#'
#' @param test \[`logical(1)`] Condition to evaluate.
#' @param yes An \R object to return when `test` evaluates to `TRUE`.
#' @param no An \R object to return when `test` evaluates to `FALSE`.
#' @noRd
ifelse_ <- function(test, yes, no) {
  if (test) {
    return(yes)
  }
  no
}

#' Return `yes` If `test` Is `TRUE`, Otherwise Return `NULL`
#'
#' @param test \[`logical(1)`] Condition to evaluate.
#' @param yes An \R object to return when `test` evaluates to `TRUE`.
#' @noRd
onlyif <- function(test, yes) {
  if (test) {
    return(yes)
  }
  NULL
}

#' Number of Unique Values
#'
#' @inheritParams data.table::uniqueN
#' @noRd
n_unique <- data.table::uniqueN

local_methods <- c("randomization_test", "t_test", "wilcoxon_test")

global_methods <- c("anova", "kruskal_test")
