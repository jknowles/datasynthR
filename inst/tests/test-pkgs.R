library(testthat)
library(MASS)


takes_less_than <- function(amount) {
  function(expr) {
    duration <- system.time(force(expr))["elapsed"]
    
    expectation(
      duration < amount,
      stringr::str_c("took ", duration, " seconds, which is more than ", amount)
    )
  }
}

is_less_than <- function(expected, label=NULL, ...) {
  find_expr <- function(name, env = parent.frame()) {
    subs <- do.call("substitute", list(as.name(name), env))
    stringr::str_c(deparse(subs, width.cutoff = 500), collapse = "\n")
  }
  if (is.null(label)) {
    label <- find_expr("expected")
  }
  else if (!is.character(label) || length(label) != 1) {
    label <- deparse(label)
  }
  function(actual) {
    less <- expected > actual
    expectation(identical(less, TRUE), stringr::str_c("not less than ", 
                                                      label, "\n", 
                                                      stringr::str_c(less, collapse = "\n")))
  }
}

is_more_than <- function(expected, label=NULL, ...) {
  find_expr <- function(name, env = parent.frame()) {
    subs <- do.call("substitute", list(as.name(name), env))
    stringr::str_c(deparse(subs, width.cutoff = 500), collapse = "\n")
  }
  if (is.null(label)) {
    label <- find_expr("expected")
  }
  else if (!is.character(label) || length(label) != 1) {
    label <- deparse(label)
  }
  function(actual) {
    more <- expected < actual
    expectation(identical(more, TRUE), stringr::str_c("not more than ", 
                                                      label, "\n", 
                                                      stringr::str_c(more, collapse = "\n")))
  }
}

expect_unequal <- function (object, expected, ..., info = NULL, label = NULL, expected.label = NULL) 
{
  find_expr <- function(name, env = parent.frame()) {
    subs <- do.call("substitute", list(as.name(name), env))
    stringr::str_c(deparse(subs, width.cutoff = 500), collapse = "\n")
  }
  if (is.null(label)) {
    label <- find_expr("object")
  }
  if (is.null(expected.label)) {
    expected.label <- find_expr("expected")
  }
  expect_that(object, unequal(expected, label = expected.label, 
                              ...), info = info, label = label)
}


unequal <- function(expected, label = NULL, ...) 
{
  if (is.null(label)) {
    label <- find_expr("expected")
  }
  else if (!is.character(label) || length(label) != 1) {
    label <- deparse(label)
  }
  function(actual) {
    unequal <- !isTRUE(all.equal(expected, actual, ...))
    expectation(identical(unequal, TRUE), stringr::str_c("not equal to ", 
                                                         label, "\n", stringr::str_c(unequal, collapse = "\n")))
  }
}

#' Expectation: is returned value equal within a tolerance to some specified value?
#'
#' This is useful for testing whether a returned value is within a user specified distance 
#' from an expected value.
#
#' @param expected Expected value
#' @param tol tolerance for returned value to be distant from expected value
#' @param label For full form, label of expected object used in error
#'   messages. Useful to override default (deparsed expected expression) when
#'   doing tests in a loop.  For short cut form, object label. When
#'   \code{NULL}, computed from deparsed object.
#' @param expected.label Equivalent of \code{label} for shortcut form.
#' @param ... other values passed to \code{\link{all.equal}}
#' @family expectations
#' @export
#' @examples
#' a <- 11
#' expect_that(a, unequal(10))
#' expect_unequal(a, 10)
#' expect_that(a, approxto(10, tol=2)) # TRUE
#' expect_approxto(a, 10, tol=2)
approxto <- function(expected, tol, label=NULL, ...) {
  find_expr <- function(name, env = parent.frame()) {
    subs <- do.call("substitute", list(as.name(name), env))
    stringr::str_c(deparse(subs, width.cutoff = 500), collapse = "\n")
  }
  if (is.null(label)) {
    label <- find_expr("expected")
  }
  else if (!is.character(label) || length(label) != 1) {
    label <- deparse(label)
  }

  function(actual, tol=tol) {
#     expected <- expected
#    tol <- tol
     app <- isTRUE((expected - tol) < actual & actual < (expected + tol))
    expectation(identical(app, TRUE), 
                stringr::str_c("not within tolerance ", label, "\n", 
                               stringr::str_c(deparse(find_expr("actual")), 
                                              " is not ", 
                                              deparse(find_expr("expected")), 
                                              collapse = "\n")))
#      expectation(
#        (expected - tol) < actual & actual < (expected + tol),
#        stringr::str_c("not within tolerance ", label, "\n", 
#                       stringr::str_c(actual, 
#                                      " is not approximate", 
#                                       collapse = "\n"))#)
#      )
     
  }
}

#' @export
#' @rdname approxto
#' @inheritParams expect_that
expect_approxto <- function (object, expected, ..., info = NULL, label = NULL, expected.label = NULL) 
{
  find_expr <- function(name, env = parent.frame()) {
    subs <- do.call("substitute", list(as.name(name), env))
    stringr::str_c(deparse(subs, width.cutoff = 500), collapse = "\n")
  }
  if (is.null(label)) {
    label <- find_expr("object")
  }
  if (is.null(expected.label)) {
    expected.label <- find_expr("expected")
  }
  expect_that(object, approxto(expected, tol, label = expected.label, 
                               ...), info = info, label = label)
}
