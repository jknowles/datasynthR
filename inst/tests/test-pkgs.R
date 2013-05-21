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