library(testthat)
source("inst/test-pkgs.R")
source("R/correlatedsamplers.R")
source("R/generators.R")

test_file("inst/test-generators.R")
test_file("inst/test-genFactor.R")
test_file("inst/test-correlatedsamplers.R")