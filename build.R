library(roxygen2)

roclet <- rd_roclet()
roc_out(roclet, "R/generators.R", ".")
