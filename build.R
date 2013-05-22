library(roxygen2)

roclet <- rd_roclet()
roc_out(roclet, "R/generators.R", ".")
roc_out(roclet, "R/correlatedsamplers.R", ".")
roc_out(roclet, "R/helpers.R", ".")
roc_out(roclet, "R/missingness.R", ".")



