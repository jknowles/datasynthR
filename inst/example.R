################################################################################
# Working example
################################################################################
source("R/correlatedsamplers.R")
source("R/generators.R")
source("R/missingness.R")
source("R/helpers.R")

struc <- list(dist=c("norm", "norm", "unif", "pois", "pois", "gamma", 
                        "weibull"), 
                 rho=c(0.7, 0.3, -0.5, 0.3, -0.8, 0.05, 0.7), 
                 names=c("test1", "test2", "noise", "daysattended", 
                         "daysOUT", "bad", "bad2"))


studat <- genNumeric(25000, pattern=struc)

cor(studat[, 1], studat[, 2])
cor(studat[, 1], studat[, 3])
cor(studat[, 1], studat[, 4])
cor(studat[, 1], studat[, 5])
cor(studat[, 1], studat[, 6])
cor(studat[, 1], studat[, 7])
cor(studat[, 1], studat[, 8])
