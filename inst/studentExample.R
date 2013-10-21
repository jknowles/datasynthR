################################################################################
#
# Simulate student data with 2 assessments
# Grade
# Attendance, discipline, mobility 
# Demographics 
# Assign them to schools
#
#  Yes / No treatment (.1 SD)
#  Dosage treatment (non-linear)
#
################################################################################

source("R/correlatedsamplers.R")
source("R/generators.R")
source("R/missingness.R")
source("R/helpers.R")

set.seed(12059)

N <- 100000

struc1 <- list(names=c("read", "math", "attend", "discipline", "mobility"),
              dist = c("norm", "norm", "gamma", "pois", "pois"), 
              rho = c(0.73, 0.8, 0.4, -0.5, -0.3))
# 
# struc2 <- list(names=c("read", "math", "attend", "discipline", "mobility"),
#                dist = c("norm", "norm", "gamma", "pois", "pois"), 
#                rho = c(0.73, 0.65, 0.3, -0.2, -0.1))


studat <- genNumeric(N, pattern=struc1)

# studat2 <- genNumeric(N, pattern=struc2, seed = studat[, 3])
# 


myFactors <- genFactor(N, 4, nlevel=4, rho=0.3)
names(myFactors) <- c("econ", "race", "sped", "ell")
studat <- cbind(studat, myFactors)


myF <- list(vars=c("read", "math", "attend", "discipline", "mobility"), 
            coefs=c(4, -2, 1.2, -7, 2))

studat$treat <- genBinomialDV(studat, form=myF, intercept=-7, type="binary")

schools <- genFactor(N, 1, 200, 0.1, seed=studat$seed, keepSeed=FALSE)
studat$schools <- schools
names(studat)[12] <- "schools"

studat$ever_treat <- studat$treat

fuzzy <- function(x) x + rnorm(1,0.1,0.1)

studat$read2 <- studat$read
studat$read2[studat$treat==1] <- sapply(studat$read2[studat$treat==1],fuzzy)

fuzzy <- function(x) x + rnorm(1,0.08,0.14)

studat$math2 <- studat$math
studat$math2[studat$treat==1] <- sapply(studat$math2[studat$treat==1],fuzzy)

studat$math <- studat$math2; studat$math2 <- NULL
studat$read <- studat$read2; studat$read2 <- NULL

studat$seed <- NULL

studat$attend2 <- 1 - studat$attend
studat$attend2 <- studat$attend2 + 0.25
studat$attend2[studat$attend2 > 1] <- 1
studat$attend2[studat$attend2 < 0] <- 0
studat$attend <- studat$attend2; studat$attend2 <- NULL


###############################################
# Attendance has too many 0 values, is 0 inflated

studat$discipline <- scale(studat$discipline)
studat$discipline[studat$discipline <0] <- 0
studat$discipline <- round(studat$discipline, digits=0)

###############################################
# Discipline variable needs more variability!

studat$mobility <- scale(studat$mobility)
studat$mobility[studat$mobility <0] <- 0
studat$mobility <- round(studat$mobility, digits=0)

studat$year <- 1

sids <- 100000:600000

studat$sid <- sample(sids, nrow(studat), replace=FALSE)
studat$year <- 2

tmp <- studat
tmp$year <- 1



################################################################################
## Add another year of data for all students
## pre-test math and reading scores affect post-treatment math and reading
## 3 years of data
## No year or grade effects (use normalized scores)


save(studat, file="simulated_treatment1.rda", compress="gzip")





