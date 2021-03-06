################################################################################
##
# Simulate student data with 2 assessments
# Grade
# Attendance, discipline, mobility 
# Demographics 
# Assign them to schools
#
#  Yes / No treatment (.1 SD)
#  Dosage treatment (non-linear)
#  
#  Simulate random assignment 
#  Biased assignment, observables
#  Biased assignment, unobservables
################################################################################

source("R/correlatedsamplers.R")
source("R/generators.R")
source("R/missingness.R")
source("R/helpers.R")

set.seed(12059)

N <- 100000

struc1 <- list(names=c("read", "math", "attend", "discipline", "mobility"),
              dist = c("norm", "norm", "gamma", "pois", "pois"), 
              rho = c(0.9, 0.75, 0.4, -0.7, 0.45))

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

studat$treat_pr <- genBinomialDV(studat, form=myF, intercept=-9, type="response")

weighted_binom <- function(x) rbinom(1, 1, x)

studat$treat <- sapply(studat$treat_pr, weighted_binom)
studat$treat_pr <- NULL

schools <- genFactor(N, 1, 200, 0.1, seed=studat$seed, keepSeed=FALSE)
studat$schools <- schools[,1]
names(studat)[ncol(studat)] <- "schools"

studat$ever_treat <- studat$treat

studat$seed <- NULL



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


sids <- 100000:600000

rm(schools)


studat$sid <- sample(sids, nrow(studat), replace=FALSE)
studat$year <- 2
row.names(studat) <- NULL

mtmp <- studat
mtmp$year <- 1
mtmp$treat <- 0
row.names(mtmp) <- 100001:200000

#studat <- rbind.fill(studat, mtmp)

studat <- rbind(studat, mtmp)

mtmp$year <- 3
mtmp$treat <- 0

studat <- rbind(studat, mtmp)

rm(mtmp, myF, sids, N, struc1)

################################################################################
# Add growth
################################################################################

growthG <- function(x) x + runif(1, -0.9, 0.9)

studat$attend2 <- studat$attend
studat$attend2[studat$year==2] <- sapply(studat$attend2[studat$year==2], growthG)
studat$attend2[studat$year==3] <- sapply(studat$attend2[studat$year==3], growthG)

growth <- function(x) x + rnorm(1, 0.1, 1)


studat$read2 <- studat$read
studat$read2[studat$year==2] <- sapply(studat$read2[studat$year==2], growth)
studat$read2[studat$year==3] <- sapply(studat$read2[studat$year==3], growth)


studat$math2 <- studat$math
studat$math2[studat$year==2] <- sapply(studat$math2[studat$year==2], growth)
studat$math2[studat$year==3] <- sapply(studat$math2[studat$year==3], growth)

growth2 <- function(x) rpois(1, x)

studat$discipline2 <- studat$discipline
studat$discipline2[studat$year==2] <- sapply(studat$discipline2[studat$year==2], growth2)
studat$discipline2[studat$year==3] <- sapply(studat$discipline2[studat$year==3], growth2)

studat$mobility2 <- studat$mobility
studat$mobility2[studat$year==2] <- sapply(studat$mobility2[studat$year==2], growth2)
studat$mobility2[studat$year==3] <- sapply(studat$mobility2[studat$year==3], growth2)


################################################################################
# Treatment

fuzzy <- function(x) x + rnorm(1,0.1,0.1)

studat$read2 <- studat$read
studat$read2[studat$ever_treat==1 & studat$year==3] <- 
  sapply(studat$read2[studat$ever_treat==1 & studat$year==3],fuzzy)

fuzzy <- function(x) x + rnorm(1,0.08,0.14)

studat$math2 <- studat$math
studat$math2[studat$ever_treat==1 & studat$year==3] <- 
  sapply(studat$math2[studat$ever_treat==1 & studat$year==3],fuzzy)

studat$math <- studat$math2; studat$math2 <- NULL
studat$read <- studat$read2; studat$read2 <- NULL
studat$attend <- studat$attend2; studat$attend2 <- NULL
studat$mobility <- studat$mobility2; studat$mobility2 <- NULL
studat$discipline <- studat$discipline2; studat$discipline2 <- NULL

studat$attend <- studat$attend + 0.3
studat$attend[studat$attend < 0] <- 0
studat$attend[studat$attend > 1] <- 1


save(studat, file="simulated_treatment_propensity.rda", compress="xz")

################################################################################
################################################################################
# Create random assignment experiment
################################################################################

set.seed(784)

N <- 100000

struc1 <- list(names=c("read", "math", "attend", "discipline", "mobility"),
               dist = c("norm", "norm", "gamma", "pois", "pois"), 
               rho = c(0.9, 0.75, 0.4, -0.7, 0.45))

# struc2 <- list(names=c("read", "math", "attend", "discipline", "mobility"),
#                dist = c("norm", "norm", "gamma", "pois", "pois"), 
#                rho = c(0.73, 0.65, 0.3, -0.2, -0.1))


studat <- genNumeric(N, pattern=struc1)

# studat2 <- genNumeric(N, pattern=struc2, seed = studat[, 3])
# 


myFactors <- genFactor(N, 4, nlevel=4, rho=0.3)
names(myFactors) <- c("econ", "race", "sped", "ell")
studat <- cbind(studat, myFactors)

schools <- genFactor(N, 1, 200, 0.1, seed=studat$seed, keepSeed=FALSE)
studat$schools <- schools[,1]
names(studat)[ncol(studat)] <- "schools"

studat$seed <- NULL



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

rm(schools)


studat$sid <- sample(sids, nrow(studat), replace=FALSE)
studat$year <- 2
row.names(studat) <- NULL

treats <- sample(studat$sid, 10000, replace=FALSE)

studat$treat <- 0
studat$treat[studat$sid %in% treats] <- 1

studat$ever_treat <- studat$treat


mtmp <- studat
mtmp$year <- 1
mtmp$treat <- 0
row.names(mtmp) <- 100001:200000

#studat <- rbind.fill(studat, mtmp)

studat <- rbind(studat, mtmp)

mtmp$year <- 3
mtmp$treat <- 0

studat <- rbind(studat, mtmp)

rm(mtmp, myF, sids, N, struc1)

################################################################################
# Add growth
################################################################################

growthG <- function(x) x + runif(1, -0.9, 0.9)

studat$attend2 <- studat$attend
studat$attend2[studat$year==2] <- sapply(studat$attend2[studat$year==2], growthG)
studat$attend2[studat$year==3] <- sapply(studat$attend2[studat$year==3], growthG)

growth <- function(x) x + rnorm(1, 0.1, 1)


studat$read2 <- studat$read
studat$read2[studat$year==2] <- sapply(studat$read2[studat$year==2], growth)
studat$read2[studat$year==3] <- sapply(studat$read2[studat$year==3], growth)


studat$math2 <- studat$math
studat$math2[studat$year==2] <- sapply(studat$math2[studat$year==2], growth)
studat$math2[studat$year==3] <- sapply(studat$math2[studat$year==3], growth)

growth2 <- function(x) rpois(1, x)

studat$discipline2 <- studat$discipline
studat$discipline2[studat$year==2] <- sapply(studat$discipline2[studat$year==2], growth2)
studat$discipline2[studat$year==3] <- sapply(studat$discipline2[studat$year==3], growth2)

studat$mobility2 <- studat$mobility
studat$mobility2[studat$year==2] <- sapply(studat$mobility2[studat$year==2], growth2)
studat$mobility2[studat$year==3] <- sapply(studat$mobility2[studat$year==3], growth2)


################################################################################
# Treatment

fuzzy <- function(x) x + rnorm(1,0.1,0.1)

studat$read2 <- studat$read
studat$read2[studat$ever_treat==1 & studat$year==3] <- 
  sapply(studat$read2[studat$ever_treat==1 & studat$year==3],fuzzy)

fuzzy <- function(x) x + rnorm(1,0.08,0.14)

studat$math2 <- studat$math
studat$math2[studat$ever_treat==1 & studat$year==3] <- 
  sapply(studat$math2[studat$ever_treat==1 & studat$year==3],fuzzy)

studat$math <- studat$math2; studat$math2 <- NULL
studat$read <- studat$read2; studat$read2 <- NULL
studat$attend <- studat$attend2; studat$attend2 <- NULL
studat$mobility <- studat$mobility2; studat$mobility2 <- NULL
studat$discipline <- studat$discipline2; studat$discipline2 <- NULL

studat$attend <- studat$attend + 0.3
studat$attend[studat$attend < 0] <- 0
studat$attend[studat$attend > 1] <- 1


save(studat, file="simulated_treatment_random.rda", compress="xz")


################################################################################
################################################################################
##
## Biased treatment assignment, unobservable
##
################################################################################
################################################################################
set.seed(487)

N <- 100000

struc1 <- list(names=c("read", "math", "luck", "attend", "discipline", "mobility"),
               dist = c("norm", "norm", "norm", "gamma", "pois", "pois"), 
               rho = c(0.9, 0.75, 0.5, 0.6, -0.7, 0.45))

# struc2 <- list(names=c("read", "math", "attend", "discipline", "mobility"),
#                dist = c("norm", "norm", "gamma", "pois", "pois"), 
#                rho = c(0.73, 0.65, 0.3, -0.2, -0.1))


studat <- genNumeric(N, pattern=struc1)

# studat2 <- genNumeric(N, pattern=struc2, seed = studat[, 3])
# 


myFactors <- genFactor(N, 4, nlevel=4, rho=0.3)
names(myFactors) <- c("econ", "race", "sped", "ell")
studat <- cbind(studat, myFactors)


myF <- list(vars=c("read", "math", "luck", "discipline", "mobility"), 
            coefs=c(0.1, 0.1, 0.5, -0.2, -0.2))

studat$treat_pr <- genBinomialDV(studat, form=myF, intercept=-9, type="response")

weighted_binom <- function(x) rbinom(1, 1, x)
studat$treat <- sapply(studat$treat_pr, weighted_binom)
studat$treat_pr <- NULL

schools <- genFactor(N, 1, 200, 0.1, seed=studat$seed, keepSeed=FALSE)
studat$schools <- schools[,1]
names(studat)[ncol(studat)] <- "schools"

studat$ever_treat <- studat$treat


studat$seed <- NULL

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

rm(schools)


studat$sid <- sample(sids, nrow(studat), replace=FALSE)
studat$year <- 2
row.names(studat) <- NULL

treats <- sample(studat$sid, 10000, replace=FALSE)

studat$treat <- 0
studat$treat[studat$sid %in% treats] <- 1

studat$ever_treat <- studat$treat


mtmp <- studat
mtmp$year <- 1
mtmp$treat <- 0
row.names(mtmp) <- 100001:200000

#studat <- rbind.fill(studat, mtmp)

studat <- rbind(studat, mtmp)

mtmp$year <- 3
mtmp$treat <- 0

studat <- rbind(studat, mtmp)

rm(mtmp, myF, sids, N, struc1)

################################################################################
# Add growth
################################################################################

growthG <- function(x) x + runif(1, -0.9, 0.9)

studat$attend2 <- studat$attend
studat$attend2[studat$year==2] <- sapply(studat$attend2[studat$year==2], growthG)
studat$attend2[studat$year==3] <- sapply(studat$attend2[studat$year==3], growthG)

growth <- function(x) x + rnorm(1, 0.1, 1)


studat$read2 <- studat$read
studat$read2[studat$year==2] <- sapply(studat$read2[studat$year==2], growth)
studat$read2[studat$year==3] <- sapply(studat$read2[studat$year==3], growth)


studat$math2 <- studat$math
studat$math2[studat$year==2] <- sapply(studat$math2[studat$year==2], growth)
studat$math2[studat$year==3] <- sapply(studat$math2[studat$year==3], growth)

growth2 <- function(x) rpois(1, x)

studat$discipline2 <- studat$discipline
studat$discipline2[studat$year==2] <- sapply(studat$discipline2[studat$year==2], growth2)
studat$discipline2[studat$year==3] <- sapply(studat$discipline2[studat$year==3], growth2)

studat$mobility2 <- studat$mobility
studat$mobility2[studat$year==2] <- sapply(studat$mobility2[studat$year==2], growth2)
studat$mobility2[studat$year==3] <- sapply(studat$mobility2[studat$year==3], growth2)


################################################################################
# Treatment

fuzzy <- function(x) x + rnorm(1,0.1,0.1)

studat$read2 <- studat$read
studat$read2[studat$ever_treat==1 & studat$year==3] <- 
  sapply(studat$read2[studat$ever_treat==1 & studat$year==3],fuzzy)

fuzzy <- function(x) x + rnorm(1,0.08,0.14)

studat$math2 <- studat$math
studat$math2[studat$ever_treat==1 & studat$year==3] <- 
  sapply(studat$math2[studat$ever_treat==1 & studat$year==3],fuzzy)

studat$math <- studat$math2; studat$math2 <- NULL
studat$read <- studat$read2; studat$read2 <- NULL
studat$attend <- studat$attend2; studat$attend2 <- NULL
studat$mobility <- studat$mobility2; studat$mobility2 <- NULL
studat$discipline <- studat$discipline2; studat$discipline2 <- NULL

studat$attend <- studat$attend + 0.3
studat$attend[studat$attend < 0] <- 0
studat$attend[studat$attend > 1] <- 1


save(studat, file="simulated_treatment_biased.rda", compress="xz")
