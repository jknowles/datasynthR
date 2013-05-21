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

myF <- list(vars=c("test1", "test2", "daysattended", "daysOUT", "bad"), 
                coefs=c(4, -2, 1.2, -7, 2))

studat$out <- genBinomialDV(studat, form=myF, intercept=7)

testGLM <- glm(out ~ test1 + test2 + daysattended + daysOUT, data=studat, 
               family="binomial")

summary(testGLM)

testGLM2 <- glm(out ~ test1 + test2 + daysattended + daysOUT + bad, data=studat, 
               family="binomial")

summary(testGLM2)

##################
# What about factor variables?

myFactors <- genFactor(25000, 6, nlevel=4, rho=0.3)
names(myFactors) <- c("econ", "race", "sped", "ell", "retention", "other")
studat <- cbind(studat, myFactors)
rm(myFactors)

testGLM3 <- glm(out ~ test1 + test2 + daysattended + daysOUT + bad + econ + 
                  race + sped + ell + retention + other, 
                data=studat, 
                family="binomial")

summary(testGLM3)

######################
# Chain together with factor variables

###########

seeds <- genNumeric(1000, 6, rho=0.3)

struc <- list(dist=c("norm", "norm", "unif", "pois", "pois", "gamma", 
                     "weibull"), 
              rho=c(0.7, 0.3, -0.5, 0.3, -0.8, 0.05, 0.7), 
              names=c("test1", "test2", "noise", "daysattended", 
                      "daysOUT", "bad", "bad2"), 
              seed = cbind(seeds[,1], seeds[,2], seeds[,3], seeds[, 4], seeds[, 5], 
                       seeds[, 6], seeds[,1]))

dat <- genNumeric(1000, pattern=struc)

cor(seeds[,1], dat[,1])
cor(seeds[,2], dat[,2])
cor(seeds[,3], dat[,3])
cor(seeds[,4], dat[,4])
cor(seeds[,5], dat[,5])
cor(seeds[,6], dat[,6])
cor(seeds[,1], dat[,7])

dat1 <- genFactor(1000, 3, nlevel=3, rho=0.8)

dat2 <- genFactor(1000, 4, nlevel=4, rho=0.3, seed=dat[,6])
dat3 <- genFactor(1000, 4, nlevel=6, rho=-0.7, seed=dat2[,4])
rm(dat3)

identical(dat2[,4], dat3[,1])

