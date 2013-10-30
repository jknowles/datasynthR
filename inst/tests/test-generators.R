################################################################################
# Test generators
################################################################################

context("Simple Random Correlation Matrices")


N <- 8000
K <- 10
RHO1 <- 0.45
covmat <- genNumeric(N, K, rho=RHO1)

test_that("Dimensions are correct", {
  expect_equal(nrow(covmat), N)
  expect_equal(ncol(covmat), K)
  expect_is(covmat, "matrix")
  
})


covE <- rep(NA, K-1)

for(i in 2:ncol(covmat)){
  covE[i-1] <- cor(covmat)[i-1, i]
}

tol <- RHO1 * K * .005
err <- sum(abs(covE - RHO1)) / K

test_that("Correlations are correct", {
  expect_that(err, is_less_than(tol))
})



context("User Controlled Random Correlation Matrices")


struc <- list(dist = c("norm", "norm", "pois"), 
                rho = c(-.05, -.4, 0.4),
              names = c("score", "score2", "accept"))


N <- 5000

covmat <- genNumeric(N, pattern=struc)

test_that("Dimensions are correct in simple case", {
  expect_equal(nrow(covmat), N)
  expect_equal(ncol(covmat), length(struc$dist))
  expect_is(covmat, "data.frame")
  
})

err <- sum(abs(cor(covmat)[2:ncol(covmat),1] - struc$rho[2:length(struc$rho)])) / length(struc$rho)
tol <- sum(abs(struc$rho[2:length(struc$rho)]) *.05)

test_that("Correlations are correct in simple case", {
  expect_that(err, is_less_than(tol))
})

test_that("Names get passed properly in simple case", {
  expect_equivalent(names(covmat), struc$names)
  expect_is(covmat, "data.frame")
  
})

context("User Controlled Random Correlation Matrices Complex")

struc2 <- list(dist = c("norm", "chisq", "pois", "norm", 
                       "weibull", "gamma"), 
              rho = c(-.05, -.4, 0.3, 0.9, .03, -.6),
              names = c("score", "accept", "score2", "days", "days2", 
                        "luck"))

covmat <- genNumeric(N, pattern=struc2)



test_that("Dimensions are correct in complex case", {
  expect_equal(nrow(covmat), N)
  expect_equal(ncol(covmat), length(struc2$dist))
  expect_is(covmat, "data.frame")
  
})

err <- sum(abs(cor(covmat)[2:ncol(covmat),1] - struc$rho[2:length(struc2$rho)])) / length(struc2$rho)
tol <- sum(abs(struc2$rho[2:length(struc2$rho)]) *.05)

test_that("Correlations are correct in complex case", {
  expect_that(err, is_less_than(tol))
})

test_that("Names get passed properly in complex case", {
  expect_equivalent(names(covmat), struc2$names)
  expect_is(covmat, "data.frame")
  
})

              
context("User specified Seed is Correct")

N <- 1000
K <- 6
RHO1 <- 0.3
RHO2 <- -0.9
S1 <- rnorm(N)
S2 <- runif(N)



struc3 <- list(dist = c("norm", "chisq", "pois", "norm", 
                        "weibull", "gamma"), 
               rho = c(-.05, -.4, 0.3, 0.9, .03, -.6),
               names = c("score", "accept", "score2", "days", "days2", 
                         "luck"),
               seed = c(runif(N), rpois(N, 7), rpois(N, 3), rgamma(N, shape=2), 
                        runif(N)))


covmat <- genNumeric(N, K, rho=RHO1, seed=S1)
covmat2 <- genNumeric(N, K, rho=RHO2, seed=S2)
covmat3 <- genNumeric(N, pattern = struc3)

err1 <- cor(covmat)[2:K, 1]
tol1 <- abs(sum(err1 - RHO1)) / K


err2 <- cor(covmat2)[2:K, 1]
tol2 <- abs(sum(err2 - RHO2)) / K

test_that("Function appropriately switches between fixed and variable seed",{
  expect_error(genNumeric(N, K, rho=-0.5, seed=3))
  expect_error(genNumeric(N, K, rho=0.3, seed=rep(7, 10)))
  expect_error(genNumeric(N, K, rho=0.3, seed=rep(7, N+10)))
  expect_error(genNumeric(N, K, rho=0.3, seed=cbind(rep(7, N+10), rep(7, N+10))))
  expect_error(genNumeric(N, K, rho=0.3, seed=cbind(rep(7, N), rep(7, N))))
})

test_that("Function gets the right answer!", {
  expect_identical(covmat[, 1], S1)
  expect_identical(covmat2[, 1], S2)
  expect_that(tol1, is_less_than(0.05))
  expect_that(tol2, is_less_than(0.1))
})

context("Test user specified seed in patterned structure")

seeds <- genNumeric(1000, 6, rho=0.3)
struc <- list(dist=c("norm", "norm", "chisq", "pois", "pois", "gamma", 
                     "weibull"), 
              rho=c(0.7, 0.3, -0.5, 0.3, -0.8, 0.05, 0.7), 
              names=c("test1", "test2", "noise", "daysattended", 
                      "daysOUT", "bad", "bad2"), 
              seed = cbind(seeds[,1], seeds[,2], seeds[,3], seeds[, 4], seeds[, 5], 
                           seeds[, 6], seeds[,1]))

dat <- genNumeric(1000, seed=TRUE, pattern=struc)

test_that("Data generates the right answer", {
  expect_that(abs(cor(seeds[,1], dat[,1]) - struc$rho[1]), is_less_than(.05))
  expect_that(abs(cor(seeds[,2], dat[,2]) - struc$rho[2]), is_less_than(.05))
  expect_that(abs(cor(seeds[,3], dat[,3]) - struc$rho[3]), is_less_than(.05))
  expect_that(abs(cor(seeds[,4], dat[,4]) - struc$rho[4]), is_less_than(.05))
  expect_that(abs(cor(seeds[,5], dat[,5]) - struc$rho[5]), is_less_than(.05))
  expect_that(abs(cor(seeds[,6], dat[,6]) - struc$rho[6]), is_less_than(.05))
  expect_that(abs(cor(seeds[,1], dat[,7]) - struc$rho[7]), is_less_than(.05))
})



context("Speed")


N <- 98756
K <- 25
P <- 0.43278456
RHO2 <- -0.24
covmat <- genNumeric(N, K, rho=RHO1)

test_that("Function executes in reasonable length of time", {
  expect_that(genNumeric(1000, 25, rho=0.3), takes_less_than(0.5))
  expect_that(genNumeric(500, 40, rho=0.3), takes_less_than(0.5))
  expect_that(genNumeric(87953, 8, rho=0.3), takes_less_than(0.5))
  expect_that(genNumeric(87953, 80, rho=0.3), takes_less_than(1))
})

context("Generate simple correlated factors")

N <- 10000
K <- 30
LEVS <- 12
RHO1 <- 0.2

test <- genFactor(N, K, nlevel=LEVS, rho=RHO1)


test_that("genFactor outputs data in expected size and shape", {
  expect_is(test, "data.frame")
  expect_equal(nrow(test), N)
  expect_equal(ncol(test), K)
  expect_equal(length(table(test[1])), LEVS)
})


context("Are the variables in the data.frame properly related?")


chiTest <- function(i, j, data) {chisq.test(data[,i], data[,j])$p.value}
chiP <- Vectorize(chiTest, vectorize.args=list("i", "j"))
Cresults <- outer(1:K, 1:K, chiP, data=test)

gTest <- function(i, j, data) {gammaGK(data[, i], data[, j])$g}
gP <- Vectorize(gTest, vectorize.args=list("i", "j"))
Gresults <- outer(1:K, 1:K, gP, data=test)

gTestSE <- function(i, j, data) {gammaGK(data[, i], data[, j])$se}
gPse <- Vectorize(gTestSE, vectorize.args=list("i", "j"))
GSEresults <- outer(1:K, 1:K, gPse, data=test)

chiSQ <- rep(NA, K-1)

for(i in 2:ncol(Cresults)){
  chiSQ[i-1] <- Cresults[i-1, i] < 0.05
}

summary(chiSQ)

gMAT <- cbind(rep(NA, K-1),rep(NA, K-1)) 

for(i in 2:ncol(Cresults)){
  gMAT[i-1, 1] <- Gresults[i-1, i] + 2* GSEresults[i-1, i]
  gMAT[i-1, 2] <- Gresults[i-1, i] - 2* GSEresults[i-1, i]
}

tol <- RHO1 * .01
test2 <- abs(round(gMAT[,1], digits=3) - RHO1) > tol


test_that("Bivariate relationships exist and magnitude is correct", {
  expect_equivalent(length(chiSQ), table(chiSQ)["TRUE"])
  expect_that(table(test2)["TRUE"][[1]], is_more_than(round(.9*length(test2))))
})


context("Speed tests")

test_that("Function is not slow", {
  expect_that(genFactor(N, K, nlevel=LEVS, rho=RHO1), takes_less_than(1))
})

context("Generate user specified correlated factors")

N <- 5000
K <- 4
LEVS <- 4
RHO1 <- -0.2

S1 <- sample(letters[1:5], N, replace=TRUE)
S2 <- rnorm(N)
test <- genFactor(N, K, nlevel=LEVS, rho=RHO1, seed=S2)

test2 <- genFactor(N, K, nlevel=LEVS, rho=RHO1, seed=S1)

tol <- 0.1

test_that("Correlations are reasonable", {
  expect_that(abs(gammaGK(test[,1], test[,5])$gamma - RHO1), is_less_than(tol))
  expect_that(abs(gammaGK(test[,1], test[,2])$gamma - RHO1), is_less_than(tol))
  expect_that(abs(gammaGK(test[,1], test[,3])$gamma - RHO1), is_less_than(tol))
  expect_that(abs(gammaGK(test[,1], test[,4])$gamma - RHO1), is_less_than(tol))
  expect_that(abs(gammaGK(test2[,1], test2[,5])$gamma - RHO1), is_less_than(tol))
  expect_that(abs(gammaGK(test2[,1], test2[,2])$gamma - RHO1), is_less_than(tol))
  expect_that(abs(gammaGK(test2[,1], test2[,3])$gamma - RHO1), is_less_than(tol))
  expect_that(abs(gammaGK(test2[,1], test2[,4])$gamma - RHO1), is_less_than(tol))
})

N <- 5000
K <- 4
LEVS <- 4
RHO1 <- -0.9
RHO2 <- 0.99

S1 <- sample(letters[1:5], N, replace=TRUE)
S2 <- rnorm(N)
test <- genFactor(N, K, nlevel=LEVS, rho=RHO1, seed=S2)

test2 <- genFactor(N, K, nlevel=LEVS, rho=RHO2, seed=S1)

test_that("Test extreme values of RHO", {
  expect_that(abs(gammaGK(test[,1], test[,5])$gamma - RHO1), is_less_than(tol))
  expect_that(abs(gammaGK(test[,1], test[,2])$gamma - RHO1), is_less_than(tol))
  expect_that(abs(gammaGK(test[,1], test[,3])$gamma - RHO1), is_less_than(tol))
  expect_that(abs(gammaGK(test[,1], test[,4])$gamma - RHO1), is_less_than(tol))
  expect_that(abs(gammaGK(test2[,1], test2[,5])$gamma - RHO2), is_less_than(tol))
  expect_that(abs(gammaGK(test2[,1], test2[,2])$gamma - RHO2), is_less_than(tol))
  expect_that(abs(gammaGK(test2[,1], test2[,3])$gamma - RHO2), is_less_than(tol))
  expect_that(abs(gammaGK(test2[,1], test2[,4])$gamma - RHO2), is_less_than(tol))
})

context("Can chain together using common starting seed")


struc <- list(dist=c("norm", "norm", "unif", "pois", "pois", "gamma", 
                     "weibull"), 
              rho=c(0.7, 0.3, -0.5, 0.3, -0.8, 0.05, 0.7), 
              names=c("test1", "test2", "noise", "daysattended", 
                      "daysOUT", "bad", "bad2"), 
              seed = cbind(seeds[,1], seeds[,2], seeds[,3], seeds[, 4], seeds[, 5], 
                           seeds[, 6], seeds[,1]))

dat1 <- genFactor(1000, 3, nlevel=3, rho=0.8)
dat2 <- genFactor(1000, 4, nlevel=4, rho=0.3, seed=rnorm(1000))
dat3 <- genFactor(1000, 4, nlevel=6, rho=-0.7, seed=dat2[,4])

identical(dat2[,4], dat3[,1])


test_that("Test that seed is preserved and identical in two data frames",{
  expect_identical(dat2[, 4], dat3[, 1])
  
})

try1 <- genFactor(25000, 4, nlevel=27, rho=0.3, seed=rnorm(25000))
try2 <- genFactor(25000, 4, nlevel=100, rho=0.3, seed=rnorm(25000))
try3 <- genFactor(25000, 4, nlevel=200, rho=0.3, seed=rnorm(25000))
try4 <- genFactor(25000, 4, nlevel=400, rho=0.3, seed=rnorm(25000))



test_that("Factor levels greater than 26 can be generated", {
  expect_equal(length(unique(try1)), 27)
  expect_equal(length(unique(try2)), 100)
  expect_equal(length(unique(try3)), 200)
  expect_equal(length(unique(try4)), 400)
})


context("Generate formulas")

set.seed(382)
seeds <- genNumeric(10000, 6, rho=0.1)

struc <- list(dist=c("norm", "norm", "unif", "pois", "pois", "gamma", 
                     "weibull"), 
              rho=c(0.7, 0.3, -0.5, 0.3, -0.8, 0.05, 0.7), 
              names=c("test1", "test2", "noise", "daysattended", 
                      "daysOUT", "bad", "bad2"), 
              seed = cbind(seeds[,1], seeds[,2], seeds[,3], seeds[, 4], seeds[, 5], 
                           seeds[, 6], seeds[,1]))

dat <- genNumeric(10000, pattern=struc)

dat1 <- genFactor(10000, 3, nlevel=3, rho=0.8)
dat2 <- genFactor(10000, 4, nlevel=4, rho= - 0.1, seed=dat[,6])
dat3 <- genFactor(10000, 4, nlevel=6, rho= -0.2, seed=dat2[,4])
identical(dat2[,4], dat3[,1])

names(dat1) <- sample(LETTERS, length(names(dat1)))
names(dat2) <- sample(letters, length(names(dat2)))
names(dat3) <- sample(paste0(letters,LETTERS), length(names(dat3)))
mdf <- cbind(dat, dat1)
mdf <- cbind(mdf, dat2)
mdf <- cbind(mdf, dat3)
#mdf <- mdf[, c(2:6, 12, 16, 19, 11)]

myF <- list(vars = sample(names(mdf), 7))

genFormula(mdf, myF$vars)

context("Generate binomial dependent variables")

myF$coefs <- rnorm(length(genFormula(mdf, myF$vars)[-1]), mean=0, sd=4)

mdf$out <- genBinomialDV(mdf, form=myF, intercept=-2)
table(mdf$out)

mod1 <- glm(out ~ ., data=mdf, family="binomial")
mod2 <- glm(out ~ bad2 + cC + F + tT + q + hH + r, data=mdf, 
            family="binomial")

context("Generate other dependent variables")

