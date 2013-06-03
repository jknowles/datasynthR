################################################################################
# Test missingness
################################################################################

context("vector NA search")

test1 <- c(rnorm(100), rep(NA, 28))
test2 <- c(rnorm(1000), rep(NA, 7))
test3 <- c(rep(NA, 400), rnorm(20))

test_that("Correct randomness is detected", {
  expect_equivalent(vecNAsearch(test1), 28/128)
  expect_equivalent(vecNAsearch(test2), 7/1007)
  expect_equivalent(vecNAsearch(test3), 400/420)
})

test_that("Function is fast!", {
  expect_that(vecNAsearch(test1), takes_less_than(.002))
  expect_that(vecNAsearch(test2), takes_less_than(.002))
  expect_that(vecNAsearch(test3), takes_less_than(.002))
  
})



context("Vector MCAR")

P <- 0.1
N <- 1000
x <- rnorm(N)
y <- MCARx(x, P)
a <- length(y[is.na(y)])
t <- .05
e <- N * P
a2 <- abs(e - a) / N

test_that("Correct random chance of missingness is created", {
  expect_that(a2, is_less_than(t))
})

P <- 0.456
N <- 100000
x <- rnorm(N)
y <- MCARx(x, P)
a <- length(y[is.na(y)])
t <- .05
e <- N * P
a2 <- abs(e - a) / N

test_that("Function is fast and accurate", {
  expect_that(MCARx(x,P), takes_less_than(.02))
  expect_that(a2, is_less_than(t))
})


context("Matrix DF")


N <- 8000
K <- 10
P <- 0.1
RHO1 <- 0.45
covmat <- genNumeric(N, K, rho=RHO1)

covmat2 <- MCAR.df(covmat, P)

a <- apply(covmat2, 2, vecNAsearch) * N
t <- .05
e <- K * N * P
a2 <- abs(sum(a) - e) / N

test_that("Function is fast and accurate", {
  expect_that(MCAR.df(covmat,P), takes_less_than(.05))
  expect_that(a2, is_less_than(t))
})

N <- 98756
K <- 25
P <- 0.43278456
RHO2 <- -0.24
covmat <- genNumeric(N, K, rho=RHO2)

covmat2 <- MCAR.df(covmat, P)

a <- apply(covmat2, 2, vecNAsearch) * N
t <- .05
e <- K * N * P
a2 <- abs(sum(a) - e) / N

test_that("Function is fast and accurate", {
  expect_that(MCAR.df(covmat,P), takes_less_than(1))
  expect_that(a2, is_less_than(t))
})


test_that("Function produces data MCAR", {
  
})

context("dimNA function")


N <- 9874
K <- 25
P <- 0.43278456
RHO2 <- -0.24
covmat <- genNumeric(N, K, rho=RHO2)

covmat2 <- MCAR.df(covmat, P)

MR <- dimNA(covmat2)


#a <- apply(covmat2, 2, vecNAsearch) * N
t <- .005
e <- K * N * P
#a2 <- abs(sum(a) - e) / N

test_that("Function returns the correct dimensions answers", {
  expect_equivalent(MR$TotalCells, N*K)
  expect_that(abs(MR$TotalMissing - e) / (N*K), is_less_than(t))
  expect_that(abs(MR$TotalProportionMissing - P), is_less_than(t))
})

test_that("Function is fast!", {
  expect_that(dimNA(covmat2), takes_less_than(.1))
})




context("MAR data")

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

myF <- list(vars = sample(names(mdf), 4))

genFormula(mdf, myF$vars)

#myF$coefs <- rnorm(length(genFormula(mdf, myF$vars)[-1]), mean=0, sd=4)

#genFormula(mdf, myF$vars)

# myF$coefs <- rnorm(length(genFormula(mdf, myF$vars)[-1]), mean=0, sd=4)
# 
# mdf$out.1 <- genBinomialDV(mdf, form=myF, intercept=-2, type="response")

"%w/o%" <- function(x, y) x[!x %in% y] #--  x without y


misslist <- sample(names(mdf)[c(-1, -22)], 5)
probs <- 0.12

MAR.df <- function(df, vars, probs){
  if(length(probs) == 1){
    probs <- rep(probs, length(vars))
  } else if(length(probs) > 1){
    if(length(probs) != length(vars)) stop("Lengths don't match.")
  }
  for(i in 1:length(vars)){
    myF <- list(vars=sample(names(df) %w/o% vars, 4))
    myF$coefs <- rnorm(length(genFormula(df, myF$vars)[-1]), mean=0, sd=1)
    eval(parse(text=paste0("df$out", i, "<- genBinomialDV(df, form=myF, intercept=0, type='response')")))
    eval(parse(text=paste0("df$",misslist[i],"[df$out",i," < df$out1[cutoff(df$out1,",probs[i],")]] <- NA")))
  }
  return(df)
}

mdf2 <- MAR.df(mdf, vars=misslist, probs=probs)



length(misslist)

names(mdf)

library(eeptools)

mdf$out1[cutoff(mdf$out1, .1)]



# 
# 
# N <- 5000
# K <- 25
# P <- 0.43278456
# RHO2 <- -0.24
# covmat <- genNumeric(N, K, rho=RHO2)
# 
# 
# covmat2 <- MAR(covmat, 1)

