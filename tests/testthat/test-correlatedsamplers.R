################################################################################
################################################################################
# Test the direction and strength of correlations for various combinations of 
#  distributions
# Start with normal seed and correlations from then
# Other distributions:
#    poisson
#    weibull
#    uniform
#    gamma
#    beta
#    binomial
#    chisq

context("Correlations starting with a normally distributed variable")

set.seed(347829)

RHO1 <- 0.7
RHO2 <- -0.7
RHO3 <- 0.01
a <- rnorm(100000)
b <- sapply(a, rnormcor, RHO1)
c <- sapply(a, rnormcor, RHO2)
d <- sapply(a, rnormcor, RHO3)

tol <- .01

test_that("Correlation is correct", {
  expect_that(abs(cor(a,b) - RHO1), is_less_than(tol))
  expect_that(abs(cor(a,c) - RHO2), is_less_than(tol))
  expect_that(abs(cor(a,d) - RHO3), is_less_than(tol))
  
})

test_that("Direction is correct", {
  expect_equal(sign(cor(a,b)), sign(RHO1))
  expect_equal(sign(cor(a,c)), sign(RHO2))
  expect_equal(sign(cor(a,d)), sign(RHO3))
})

context("Normal Vector")

RHO1 <- 0.7
RHO2 <- -0.7
RHO3 <- 0.01
a <- rnorm(100000)
b <- rnormcorV(a, RHO1)
c <- rnormcorV(a, RHO2)
d <- rnormcorV(a, RHO3)

tol <- .01

test_that("Correlation is correct", {
  expect_that(abs(cor(a,b) - RHO1), is_less_than(tol))
  expect_that(abs(cor(a,c) - RHO2), is_less_than(tol))
  expect_that(abs(cor(a,d) - RHO3), is_less_than(tol))
  
})

test_that("Direction is correct", {
  expect_equal(sign(cor(a,b)), sign(RHO1))
  expect_equal(sign(cor(a,c)), sign(RHO2))
  expect_equal(sign(cor(a,d)), sign(RHO3))
})


context("Chi-square")

a <- rnorm(5000)
RHO1 <- 0.7
RHO2 <- -0.7
RHO3 <- 0.01
b <- rchisqcor(a, RHO1)
c <- rchisqcor(a, RHO2)
d <- rchisqcor(a, RHO3)
tol <- .08

test_that("Correlation of chi-square is correct", {
  expect_that(abs(cor(a,b) - RHO1), is_less_than(tol))
  expect_that(abs(cor(a,c) - RHO2), is_less_than(tol))
  expect_that(abs(cor(a,d) - RHO3), is_less_than(tol))
})

test_that("Direction is correct", {
  expect_equal(sign(cor(a,b)), sign(RHO1))
  expect_equal(sign(cor(a,c)), sign(RHO2))
  expect_equal(sign(cor(a,d)), sign(RHO3))
})


context("Poisson")


a <- rnorm(5000)
RHO1 <- 0.7
RHO2 <- -0.7
RHO3 <- 0.1

b <- rpoiscor(a, RHO1)
c <- rpoiscor(a, RHO2)
d <- rpoiscor(a, RHO3)
tol <- 0.05

test_that("Correlation of poisson is correct", {
  expect_that(abs(cor(a,b) - RHO1), is_less_than(tol))
  expect_that(abs(cor(a,c) - RHO2), is_less_than(tol))
  expect_that(abs(cor(a,d) - RHO3), is_less_than(tol))
})

test_that("Direction is correct", {
  expect_equal(sign(cor(a,b)), sign(RHO1))
  expect_equal(sign(cor(a,c)), sign(RHO2))
  expect_equal(sign(cor(a,d)), sign(RHO3))
})

test_that("Result is poisson distributed", {
  expect_equal(mean(b), sd(b), tol = 0.4)
  expect_equal(mean(c), sd(c), tol = 0.4)
  expect_equal(mean(d), sd(d), tol = 0.4)
})

context("Uniform")

RHO1 <- 0.7
RHO2 <- -0.7
RHO3 <- 0.05

a <- runif(1000)
b <- runifcor.cor(a, RHO1)
c <- runifcor.cor(a, RHO2)
d <- runifcor.cor(a, RHO3)

# 
# test_that("Correlation is correct", {
#   expect_that(abs(cor(a,b) - RHO1), is_less_than(tol))
#   expect_that(abs(cor(a,c) - RHO2), is_less_than(tol))
#   expect_that(abs(cor(a,d) - RHO3), is_less_than(tol))
#   
# })
# 
# test_that("Direction is correct", {
#   expect_equal(sign(cor(a,b)), sign(RHO1))
#   expect_equal(sign(cor(a,c)), sign(RHO2))
#   expect_equal(sign(cor(a,d)), sign(RHO3))
# })
# 
# test_that("Result is uniformly distributed", {
#   #expect_equivalent(mean(b), sd(b))
#   
# })

context("Weibull")

library(MASS)
RHO1 <- 0.7
RHO2 <- -0.7
RHO3 <- 0.05
a <- rnorm(1000)
b <- rweibullcor(a, RHO1)
c <- rweibullcor(a, RHO2)
d <- rweibullcor(a, RHO3)


test_that("Correlation is correct", {
  expect_that(abs(cor(a,b) - RHO1), is_less_than(tol))
  expect_that(abs(cor(a,c) - RHO2), is_less_than(tol))
  expect_that(abs(cor(a,d) - RHO3), is_less_than(tol))
  
})


test_that("Direction is correct", {
  expect_equal(sign(cor(a,b)), sign(RHO1))
  expect_equal(sign(cor(a,c)), sign(RHO2))
  expect_equal(sign(cor(a,d)), sign(RHO3))
})

test_that("Result is uniformly distributed", {
  #expect_equivalent(mean(b), sd(b))
  
})

context("Gamma")


RHO1 <- 0.3
RHO2 <- -0.3
RHO3 <- 0.05
a <- rnorm(1000)
b <- rgammacor(a, RHO1)
c <- rgammacor(a, RHO2)
d <- rgammacor(a, RHO3)


test_that("Correlation is correct", {
  expect_that(abs(cor(a,b) - RHO1), is_less_than(tol))
  expect_that(abs(cor(a,c) - RHO2), is_less_than(tol))
  expect_that(abs(cor(a,d) - RHO3), is_less_than(tol))
  
})


test_that("Direction is correct", {
  expect_equal(sign(cor(a,b)), sign(RHO1))
  expect_equal(sign(cor(a,c)), sign(RHO2))
  expect_equal(sign(cor(a,d)), sign(RHO3))
})

test_that("Result is gamma distributed", {
  #expect_equivalent(mean(b), sd(b))
  
})

context("Negative Binomial")


RHO1 <- 0.7
RHO2 <- -0.7
RHO3 <- 0.05
a <- rnorm(1000)
b <- rnegbinomcor(a, RHO1)
c <- rnegbinomcor(a, RHO2)
d <- rnegbinomcor(a, RHO3)

# cor1 <- glm(b ~ a, family="binomial")
# cor2 <- glm(c ~ 0 + a, family="binomial")
# cor3 <- glm(d ~ 0 + a, family="binomial")

test_that("Correlation is correct", {
  expect_that(abs(cor(a,b) - RHO1), is_less_than(tol))
  expect_that(abs(cor(a,c) - RHO2), is_less_than(tol))
  expect_that(abs(cor(a,d) - RHO3), is_less_than(tol))
  
})


test_that("Direction is correct", {
  expect_equal(sign(cor(a,b)), sign(RHO1))
  expect_equal(sign(cor(a,c)), sign(RHO2))
  expect_equal(sign(cor(a,d)), sign(RHO3))
})

tol <- 25

test_that("Result is negative binomial distributed", {
  expect_that(length(table(b)), is_less_than(tol))
  expect_that(length(table(c)), is_less_than(tol))
  expect_that(length(table(d)), is_less_than(tol))
})

context("Binomial")


RHO1 <- 0.7
RHO2 <- -0.7
RHO3 <- 0.05
a <- rnorm(1000)
b <- rbinomcor(a, RHO1)
c <- rbinomcor(a, RHO2)
d <- rbinomcor(a, RHO3)

tol <- 0.25

test_that("Correlation is correct", {
  expect_that(abs(cor(a,b) - RHO1), is_less_than(tol))
  expect_that(abs(cor(a,c) - RHO2), is_less_than(tol))
  expect_that(abs(cor(a,d) - RHO3), is_less_than(tol))
  
})


test_that("Direction is correct", {
  expect_equal(sign(cor(a,b)), sign(RHO1))
  expect_equal(sign(cor(a,c)), sign(RHO2))
  expect_equal(sign(cor(a,d)), sign(RHO3))
})

test_that("Result is binomial distributed", {
  expect_equivalent(length(table(b)), 2)
  expect_equivalent(length(table(c)), 2)
  expect_equivalent(length(table(d)), 2)
  
})



#context("Correlations starting with a normally distributed variable")
