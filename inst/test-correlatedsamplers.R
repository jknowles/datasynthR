################################################################################

context("Correlations starting with a normally distributed variable")

RHO1 <- 0.7
RHO2 <- -0.7
RHO3 <- 0.01
a <- rnorm(100000)
b <- sapply(a, rnormcor, RHO1)
c <- sapply(a, rnormcor, RHO2)
d <- sapply(a, rnormcor, RHO3)
cor(a,b)

test_that("Correlation is correct", {
  expect_equivalent(cor(a,b), RHO1)
  expect_equivalent(cor(a,c), RHO2)
  expect_equivalent(cor(a,d), RHO3)
  
})

test_that("Direction is correct", {
  expect_equal(sign(cor(a,b)), sign(RHO1))
  expect_equal(sign(cor(a,c)), sign(RHO2))
  expect_equal(sign(cor(a,d)), sign(RHO3))
})



context("Chi-square")

RHO1 <- -.5
a <- rchisq(1000, df=12)
b <- sapply(a, rchisqcor,RHO1)
cor(a,b)

test_that("Correlation of chi-square is correct", {
  expect_equivalent(cor(a,b), RHO1)
})

test_that("Direction is correct", {
  expect_equal(sign(cor(a,b)), sign(RHO1))
  expect_equal(sign(cor(a,c)), sign(RHO2))
  expect_equal(sign(cor(a,d)), sign(RHO3))
})


context("Poisson")


a <- rnorm(1000)
tmp <- sapply(a, rnormcor, rho=0.7)
b <- round(tmp,digits=0)
b <- b + abs(min(b))
b <- b^3

mean(b)
sd(b)
qplot(b,a)
cor(a,b)

test_that("Correlation of poisson is correct", {
  expect_equivalent(cor(a,b), RHO1)
})

test_that("Direction is correct", {
  expect_equal(sign(cor(a,b)), sign(RHO1))
  expect_equal(sign(cor(a,c)), sign(RHO2))
  expect_equal(sign(cor(a,d)), sign(RHO3))
})

test_that("Result is poisson distributed", {
  expect_equivalent(mean(b), sd(b))
  
})

context("Uniform")

RHO1 <- 0.7
RHO2 <- -0.7
RHO3 <- 0.05

a <- runif(1000)
b <- runifcor.cor(a, RHO1)
c <- runifcor.cor(a, RHO2)
d <- runifcor.cor(a, RHO3)


test_that("Correlation is correct", {
  expect_equivalent(cor(a,b), RHO1)
  expect_equivalent(cor(a,c), RHO2)
  expect_equivalent(cor(a,d), RHO3)
  
})

test_that("Direction is correct", {
  expect_equal(sign(cor(a,b)), sign(RHO1))
  expect_equal(sign(cor(a,c)), sign(RHO2))
  expect_equal(sign(cor(a,d)), sign(RHO3))
})

test_that("Result is uniformly distributed", {
  #expect_equivalent(mean(b), sd(b))
  
})

context("Weibull")

library(MASS)
RHO1 <- 0.7
RHO2 <- -0.7
RHO3 <- 0.05
a <- rnorm(1000)
b <- rweibullcor(x, RHO1)
c <- rweibullcor(x, RHO2)
d <- rweibullcor(x, RHO3)


test_that("Correlation is correct", {
  expect_equivalent(cor(a,b), RHO1)
  expect_equivalent(cor(a,c), RHO2)
  expect_equivalent(cor(a,d), RHO3)
  
})


test_that("Direction is correct", {
  expect_equal(sign(cor(a,b)), sign(RHO1))
  expect_equal(sign(cor(a,c)), sign(RHO2))
  expect_equal(sign(cor(a,d)), sign(RHO3))
})

test_that("Result is uniformly distributed", {
  #expect_equivalent(mean(b), sd(b))
  
})
