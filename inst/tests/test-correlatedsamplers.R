################################################################################

context("Correlations starting with a normally distributed variable")

RHO1 <- 0.7
RHO2 <- -0.7
RHO3 <- 0.01
a <- rnorm(100000)
b <- sapply(a, rnormcor, RHO1)
c <- sapply(a, rnormcor, RHO2)
d <- sapply(a, rnormcor, RHO3)


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

a <- rchisq(1000, df=12)
RHO1 <- 0.7
RHO2 <- -0.7
RHO3 <- 0.01
b <- sapply(a, rchisqcor, RHO1)
c <- sapply(a, rchisqcor, RHO2)
d <- sapply(a, rchisqcor, RHO3)

test_that("Correlation of chi-square is correct", {
  expect_equivalent(cor(a,b), RHO1)
  expect_equivalent(cor(a,c), RHO2)
  expect_equivalent(cor(a,d), RHO3)
})

test_that("Direction is correct", {
  expect_equal(sign(cor(a,b)), sign(RHO1))
  expect_equal(sign(cor(a,c)), sign(RHO2))
  expect_equal(sign(cor(a,d)), sign(RHO3))
})


context("Poisson")

a <- rnorm(1000)
RHO1 <- 0.7
RHO2 <- -0.7
RHO3 <- 0.01
b <- sapply(a, rpoiscor, RHO1)
c <- sapply(a, rpoiscor, RHO2)
d <- sapply(a, rpoiscor, RHO3)


test_that("Correlation of poisson is correct", {
  expect_equivalent(cor(a,b), RHO1)
  expect_equivalent(cor(a,c), RHO2)
  expect_equivalent(cor(a,d), RHO3)
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
b <- rweibullcor(a, RHO1)
c <- rweibullcor(a, RHO2)
d <- rweibullcor(a, RHO3)


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

context("Gamma")


RHO1 <- 0.7
RHO2 <- -0.7
RHO3 <- 0.05
a <- rnorm(1000)
b <- rgammacor(a, RHO1)
c <- rgammacor(a, RHO2)
d <- rgammacor(a, RHO3)


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

test_that("Result is gamma distributed", {
  #expect_equivalent(mean(b), sd(b))
  
})

context("Binomial")


RHO1 <- 0.7
RHO2 <- -0.7
RHO3 <- 0.05
a <- rnorm(1000)
b <- rbinomcor(a, RHO1)
c <- rbinomcor(a, RHO2)
d <- rbinomcor(a, RHO3)

# cor1 <- glm(b ~ a, family="binomial")
# cor2 <- glm(c ~ 0 + a, family="binomial")
# cor3 <- glm(d ~ 0 + a, family="binomial")

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

test_that("Result is binomial distributed", {
  expect_equivalent(length(table(b)), 2)
  expect_equivalent(length(table(c)), 2)
  expect_equivalent(length(table(d)), 2)

})
