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


contxt("Chi-square")

RHO1 <- -.5
a <- rchisq(1000, df=12)
b <- sapply(a, rchisqcor,RHO1)
cor(a,b)

test_that("Correlation of chi-square is correct", {
  expect_equivalent(cor(a,b), RHO1)
})




a <- rnorm(1000)
b <- sapply(a, rnormcor, 0.7)
cor(a,b)



a <- rpois(1000, lambda=4)
b <- sapply(a, rchisqcor, rho=-.5)
cor(a,b)

qplot(a,b) + geom_smooth()