datasynthR
==========

Functions to procedurally generate synthetic data in R for testing and collaboration.

`datasynthR` allows the user to generate data of known distributional properties with 
known correlation structures. This is useful for testing statistical model data, building functions to operate on very large datasets, or training others in using R!

`datasynthR` is built around a simple architecture of correlated samplers, dataset generators, and missingness functions. This work was inspired by a 
[list-serv](https://stat.ethz.ch/pipermail/r-help/2007-October/142556.html) post by Joris Dewolf detailing how to generate random correlated 
variables from a normal distribution. 

- *correlated samplers* take a vector V of length n and generate another vector of length n correlated with V at a user specified level of *rho* following the 
distributional form specified in the sampler
- *generators* combine *correlated samplers* together to generate datasets of 
known structure, with variables of known distributions, with known correlations 
among them
- *missingness functions* apply specifc types of missingness to datasets in order 
to simulate real world data

## Correlated samplers

These are the workhorse functions of `datasynthR`. These functions allow for the 
construction of bivariate correlated relationships. Creating correlated vectors 
is as simple as calling the function on the vector you wish to correlate the 
data with and specifying the desired correlation: 

```
a <- rnorm(100)
b <- rnormcorV(a, rho=0.2)

cor(a, b)

[1] 0.19146809
```

The following distributions are impelmented in `datasynthR`: 

- normal (`rnormcorV`)
- chi-squared (`rchisqcor`)
- Poisson (`rpoiscor`)
- Gamma (`rgammacor`)
- uniform (`runifcor`)
- Weibull (`rweibullcor`)
- Negative Binomial (`rnegbinomcor`)
- Binomial (`rbinomcor`)

These functions are able to provide close approximations of `rho` specified when 
the intial vector fed to them is approximately normal and of sufficient length. 
Currently mileage may vary when coming from other distributions. Additionally, the 
distributional properties of these distributions is currently not allowed to 
be manipulated by the user. This may result in non-canonical shapes of these 
distributions currently.

The following distributions will be added in the future: 

- t 
- beta
- geometric


## Generators

Generators build on the correlated samplers by allowing the user to describe 
a data frame structure as a named list and then generate data that meets that 
structure. Currently the two main generators are to generate numeric (`genNumeric`)
or factor (`genFactor`) data. These generators are designed to take a variety 
of high level descriptions of the data structure from the user and generate 
a data frame that approximates that data. 

```
# Generate 1000 cases 3 variables with 0.3 bivariate correlations
dat1 <- genNumeric(1000, 3, rho=0.3)
cor(dat1)

# Note that var1 and var 3 are not correlated at 0.3

# Generate 1000 cases with user specified correlations

struc <- list(dist=c("norm", "pois", "gamma"), rho=c(0.2, -.5, .5),
  names=c("super", "cool", "data"))

dat2 <- genNumeric(1000, pattern=struc)

cor(dat2[, 1], dat2[, 2])
cor(dat2[, 2], dat2[, 3])
cor(dat2[, 3], dat2[, 4])

# Note correlations are not perfect, but the general structure is preserved

# A bit more complicated case passing different seeds

struc3 <- list(dist = c("norm", "chisq", "pois", "norm", 
  "weibull", "gamma"), 
  rho = c(-.05, -.4, 0.3, 0.9, .03, -.6),
  names = c("score", "accept", "score2", "days", "days2", "luck"),
  seed = c(runif(100), rpois(100, 7), rpois(100, 3), 
  rgamma(100, shape=2), runif(100)))

covmat3 <- genNumeric(100, pattern = struc3)

cor(covmat3)
```

Notice that we can also pass names of variables to `genNumeric`. Our next 
step is to generate correlated factors. For now this function is restricted to 
constructing **unordered** factors that have randomly generated level names made 
up of sampling from combinations of the R `letters` and `LETTERS` vectors. 

```
N <- 1000
K <- 4
LEVS <- 5
RHO1 <- -0.2

S1 <- sample(letters[1:5], N, replace=TRUE)
S2 <- rnorm(N)

test  <- genFactor(N, K, nlevel=LEVS, rho=0.3)
test2 <- genFactor(N, K, nlevel=LEVS, rho=RHO1, seed=S1)
test3 <- genFactor(N, K, nlevel=LEVS, rho=RHO1, seed=S2)

```

Assessing `rho` for correlations among unordered factors can be difficult. 
`datasynthR` includes a function for estimating the Goodman and Kruskal gamma statistic adapted from [Simon Jackman](http://jackman.stanford.edu/classes/151B/06/class0517.r). Using the `gammaGK` function we can calculate whether or not 
our factors are related to one another in the desired ways:

```

gammaGK(test[,1], test[,2])
gammaGK(test2[,1], test2[,4])
gammaGK(test3[,1], test3[,5])

```

There is currently no way to assign an ordered structure to these factors. Also the 
level of precision we are able to achieve with our `rho` depends on the number 
of observations, the number of levels in the factor, and underlying distribution 
of the seed variable. 

## Missingness

Often we want to test our method, procedure, or analysis for robustness against 
missing data. `datasynthR` provides a number of convenient methods to simulate 
missingness in data -- whether it was generated in `datasynthR` or not. 

Missing data can come in three forms:

1. Missing completely at random (MCAR)
2. Missing at random (MAR)
3. Missing not at random (MNAR)

Currently MCAR is well implemented in `datasynthR` and MAR is experimentally 
implemented though not robust to most underlying data yet. 

```

covmatMISSING <- MCAR.df(covmat, 0.1)
MCARtest <- MCARcheck.df(covmatMISSING)
summary.MCARcheck(MCARtest)

```

This allows us to ensure that missingness was generated completely at random. These cehcks forces each pair of columns into factors with levels for missing and 
non-missing elements and then performs a Goodman and Kruskal gamma test on those
factors to determine if they are correlated, to what degree, and if this is 
statistically significant. 


