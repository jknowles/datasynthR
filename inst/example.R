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

seeds <- genNumeric(10000, 6, rho=0.3)

struc <- list(dist=c("norm", "norm", "unif", "pois", "pois", "gamma", 
                     "weibull"), 
              rho=c(0.7, 0.3, -0.5, 0.3, -0.8, 0.05, 0.7), 
              names=c("test1", "test2", "noise", "daysattended", 
                      "daysOUT", "bad", "bad2"), 
              seed = cbind(seeds[,1], seeds[,2], seeds[,3], seeds[, 4], seeds[, 5], 
                       seeds[, 6], seeds[,1]))

dat <- genNumeric(10000, pattern=struc)

# cor(seeds[,1], dat[,1])
# cor(seeds[,2], dat[,2])
# cor(seeds[,3], dat[,3])
# cor(seeds[,4], dat[,4])
# cor(seeds[,5], dat[,5])
# cor(seeds[,6], dat[,6])
# cor(seeds[,1], dat[,7])

dat1 <- genFactor(10000, 3, nlevel=3, rho=0.8)
dat2 <- genFactor(10000, 4, nlevel=4, rho=0.3, seed=dat[,6])
dat3 <- genFactor(10000, 4, nlevel=6, rho=-0.7, seed=dat2[,4])

identical(dat2[,4], dat3[,1])

names(dat1) <- sample(LETTERS, length(names(dat1)))
names(dat2) <- sample(letters, length(names(dat2)))
names(dat3) <- sample(paste0(letters,LETTERS), length(names(dat3)))

mdf <- cbind(dat, dat1)
mdf <- cbind(mdf, dat2)
mdf <- cbind(mdf, dat3)

mdf <- mdf[, c(2:6, 12, 16, 19, 11)]

myF <- list(vars = sample(names(mdf), 5))

genFormula(mdf, myF$vars)

myF$coefs <- rnorm(length(genFormula(mdf, myF$vars)[-1]))

mdf$out <- genBinomialDV(mdf, form=myF, intercept=-50)
table(mdf$out)

mod1 <- glm(out ~ ., data=mdf, family="binomial")
mod2 <- glm(out ~ noise + c + test1 + out + aA, data=mdf, 
            family="binomial")


exp <- paste0(myF$coefs, "*","mod$", genFormula(mdf, myF$vars)[-1], collapse=" + ")
myF$exp <- parse(text=exp)
mm <- eval(parse(text=paste0("terms(~", paste0(myF$vars, collapse=' + '),")")))
mod <- as.data.frame(model.matrix(mm, mdf))
mod$z <- 150 + eval(myF$exp) + runif(dim(mod)[1])
#z <- ifelse(scale==TRUE, scale(z, center=TRUE)[,1], z)

mod$pr = 1/(1+exp(-mod$z))         # pass through an inv-logit function
summary(mod$pr)
mod$y = rbinom(dim(mod)[1],1,mod$pr) 

summary(glm(y ~ . -z -pr , data=mod, family="binomial"))

table(mod$y)
summary(mod)

mod$z <- scale(mod$z)
mod$b <- pnorm(mod$z)
mod$b <- rbinom(dim(mod)[1], 1, pr=mod$b)
  

mod$pr <- rnormcorV(mod$z, 0.1)
mod$pr <- rbinomcor(scale(mod$z), 0.99)

mod$pr <- 1/(1+exp(-mod$z))
y = rbinom(dim(mod)[1], 1, pr)
return(y)

###########################################################
exp <- paste0(myF$coefs, "*","mod$", genFormula(mdf, myF$vars)[-1], collapse=" + ")
myF$exp <- parse(text=exp)


genFormula(mdf, myF$vars)
z <- eval(parse(text=paste0("terms(~", paste0(myF$vars, collapse=' + '),")")))
test <- model.matrix(z, mdf)




exp <- paste0(myF$coefs, "*","mod$", genFormula(mdf, myF$vars)[-1], collapse=" + ")
myF$exp <- parse(text=exp)
mm <- eval(parse(text=paste0("terms(~", paste0(myF$vars, collapse=' + '),")")))
mod <- as.data.frame(model.matrix(mm, mdf))

z <- 2 + eval(myF$exp) + runif(dim(mod)[1])
pr <- 1/(1+exp(scale(-z, center=TRUE)))
y = rbinom(dim(mod)[1], 1, pr)
return(y)

coefs <- runif(length(genFormula(mdf, myF$vars)), min=-2, max=2)
myF$coefs <- coefs


vars <- names(as.data.frame(test))
coefs <- runif(length(vars), min=-2, max=2)

out <- coefs * test
out <- apply(out, 1, sum)

exp <- paste0(myF$coefs, "*","df$", myF$vars, collapse="+")


z <- eval(parse(text=paste0("terms(.~", paste0(myF$vars, collapse=' + '),")")))



terms( . ~mdf[,9] + mdf[,4])

mdf$out <- genBinomialDV(mdf, form=myF, intercept=3)


