##' Assign NAs to a vector randomly
##'
##' Assign NAs to a vector randomly based on a user specified probability
##' 
##' @param x a vector of any class that NA values should be supplied to
##' @param p A value from 0 to 1 representing the proportion of x that should be replaced with NAs
##' @return x with proportion p of values replaced with NA
##' @note Yadda yadda yadda
##' @export
##' @author Jared E. Knowles
MCARx <- function(x, p){
  z <- rbinom(length(x), 1, prob=p)
  x[z==1] <- NA
  return(x)
}

##' Assign NAs to columns in a dataframe completely at random
##'
##' Create a dataframe with data missing completely at random
##' 
##' @param df a dataframe that missing data should be supplied to
##' @param p A value from 0 to 1 representing the proportion of df that should be replaced with NAs
##' @return df with proportion p of values replaced with NA
##' @note Yadda yadda yadda
##' @export
##' @author Jared E. Knowles
MCAR.df <- function(df, p){
  if(length(p) == 1){
    df <- apply(df, 2, MCARx, p)
  } else if(length(p) > 1) {
  df <- apply(df, 2, MCARx, sample(p, 1))
  }
  df <- as.data.frame(df)
  return(df)
}

##' Calculate missingness in a data.frame
##'
##' Figure out the level and type of missingness in a dataframe
##' 
##' @param df a dataframe with missing values
##' @return a list indicating the total size of df, a vector of counts of NA terms by columns
##' a count of al missing values in the df, and the total proportion of data missing
##' @note Yadda yadda yadda
##' @export
##' @author Jared E. Knowles
dimNA <- function(df){
  dims <- dim(df)[1] * dim(df)[2]
  propNA <- apply(df, 2, vecNAsearch)
  countNA <- propNA * dim(df)[1]
  total <- sum(countNA)
  totalP <- total / dims
  return(list("TotalCells" = dims, "MissingbyColumn" = countNA, 
              "TotalMissing" = total, "TotalProportionMissing" = totalP))
  
}

##' Count the proportion of NA values in a vector
##'
##' Return the proportion of vector x that are NA
##' 
##' @param x a vector of any class that contains NA values
##' @return the proportion of x that are NA values
##' @note Yadda yadda yadda
##' @export
##' @author Jared E. Knowles
vecNAsearch <- function(x){
  l <- length(x)
  lNA <- length(x[is.na(x)])
  return(lNA / l)
}

##' Are missing values in two vectors correlated?
##'
##' Return the GK test
##' 
##' @param x a vector of any class that contains NA values
##' @param y a vector of any class that contains NA values
##' @return the proportion of x that are NA values
##' @note Yadda yadda yadda
##' @details This function forces x and y to factors of missing and non-missing 
##' elements and then performs a GK test on those factors to determine if they are 
##' correlated, to what degree, and if this is statistically significant. 
##' @export
##' @author Jared E. Knowles
MCARcheck <- function(x, y){
  x1 <- rep("Missing", length(x))
  x1[is.na(x)==FALSE] <- "Not Missing"
  x1 <- factor(x1)
  y1 <- rep("Missing", length(y))
  y1[is.na(y)==FALSE] <- "Not Missing"
  y1 <- factor(y1)
  return(gammaGK(x1, y1, print=FALSE))
}


##' Test correlation of missingness for all pairwise combinations in a dataframe
##'
##' Return the GK test for all pairs of columns in a dataframe on missing values
##' 
##' @param df a dataframe
##' @return a list with three matrices containing the pairwise Gamma statistic for 
##' all columns in the dataframe, the standard error of the Gamma statistic, and 
##' the significance test of the Gamma statistic
##' @note Yadda yadda yadda
##' @details This function forces each pair of columns to factors of missing and 
##' non-missing elements and then performs a GK test on those factors to 
##' determine if they are correlated, to what degree, and if this is 
##' statistically significant. 
##' @export
##' @author Jared E. Knowles
MCARcheck.df <- function(df){
  # Add a check to subset on data that is missing
  # significance test
  missingTest1 <- function(i, j, data) {MCARcheck(data[,i], data[,j])$sig}
  mTP1 <- Vectorize(missingTest1, vectorize.args=list("i", "j"))
  K <- ncol(df)
  results1 <- outer(1:K, 1:K, mTP1, data=df)
  # gamma stat
  missingTest2 <- function(i, j, data) {MCARcheck(data[,i], data[,j])$gamma}
  mTP2 <- Vectorize(missingTest2, vectorize.args=list("i", "j"))
  results2 <- outer(1:K, 1:K, mTP2, data=df)
  # standard errors
  missingTest3 <- function(i, j, data) {MCARcheck(data[,i], data[,j])$se}
  mTP3 <- Vectorize(missingTest3, vectorize.args=list("i", "j"))
  results3 <- outer(1:K, 1:K, mTP3, data=df)
  return(list(gammas = results2, se = results3, sig.test = results1, 
              names = names(df)))
}

##' Summary for MCAR checks
##' 
##' @param results results from MCARcheck.df
##' @param p threshold for statistical significance on the gamma
##' @param print Option to print out results
##' @return summary of dependencies among missing values in dataset
##' @export
##' @author Jared E. Knowles
summary.MCARcheck <- function(results, p, print){
  library(reshape)
  m1 <- melt(results$gammas)
  names(m1)[3] <- "gamma"
  m2 <- melt(results$se)
  names(m2)[3] <- "se"
  m3 <- melt(results$sig)
  names(m3)[3] <- "sig"
  df <- merge(m1, m2)
  df <- merge(df, m3)
  df$X1 <- factor(df$X1, labels=results$names)
  df$X2 <- factor(df$X2, labels=results$names)
  names(df)[1] <- "Var1"
  names(df)[2] <- "Var2"
  df <- subset(df, Var1 != Var2)
  df <- cull(df)
  # Build summary items
  p <- .05
  df.tmp <- df[df$sig < p, ]
  
  if(print==TRUE){
    cat("Goodman-Kruskal Gamma Statistics for Missingness:\n")
    cat(paste("Total Pairs of Variables:", nrow(df), "\n"))
    cat(paste("Variable Pairs With Correlated Missingness:",nrow(df.tmp),"\n"))
    cat(paste("Variable Pairs Without Correlated Missingness:",
              nrow(df) - nrow(df.tmp),"\n\n"))
    cat(paste("Values of statistically significant gamma: ",
              unique(df.tmp$gamma),"\n"))
  }
  return(list(totalpairs = nrow(df), correlatedmissingpairs = nrow(df.tmp), 
              noncorrelatedmissing = nrow(df) - nrow(df.tmp), 
              gammas = unique(df.tmp$gamma), 
              percentcorrelated = nrow(df.tmp)/ nrow(df), 
              data = df))
}

##' MCAR Plot
##' 
##' @param results results from MCARcheck.df
##' @return ggplot2 object 
##' @export
##' @author Jared E. Knowles
MCARplot <- function(results){
  z.m <- melt(results$gammas)
  ggplot(z.m, aes(X1, X2, fill = value)) + geom_tile() + 
    scale_fill_gradient2(low = "blue",  high = "yellow") + 
    coord_cartesian(xlim=c(1, max(z.m$X1)), ylim=c(1, max(z.m$X2))) +
    labs(x="Var1", y="Var2", 
         title="Gamma Coefficients for Missingness Among \n Variables in Dataset")
}

##' Assign NAs to columns in a dataframe at random
##'
##' Create a dataframe with data missing at random by generating probability models
##' of missing data from observable characteristics and then eliminating data based 
##' on those results.
##' 
##' @param df a dataframe that missing data should be supplied to
##' @param probs A single value from 0 to 1 representing the proportion of each column 
##' that should be replaced with NAs, or a list of such values length of \code{vars}
##' @param vars A list of variable names to be used in generating the probability models 
##' for the missing data
##' @return df with proportion p of values replaced with NA
##' @note Yadda yadda yadda
##' @export
##' @author Jared E. Knowles
MAR.df <- function(df, vars, probs){
  "%w/o%" <- function(x, y) x[!x %in% y] #--  x without y
  if(length(probs) == 1){
    probs <- rep(probs, length(vars))
  } else if(length(probs) > 1){
    if(length(probs) != length(vars)) stop("Lengths don't match.")
  }
  for(i in 1:length(vars)){
    myF <- list(vars=sample(names(df) %w/o% vars, 4))
    myF$coefs <- rnorm(length(genFormula(df, myF$vars)[-1]), mean=0, sd=1)
    eval(parse(text=paste0("df$out", i, "<- genBinomialDV(df, form=myF, intercept=0, type='response')")))
  }
  for(j in 1:ncol(df)){
    s <- sample(1:length(vars), 1)
    eval(parse(text=paste0("q <- quantile(df$out",s ," , probs[",s,"], na.rm=T)")))
    eval(parse(text=paste0("df[,", j, "][df$out", s," < ", q, "] <- NA")))
  }
  return(df)
}

##' Missing Not at Random


