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
  df <- as.data.frame(df)
  }
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
  return(list(gammas = results2, se = results3, sig.test = results1))
}

# MAR




# # Consider building a probability model into the function to calculate missingness
# MAR <- function(df, scale){
#   s <- scale
#   n <- ncol(df) - 2
#   for(i in 1:n){
#     cor <- cor(df[, i], df[, i+1])
#     cor <- abs(cor)
#     df[, i:i+1] <- apply(df[, i:i+1], 2, MCARx, s*cor)
#   }  
# }
