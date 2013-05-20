##'
##'
##'


# Genereate missingness 


MCARx <- function(x, p){
  z <- rbinom(length(x), 1, prob=p)
  x[z==1] <- NA
  return(x)
}

MCAR.df <- function(df, p){
  if(length(p) == 1){
    df <- apply(df, 2, MCARx, p)
  } else if(length(p) > 1) {
  df <- apply(df, 2, MCARx, sample(p, 1))
  }
  return(df)
}

dimNA <- function(df){
  dims <- dim(df)[1] * dim(df)[2]
  propNA <- apply(df, 2, vecNAsearch)
  countNA <- propNA * dim(df)[1]
  total <- sum(countNA)
  totalP <- total / dims
  return(list("TotalCells" = dims, "MissingbyColumn" = countNA, 
              "TotalMissing" = total, "TotalProportionMissing" = totalP))
  
}

vecNAsearch <- function(x){
  l <- length(x)
  lNA <- length(x[is.na(x)])
  return(lNA / l)
}


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
