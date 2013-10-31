##' Estimate the Goodman and Kruskal gamma statistic
##'
##' Estimate the correlation between two unordered factor variables using the Goodman and Kruskal gamma statistic 
##' 
##' @param x an unordered factor variable
##' @param y an unordered factor variable
##' @param print a logical vector indicating whether results should be printed to the console
##' @return A named list with gamma, standard error of gamma, p-value of gamma, and statistical significance
##' @note Yadda yadda yadda
##' @export
##' @author Jared E. Knowles
##' @references Adapted from Simon Jackman from: \url{http://jackman.stanford.edu/classes/151B/06/class0517.r}
gammaGK <- function(x, y=NULL, print=FALSE){
  concordant <- function(x){ 
    ## get sum(matrix values > r AND > c) 
    ## for each matrix[r, c] 
    mat.lr <- function(r,c){ 
      lr <- x[(r.x > r) & (c.x > c)] 
      sum(lr) 
    } 
    
    ## get row and column index for each 
    ## matrix element 
    r.x <- row(x) 
    c.x <- col(x) 
    
    ## return the sum of each matrix[r, c] * sums 
    ## using mapply to sequence thru each matrix[r, c] 
    sum(x * mapply(mat.lr, r = r.x, c = c.x)) 
  } 
  
  discordant <- function(x){ 
    ## get sum(matrix values > r AND < c) 
    ## for each matrix[r, c] 
    mat.ll <- function(r,c){ 
      ll <- x[(r.x > r) & (c.x < c)] 
      sum(ll) 
    } 
    
    ## get row and column index for each 
    ## matrix element 
    r.x <- row(x) 
    c.x <- col(x) 
    
    ## return the sum of each matrix[r, c] * sums 
    ## using mapply to sequence thru each matrix[r, c] 
    sum(x * mapply(mat.ll, r = r.x, c = c.x)) 
  } 
  
  if(is.table(x) | is.matrix(x)){
    c <- concordant(x) 
    d <- discordant(x)
    n <- sum(x)
  }
  else{
    tab <- table(x,y)
    c <- concordant(tab) 
    d <- discordant(tab)
    n <- sum(tab)
  }
  gamma <- (c - d) / (c + d)
  
  arg <- (c+d)/(n*(1-(gamma^2)))
  stdError <- 1/sqrt(arg)
  z <- gamma/stdError
  if(print==TRUE){
    cat("Goodman-Kruskal gamma statistic:\n")
    cat(paste("Concordant Pairs",c,"\n"))
    cat(paste("Discordant Pairs",d,"\n\n"))
    cat(paste("Estimate of gamma:",
              signif(gamma,.Options$digits),
              "Standard error:",
              signif(stdError,.Options$digits),
              "\n\n"))
    
    cat(paste("H0: gamma = 0 vs HA: two-sided\n"))
    cat(paste("z:",
              signif(z, .Options$digits),
              "p-value:",
              signif(2*(1-pnorm(abs(z))), .Options$digits),
              "\n\n"))
    if(c<51 | d<51){
      cat("Warning: p-values are based on a normal approximation to the\n")
      cat("sampling distribution of the z test statistic, which is commonly\n")
      cat("considered to be good only if C and D are both > 50.\n")
    }
  }
  return(list(gamma = signif(gamma,.Options$digits), 
              se = signif(stdError,.Options$digits), 
              z =signif(z, .Options$digits), 
              sig = signif(2*(1-pnorm(abs(z))), .Options$digits)))
  invisible(NULL)
}

##' Cut out duplicate values in a dataframe
##'
##' Estimate the correlation between two unordered factor variables using the Goodman and Kruskal gamma statistic 
##' 
##' @param tmpdf a dataframe from which to cull
##' @return The dataframe with duplciated values culled from it
##' @note Yadda yadda yadda
##' @author Jared E. Knowles
cull <- function(tmpdf){
  tmpdf$concat_fw <- paste(tmpdf$Var1, tmpdf$Var2, sep="")
  tmpdf$concat_bw <- paste(tmpdf$Var2, tmpdf$Var1, sep="")
  
  ## Create a variable indicating the concatenation is unique
  tmpdf$unique <- 0
  
  ## Create an empty data frame to hold all the unique concatenations
  uniques <- as.data.frame(matrix(ncol=2, nrow=0))
  colnames(uniques)[1] <- "concat_fw"
  colnames(uniques)[2] <- "concat_bw"
  
  ## Loop through all the records in tmpdf
  for(i in 1:length(tmpdf$concat_fw)) {
    
    ## Isolate the current record's forward and backward concatenation
    temp <- tmpdf[i, c("concat_fw", "concat_bw")]
    
    ## If it hasn't been observed before, store it and set unique=1.
    ## Otherwise, ignore it and leave unique=0.
    if(!temp$concat_fw %in% uniques$concat_bw &
         !temp$concat_bw %in% uniques$concat_fw)
    {
      uniques <- rbind(uniques, temp)
      tmpdf$unique[i] <- 1
    }
    
  }
  
  ## Change Var1 and Var2 to NA where unique==0
  tmpdf$Var1[tmpdf$unique==0] <- NA
  tmpdf$Var2[tmpdf$unique==0] <- NA
  
  ## Clean up
  tmpdf$unique <- NULL; tmpdf$concat_fw <- NULL; tmpdf$concat_bw <- NULL;
  tmpdf <- na.omit(tmpdf)
  return(tmpdf)
}