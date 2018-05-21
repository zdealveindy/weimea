#' Weighted mean of sample attributes - species niche centroids (SNC)
#' 
#' NEEDS TO REWRITE!!!
#' Function \code{snc} calculates weighted mean of sample attributes, using matrix of species composition (\code{com}) and one or several sample attributes (\code{env}). 
#' 
#' @param com Matrix or data.frame with species composition matrix (samples x species)
#' @param env Vector, matrix or data.frame with env. variables (samples x env). This can be \code{numeric} or \code{factor}.
#' @param wstand Logical. Should be the values in \code{env} weighted-standardized prior to calculation of SNC? Weights are row sums in species composition matrix \code{com}.
#' @param object,x Object of the class \code{snc}
#' @param i,j Subscripts of the matrix of the class "snc" (rows and columns).
#' @param ... Other arguments passed into \code{as.matrix}, \code{summary} or \code{print} function. Currently not supported.
#' @details
#' Function \code{[.snc]} is for extracting specified rows and columns from matrix of class \code{snc}. As a side effect, resulting object will have concatenated the \code{com} and \code{env} attributes to match the dimension of the resulting matrix. This function is only for extracting the parts of \code{snc} object, not for replacing; attempt to replace will work, but will break the inner structure.
#' 
#' Generic function \code{extract} with argument \code{what} applied on object of class \code{"snc"} extracts the original species composition matrix or sample attribute matrix, respectively, which were used to calculate weighted means.
#' 
#' @return Object of class \code{"snc"}, which has \code{print}, \code{summary}, \code{as.matrix} and \code{'[.'} methods. Object of \code{"snc"} class contains the matrix of calculated weighted means of sample attributes (species niche centroids) for each species (species x snc) and three attributes: \code{com}, species x sample matrix from which the weighted mean was calculated, \code{env}, samples x env. variables matrix with sample attributes, and \code{wstand}, logical value specifying whether env.variables were weighted-standardized prior to calculation of SNC (weights are row sums of \code{com}). All weighted means of sample attributes must be based on the same species x sample matrix with the same number of species.
#' @examples
#' # Calculation of weighted mean of species Ellenberg indicator values using dataset Vltava
#' data (vltava)
#' mean.env <- snc(com = vltava$spe, env = vltava$env[, c('ELEVATION', 'pH')])
#' 
#' summary (mean.env)
#' 
#' # Extracting values from the object of \code{snc} class
#' mean.env[,1]
#' mean.env[1:10, 1]
#' @author David Zeleny (zeleny.david@@gmail.com)
#' @useDynLib weimea
#' @importFrom Rcpp evalCpp
#' @import RcppArmadillo
#' @export
snc <- function (com, env, wstand = FALSE)
{
  dummy <- function(df) {  
    NUM <- function(dataframe)dataframe[,sapply(dataframe,is.numeric), drop = F]
    FAC <- function(dataframe)dataframe[,sapply(dataframe,is.factor), drop = F]
    if (is.null(ncol(FAC(df))) || ncol(FAC(df)) == 0)
      DF <- df else {
        if (is.null(ncol(NUM(df))) || ncol(NUM(df)) == 0) {
          DF <- data.frame(NUM(df), ade4::acm.disjonctif(FAC(df)))
          names(DF)[1] <- colnames(df)[which(sapply(df, is.numeric))]
        } else {
          DF <- data.frame(NUM(df), ade4::acm.disjonctif(FAC(df)))
        }
      }
    return(DF)
  } 
  com <- as.matrix (com)
  env <- as.matrix (dummy (as.data.frame (env)))
  if (is.null (colnames (env))) colnames (env) <- paste ('env', 1:ncol (env), sep = '_')
  if (any (is.na(colnames(env)))) colnames (env)[is.na (colnames (env))] <- paste ('env', seq (1, sum (is.na (colnames (env)))), sep = '_')
  snc.temp <- as.data.frame (sncCpp (com, env, wstand))   # Rcpp version of snc function
  names (snc.temp) <- colnames (env)
  attr (snc.temp, 'com') <- as.data.frame (com)
  attr (snc.temp, 'env') <- as.data.frame (env)
  attr (snc.temp, 'wstand') <- wstand
  attr(snc.temp, 'class') <- c('snc', 'data.frame')
  return (snc.temp)
}


#' @rdname snc
#' @export
is.snc <- function (object)
{
  if (any (class (object) == 'snc') & !is.null (attr (object, 'com'))  & !is.null (attr (object, 'env'))) TRUE else FALSE
}


#' @rdname snc
#' @export
"[.snc" <- function (object, i, j)
{
  if (missing (i)) i <- seq_len (nrow (object))
  if (missing (j)) j <- seq_len (ncol (object))
  com <- attr (object, 'com')
  env <- attr (object, 'env')
  class (object) <- 'data.frame'
  res <- object[i, j, drop = F]
  attr (res, 'com') <- com [i,, drop = F]
  attr (res, 'env') <- env [,j, drop = F]
  class (res) <- c('snc', 'data.frame')
  res
}

#' @rdname snc
#' @export
as.matrix.snc <- function (x, ...)
{
  attr (x, 'com') <- NULL
  attr (x, 'env') <- NULL
  class (x) <- 'data.frame'
  x <- as.matrix (x)
  return (x)
}

#' @rdname snc
#' @param long should summary return long output? (TRUE vs FALSE)
#' @export
summary.snc <- function (object, long = F, ...)
{
  com <- attr (object, 'com')
  env <- attr (object, 'env')
  cat ("Object of the class 'snc'\n")
  cat ("\nWeighted means:     \t\t", dim (object), "\t(species x env)")
  cat ("\nSpecies composition:\t\t", dim (com), "\t(samples x species), \trange of values:", range (com, na.rm = T))
  cat ("\nSample attributes: \t\t", dim (env), "\t(samples x env)")
  na <- apply (env, 2, FUN = function (x) sum (is.na (x)))
  if (sum (na) > 0)  cat ("\n\t\tMissing values of env:\t", paste (names (na), na, sep = ' ', col = '\t'))
  if (long)
  {
    cat ("\n\nSummary of snc matrix:\n\n")
    print (summary (as.matrix (object)))
    cat ("\nSummary of community matrix\n\n")
    print (summary (attr (object, 'com')))
    cat ("\nSummary of sample attributes\n\n")
    print (summary (attr (object, 'env')))
  }
}

#' @rdname snc
#' @export
print.snc <- function (x, ...)
{
  print (as.matrix (x))
}

#' @rdname snc
#' @param what Attributes extracted from the object of class \code{snc}; either \code{env} for the matrix of sample attributes (samples x env), or \code{com} for matrix of species composition (samples x species).
#' @export
extract.snc <- function (x, what = 'env', ...)
{
  WHAT <- c('env', 'com')
  what <- match.arg (what, WHAT)
  if (what == 'env') x.out <- attr (x, 'env')
  if (what == 'com') x.out <- attr (x, 'com')
  return (x.out)
}