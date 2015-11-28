#' Weighted mean of species attributes
#' 
#' Function \code{wm} calculates weighted mean of species attributes, using matrices of species composition and species attributes. Other functions are for handling the objects of class \code{wm}.
#' 
#' @param sitspe Matrix or data.frame with community data matrix (sites x species)
#' @param speatt Vector, matrix or data.frame with species attributes (species x attributes)
#' @param object,x Object of the class \code{wm}
#' @param drop  In function \code{[.wm]} if \code{TRUE} coerces the result to the lowest possible dimension (i.e. vector if matrix has only one column). Currently not implemented, change into \code{TRUE} will have no effect.
#' @param i,j Subscripts of the matrix of the class "Wm" (rows and columns).
#' @param ... Other arguments passed into \code{as.matrix}, \code{summary} or \code{print} function. Currently not supported.
#' @details
#' Function \code{[.wm]} is for extracting specified rows and columns from matrix of class \code{wm}. As a side effect, resulting object will have concatenated the \code{sitspe} and \code{speatt} attributes to match the dimension of the resulting matrix. This function is only for extracting the parts of \code{wm} object, not for replacing! (attempt to replace will work, but will break the inner structure).
#' 
#' Generic functions \code{sitspe} and \code{speatt} applied on object of class \code{"wm"} extracts the original species composition matrix and species attribute matrix, respectively, which were used to calculate weighted means.
#' 
#' @return Object of class \code{"wm"}, which has \code{print}, \code{summary}, \code{as.matrix} and \code{'[.'} methods. Object of \code{"wm"} class contains the matrix of calculated weighted means of species attributes for each sample (sample x weighted mean) and two attributes: \code{sitspe} species x sample matrix from which the weighted mean was calculated, and \code{speatt} species x attributes matrix with species attributes. All weighted means of species attributes must be based on the same species x sample matrix with the same number of samples.
#' @examples
#' # Calculation of weighted mean of species Ellenberg indicator values using dataset Vltava
#' data (vltava)
#' mean.eiv <- wm (sitspe = vltava$spe, speatt = vltava$ell)
#' 
#' summary (mean.eiv)
#' 
#' # Extracting values from the object of \code{wm} class
#' mean.eiv[,1]
#' mean.eiv[1:10, 2:3]
#' @author David Zeleny (zeleny.david@@gmail.com)
#' @useDynLib weimea
#' @importFrom Rcpp evalCpp
#' @import RcppArmadillo
#' @export
wm <- function (sitspe, speatt)
{
  sitspe <- as.matrix (sitspe)
  speatt <- as.matrix (speatt)
  if (is.null (colnames (speatt)))  colnames (speatt) <- paste ('speatt', 1:ncol (speatt), sep = '_')
  if (any (is.na(colnames(speatt)))) colnames (speatt)[is.na (colnames (speatt))] <- paste ('speatt', seq (1, sum (is.na (colnames (speatt)))), sep = '_')
  wm.temp <- wm_rcpp (sitspe, speatt)   # Rcpp version of wm function
  colnames (wm.temp) <- colnames (speatt)
  attr (wm.temp, 'sitspe') <- as.data.frame (sitspe)
  attr (wm.temp, 'speatt') <- as.data.frame (speatt)
  attr(wm.temp, 'class') <- c('wm')
  wm.temp
}

#' @export
wmR <- function (sitspe, speatt)
{
  sitspe <- as.matrix (sitspe)
  speatt <- as.data.frame (as.matrix (speatt))
  wm.temp <- apply (speatt, 2, FUN = function (x) vegan::decostand (sitspe[,!is.na (x)], 'total') %*% x[!is.na(x)])
  attr (wm.temp, 'sitspe') <- sitspe
  attr (wm.temp, 'speatt') <- speatt
  attr(wm.temp, 'class') <- c('wm')
  wm.temp
}

#' @rdname wm
#' @export
is.wm <- function (object)
{
  if (any (class (object) == 'wm') & !is.null (attr (object, 'sitspe'))  & !is.null (attr (object, 'speatt'))) TRUE else FALSE
}


#' @rdname wm
#' @export
"[.wm" <- function (object, i, j, drop = F)
{
  if (missing (i)) i <- 1:nrow (object)
  if (missing (j)) j <- 1:ncol (object)
  sitspe <- attr (object, 'sitspe')
  speatt <- attr (object, 'speatt')
  object <- as.matrix (object)
  res <- object[i, j, drop = F]
  attr (res, 'sitspe') <- sitspe [i,, drop = F]
  attr (res, 'speatt') <- speatt [,j, drop = F]
  class (res) <- c('wm')
  res
}

#' @rdname wm
#' @export
as.matrix.wm <- function (x, ...)
{
  attr (x, 'sitspe') <- NULL
  attr (x, 'speatt') <- NULL
  x <- unclass (x)
  x <- as.matrix (x)
  return (x)
}

#' @rdname wm
#' @param long should summary return long output? (TRUE vs FALSE)
#' @export
summary.wm <- function (object, long = F, ...)
{
  sitspe <- attr (object, 'sitspe')
  speatt <- attr (object, 'speatt')
  cat ("Object of the class 'wm'\n")
  cat ("\nWeighted means:     \t\t", dim (object), "\t(sites x variables)")
  cat ("\nSpecies composition:\t\t", dim (sitspe), "\t(sites x species), \trange of values:", range (sitspe, na.rm = T))
  cat ("\nSpecies attributes: \t\t", dim (speatt), "\t(species x attributes)")
  na <- apply (speatt, 2, FUN = function (x) sum (is.na (x)))
  if (sum (na) > 0)  cat ("\n\t\tMissing values of sp. attributes:\t", paste (names (na), na, sep = ' ', col = '\t'))
  if (long)
  {
    cat ("\n\nSummary of wm matrix:\n\n")
    print (summary (as.matrix (object)))
    cat ("\nSummary of community matrix\n\n")
    print (summary (attr (object, 'sitspe')))
    cat ("\nSummary of species attributes\n\n")
    print (summary (attr (object, 'speatt')))
  }
}

#' @rdname wm
#' @export
print.wm <- function (x, ...)
{
  attr (x, 'sitspe') <- NULL
  attr (x, 'speatt') <- NULL
  class (x) <- 'matrix'
  print (as.matrix (x))
}
