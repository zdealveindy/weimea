#' Community weighted mean of species attributes
#' 
#' Function \code{cwm} calculates weighted mean of species attributes, using matrix of species composition and one or several species attributes.
#' 
#' @param com Matrix or data frame with species composition data (samples x species).
#' @param traits Vector, matrix or data frame with species attributes (species x attributes). This can be \code{numeric} or \code{factor}.
#' @param wstand Logical. Should be the values in \code{traits} weighted-standardized prior to calculation of CWM? Weights are column sums in community data matrix \code{com}.
#' @param object,x Object of the class \code{cwm}
#' @param i,j Subscripts of the matrix of the class "cwm" (rows and columns).
#' @param drop Argument in the subsetting function, but currently not implemented (\code{[.cwm} always return data frame with at least one variable).
#' @param ... Other arguments passed into \code{as.matrix}, \code{summary} or \code{print} function. Currently not supported.
#' @details
#' Function \code{[.cwm]} is for extracting specified rows and columns from matrix of class \code{cwm}. As a side effect, resulting object will have concatenated the \code{com} and \code{traits} attributes to match the dimension of the resulting matrix. This function is only for extracting the parts of \code{cwm} object, not for replacing - the attempt to replace will work, but will break the inner structure of the object.
#' 
#' Generic function \code{extract} with argument \code{what} applied on object of class \code{"cwm"} extracts the original species composition matrix or species attribute matrix, respectively, which were used to calculate weighted means.
#' 
#' Function \code{range_cwm} changes the range of species attributes and recalculates the CWM values.
#' @return Object of class \code{"cwm"}, which has \code{print}, \code{summary}, \code{as.matrix} and \code{'[.'} methods. Object of \code{"cwm"} class contains the matrix of calculated weighted means of species attributes for each sample (sample x weighted mean) and three attributes: \code{com}, species x sample matrix from which the weighted mean was calculated, \code{traits}, species x attributes matrix with species attributes, and \code{wstand}, logical value specifying whether traits were weighted-standardized prior to calculation of CWM (weights are column sums of \code{com}). All weighted means of species attributes must be based on the same species x sample matrix with the same number of samples.
#' @examples
#' # Calculation of weighted mean of species Ellenberg indicator values using dataset Vltava
#' data (vltava)
#' mean.eiv <- cwm (com = vltava$spe, traits = vltava$eiv)
#' 
#' summary (mean.eiv)
#' 
#' # Extracting values from the object of cwm class
#' mean.eiv[,1]
#' mean.eiv[1:10, 2:3]
#' @author David Zeleny (zeleny.david@@gmail.com)
#' @useDynLib weimea
#' @importFrom Rcpp evalCpp
#' @export
cwm <- function (com, traits, wstand = FALSE)
{
   
  dummy <- function(df) {  
    as.dummy <- function (var) 
    {
      res <- model.matrix (~ as.matrix  (var) - 1)
      lev <- levels (var)
      colnames (res) <- lev
      return (res)
    }
    temp_res <- lapply (df, FUN = function (column) if (is.factor (column)) as.dummy (column) else column)
    res <- do.call (cbind.data.frame, temp_res)
    return (res)
  }

  com <- as.matrix (com)
  traits <- as.matrix (dummy (as.data.frame (traits)))
  if (ncol (com) != nrow (traits)) stop ("The number of species in 'traits' does not match the number of species in 'com'!")
  if (is.null (colnames (traits))) colnames (traits) <- paste ('trait', 1:ncol (traits), sep = '_')
  if (any (is.na(colnames(traits)))) colnames (traits)[is.na (colnames (traits))] <- paste ('trait', seq (1, sum (is.na (colnames (traits)))), sep = '_')
  cwm.temp <- as.data.frame (cwmCpp (com, traits, wstand))   # Rcpp version of cwm function
  names (cwm.temp) <- colnames (traits)
  attr (cwm.temp, 'class') <- c('cwm', 'data.frame')
  attr (cwm.temp, 'com') <- as.data.frame (com)
  attr (cwm.temp, 'traits') <- as.data.frame (traits)
  attr (cwm.temp, 'wstand') <- wstand

  return (cwm.temp)
}


#' @rdname cwm
#' @export
is.cwm <- function (object)
{
  if (any (class (object) == 'cwm') & !is.null (attr (object, 'com'))  & !is.null (attr (object, 'traits'))) TRUE else FALSE
}


#' @rdname cwm
#' @export
"[.cwm" <- function (object, i, j, drop = FALSE)
{
  if (drop) warning ("The argument 'drop' will be ignored (the object of 'cwm' class is always data frame).")
  if (missing (i)) i <- seq_len (nrow (object))
  if (missing (j)) j <- seq_len (ncol (object))
  com <- attr (object, 'com')
  traits <- attr (object, 'traits')
  class (object) <- 'data.frame'
  res <- object[i, j, drop = FALSE]
  attr (res, 'com') <- com [i,, drop = F]
  attr (res, 'traits') <- traits [,j, drop = F]
  class (res) <- c('cwm', 'data.frame')
  res
}

#' #' @rdname cwm
#' #' @export
#' as.matrix.cwm <- function (x, ...)
#' {
#'   attr (x, 'com') <- NULL
#'   attr (x, 'traits') <- NULL
#'   attr (x, 'wstand') <- NULL
#'   x <- as.matrix (unclass (x))
#'   return (x)
#' }

#' @rdname cwm
#' @param long should summary return long output? (TRUE vs FALSE)
#' @export
summary.cwm <- function (object, long = FALSE, ...)
{
  com <- attr (object, 'com')
  traits <- attr (object, 'traits')
  cat ("Object of the class 'cwm'\n")
  cat ("\nWeighted means:     \t\t", dim (object), "\t(samples x traits)")
  cat ("\nSpecies composition:\t\t", dim (com), "\t(samples x species), \trange of values:", range (com, na.rm = T))
  cat ("\nSpecies attributes: \t\t", dim (traits), "\t(species x traits)")
  na <- apply (traits, 2, FUN = function (x) sum (is.na (x)))
  if (sum (na) > 0)  cat ("\n\t\tMissing values of sp. attributes:\t", paste (names (na), na, sep = ' ', col = '\t'))
  if (long)
  {
    cat ("\n\nSummary of cwm matrix:\n\n")
    print (summary (as.matrix (object)))
    cat ("\nSummary of species composition matrix\n\n")
    print (summary (attr (object, 'com')))
    cat ("\nSummary of species attributes\n\n")
    print (summary (attr (object, 'traits')))
  }
}

#' @rdname cwm
#' @export
print.cwm <- function (x, ...)
{
  print (as.matrix (x))
}

#' @rdname cwm
#' @export
extract <- function (x, ...) UseMethod ('extract')

#' @rdname cwm
#' @export
#' @param what Attributes extracted from the object of class \code{cwm}; either \code{traits} for the matrix of species attributes (species x traits), or \code{com} for matrix of species composition (samples x species).
#' 
extract.cwm <- function (x, what = 'traits', ...)
{
  WHAT <- c('traits', 'com')
  what <- match.arg (what, WHAT)
  if (what == 'traits') x.out <- attr (x, 'traits')
  if (what == 'com') x.out <- attr (x, 'com')
  return (x.out)
}

#' @rdname cwm
#' @param max upper bound of the rescaled values
#' @param recalculate should the mean values in \code{object} be recalculated?
#' @export
range_cwm <- function (object, max = 9, recalculate = TRUE)
{
  if (!is.cwm (object)) stop ("object is not of the class 'cwm'")
  traits <- attr (object, 'traits')
  res <- ceiling((traits - min (traits))/(diff (range (traits))/max))
  res[res == 0] <- 1
  if (recalculate) object <- cwm (com = attr (object, 'com'), traits = res) else attr (object, 'traits') <- res
  object
}

#' #' @rdname cwm
#' #' @export
#' `[[.cwm` <- function (x, i, exact = TRUE) 
#' {
#'   na <- nargs() - !missing(exact)
#'   if (!all(names(sys.call()) %in% c("", "exact"))) 
#'     warning("named arguments other than 'exact' are discouraged")
#'   out <- as.data.frame (.subset2(x, i, exact = exact), col.names = i)
#'   attr (out, 'com') <- attr (x, 'com')
#'   attr (out, 'traits') <- attr (x, 'traits')
#'   attr (out, 'wstand') <- attr (x, 'wstand')
#'   attr (out, 'class') <- c('cwm')
#'   return (out)
#' }
#' 
#' #' @rdname cwm
#' #' @export
#' `$.cwm` <- function (x, name) 
#' {
#'   a <- x[[name]]
#'   if (!is.null(a)) 
#'     return(a)
#'   a <- x[[name, exact = FALSE]]
#'   if (!is.null(a) && getOption("warnPartialMatchDollar", default = FALSE)) {
#'     names <- names(x)
#'     warning(gettextf("Partial match of '%s' to '%s' in data frame", 
#'                      name, names[pmatch(name, names)]))
#'   }
#'   return(a)
#' }