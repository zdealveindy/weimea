#' Weighted mean calculated from randomized species attributes
#' 
#' Function applied on the object of class \code{cwm}, which returns values of weighted mean calculated from randomized species attributes.
#' @param x object of the class \code{cwm}
#' @param permutations number of randomizations
#' @param FUN function to be applied on the column of calculated weighted mean values
#' @param parallel integer or NULL (default). If integer, calculation will be conducted in parallel on number of cores given by \code{parallel}
#' @param library.par character vector with libraries needed for application of function defined by \code{FUN} argument (to be exported into parallel process). Not necessary if the functions in \code{FUN} are stated explicitly (e.g. \code{vegan:::cca})
#' @param export.cl ??
#' @param ... parameters passed to other functions.
#' @author David Zeleny (zeleny.david@@gmail.com)
#' @export
randomize <- function (x, ...) UseMethod ('randomize')

#' @rdname randomize
#' @export
randomize.cwm <- function (x, permutations = 1, FUN = function (y) y, parallel = NULL, library.par = NULL, export.cl = NULL, ...)
{
  if (!is.cwm (x)) stop ("Object x is not of class 'cwm'")
  sitspe <- attr (x, 'sitspe')
  speatt <- attr (x, 'speatt')
  cwm.rand <- function (sitspe, speatt)
  {
    speatt.rand <- apply (speatt, 2, FUN = function (x) {pointer <- which (!is.na (x)); x[pointer] <- sample (x[pointer]); x})
    cwm (sitspe, speatt.rand) 
  }
  if (is.null (parallel))
  {
    temp.result <- list ()
    for (perm in seq (1, permutations))
    {
      temp.result[[perm]] <- lapply (list (cwm.rand (sitspe = sitspe, speatt = speatt)), FUN = FUN)[[1]]
    }
  }  
  
  if (!is.null (parallel))
  {
    #require (parallel)
    cl <- parallel::makeCluster(parallel)
    parallel::clusterExport (cl, varlist = c("FUN", "speatt", "sitspe", "library.par", "cwm.rand", "cwm"), envir = environment ())
    #parallel::clusterExport (cl, eval (call ('library', 'weimea')))
    if (!is.null (export.cl)) parallel::clusterExport (cl, export.cl)
    if (!is.null (library.par)) parallel::clusterEvalQ (cl, eval (call ('library', library.par)))
    temp.result <- parallel::parLapply (cl, seq (1, permutations), fun = function (x)
    {
      lapply (list (cwm.rand (sitspe, speatt)), FUN = FUN)[[1]]
    })
    parallel::stopCluster (cl)
  }
  if (permutations == 1) temp.result <- temp.result[[1]]
  return (temp.result)
}
