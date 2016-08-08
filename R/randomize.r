#' Weighted mean calculated from randomized species attributes
#' 
#' Function applied on the object of class \code{wm}, which returns values of weighted mean calculated from randomized species attributes.
#' @param x object of the class \code{wm}
#' @param permutations number of randomizations
#' @param FUN function to be applied on the column of calculated weighted mean values
#' @param progress.bar logical value, should be the progress bar indicating the progress of the calculation launched?
#' @param parallel integer or NULL (default). If integer, calculation will be conducted in parallel on number of cores given by \code{parallel}
#' @param library.par character vector with libraries needed for application of function defined by \code{FUN} argument (to be exported into parallel process). Not necessary if the functions in \code{FUN} are stated explicitly (e.g. \code{vegan:::cca})
#' @param export.cl ??
#' @param ... parameters passed to other functions.
#' @author David Zeleny (zeleny.david@@gmail.com)
#' @importFrom utils setWinProgressBar winProgressBar
#' @export
randomize <- function (x, ...) UseMethod ('randomize')

#' @rdname randomize
#' @export
randomize.wm <- function (x, permutations = 1, FUN = function (y) y, progress.bar = F, parallel = NULL, library.par = NULL, export.cl = NULL, ...)
{
  if (!is.wm (x)) stop ("Object x is not of class 'wm'")
  if (progress.bar & is.null (parallel)) win.pb <- winProgressBar(title = "Permutation progress bar", label = "", min = 0, max = permutations, initial = 0, width = 300)
  sitspe <- attr (x, 'sitspe')
  speatt <- attr (x, 'speatt')
  wm.rand <- function (sitspe, speatt)
  {
    speatt.rand <- apply (speatt, 2, FUN = function (x) {pointer <- which (!is.na (x)); x[pointer] <- sample (x[pointer]); x})
    wm (sitspe, speatt.rand) 
  }
  if (is.null (parallel))
  {
    temp.result <- list ()
    for (perm in seq (1, permutations))
    {
      if (progress.bar) setWinProgressBar (win.pb, perm)
      temp.result[[perm]] <- lapply (list (wm.rand (sitspe = sitspe, speatt = speatt)), FUN = FUN)[[1]]
    }
  }  
  
  if (!is.null (parallel))
  {
    #require (parallel)
    cl <- parallel::makeCluster(parallel)
    parallel::clusterExport (cl, varlist = c("FUN", "speatt", "sitspe", "library.par", "wm.rand", "wm"), envir = environment ())
    #parallel::clusterExport (cl, eval (call ('library', 'weimea')))
    if (!is.null (export.cl)) parallel::clusterExport (cl, export.cl)
    if (!is.null (library.par)) parallel::clusterEvalQ (cl, eval (call ('library', library.par)))
    temp.result <- parallel::parLapply (cl, seq (1, permutations), fun = function (x)
    {
      lapply (list (wm.rand (sitspe, speatt)), FUN = FUN)[[1]]
    })
    parallel::stopCluster (cl)
  }
  if (progress.bar & is.null (parallel)) close (win.pb)
  if (permutations == 1) temp.result <- temp.result[[1]]
  return (temp.result)
}
