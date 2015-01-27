#' Changes the range of species attributes
#' @param object object of the class \code{wm}
#' @param max upper bound of the rescaled values
#' @param recalculate should the mean values in \code{object} be recalculated?
#' @export

rangeWM <- function (object, max = 9, recalculate = T)
{
  if (!is.wm (object)) stop ("object is not of the class 'wm'")
  speatt <- attr (object, 'speatt')
  res <- ceiling((speatt - min (speatt))/(diff (range (speatt))/max))
  res[res == 0] <- 1
  if (recalculate) object <- wm (speatt = res, sitspe = attr (object, 'sitspe')) else attr (object, 'speatt') <- res
  object
}