#' Draws ecospace with species response surface for two dimensional simulated data
#' @author David Zeleny (zeleny.david@@gmail.com)
#' @param sampled.comm Object with simulated community data created by \code{sample.comm.2} function.
#' @param resolution Number of items drawn vertically and horizontally. Default = 200.
#' @param colors Vector of hexadecimal strings of colors used to draw individual species response surfaces. Each element should be 7 characters long, starting with \code{"#"} (e.g. \code{"#8DD3C7"}). Default is a vector of 12 colors created by pallete "Set3" from \code{library (RColorBrewer)}. 
#' @param species Vector with ID's of species drawn into plot. Default = NULL, which means all species are drawn. 
#' @param sort Should be species sorted before drawing, so as the species with largest area in the plot are drawn first? Default = TRUE.
#' @param asp Aspect ratio of the plotted region (see argument \code{asp} in \code{\link{plot.window}}).
#' @param sample.pch,sample.cex Plotting character and its size for points showing location of samples in ecospace. See \code{pch} and \code{cex} in function \code{\link{points}}.
#' @param species.pch,species.cex Plotting character and its size for points forming the species response surfase in the ecospace. See \code{pch} and \code{cex} in function \code{\link{points}}.
#' @examples
#' library (weimea)
#' draw.ecospace (sample.comm.2 (simul.comm.2 (totS = 100)), resolution = 100)
#' @rdname draw.ecospace
#' @export
draw.ecospace <- function (sampled.comm, resolution = 200, colors = NULL, species = NULL, sort = TRUE, asp = 1, sample.pch = 3, sample.cex = .5, species.pch = 16, species.cex = 1)
{
  if (is.null (colors)) colors <- RColorBrewer::brewer.pal (12, 'Set3')
  if (is.null (species)) species <- 1:sampled.comm$simul.comm$totS
  sc <- sampled.comm$simul.comm
  #windows ()
  plot.new ()
  plot.window (xlim = c(0, sc$gr1.length), ylim = c(0, sc$gr2.length), asp = asp)
  axis (1)
  axis (2)
  box ()
  gr1gr2 <- expand.grid (gr1 = round (seq (1, sc$gr1.length, length = resolution)), gr2 = round (seq (1, sc$gr2.length, length = resolution)))
  spec.prob <- sapply (1:ncol (sc$A1.all), FUN = function (sp) sc$A1.all[gr1gr2$gr1,sp]*sc$A2.all[gr1gr2$gr2, sp])
  sorted.species <- order (colSums (spec.prob[,species] > 0), decreasing = T)
  colors.all <- rep (colors, length.out = length (species))
  for (sp in sorted.species)
  {
    if (sum (spec.prob[,sp] > 0)>0)
    {
      lower.bound <- quantile (spec.prob[spec.prob[,sp]>0,sp], prob = .1)
      points.coords <- gr1gr2[spec.prob[,sp] > lower.bound,]
      species.prob <- spec.prob[spec.prob[,sp] > lower.bound,sp]
      species.prob.ranged <- (species.prob - min (species.prob))/(max (species.prob) - min (species.prob))*255
      colors.temp <- paste (colors.all[sp], format (as.hexmode (round (species.prob.ranged)), upper = T), sep = '')
      points (points.coords, col = colors.temp, pch = species.pch, cex = species.cex) 
    }
    
  }
  points (sampled.comm$sample.x1, sampled.comm$sample.x2, pch = sample.pch, cex = sample.cex)
}
