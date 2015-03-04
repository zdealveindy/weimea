#' Whittaker's index of association
#' 
#' Creates distance matrix based on Whittaker's index of association, which reflects the dissimilarities between weighted means of species attributes calculated for individual samples, which can be attributed to differences in species composition. The metric is Manhattan-type distance calculated on species profiles (i.e. on community matrix with rows standardized to row-totals).
#' @param sitspe Compositional matrix (samples x species).
#' @author David Zeleny (zeleny.david@@gmail.com)
#' @references
#' Legendre P. & Legendre L. 2012. Numerical Ecology. 3rd edn. Elsevier, Oxford, UK. 
#' @export
ia <- function (sitspe) vegan::vegdist (vegan::decostand (sitspe, 'total'), 'manhattan')/2

#' @rdname ia
#' @export
betadiv.ia <- function (sitspe) sum (ia (sitspe)^2)/(nrow (sitspe)*(nrow (sitspe)-1))*2

#' @rdname ia
#' @export
whittaker <- function (sitspe) sum (colSums (sitspe) > 0)/ mean (rowSums (sitspe > 0))
