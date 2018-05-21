#' Fits community weighted mean of species attributes (CWM) onto an ordination
#' 
#' The function fits community weighted mean of species attributes (e.g. traits), calculated for individual samples, onto an ordination, using multiple regression. The relationship is tested by permutation test (defaults to modified permutation test sensu Zeleny & Schaffers 2012).
#' @param ord An ordination object or other structure from which the ordination scores can be extracted (including a data frame or matrix of scores).
#' @param cwm Community weighted means calculated by function \code{cwm}; object of the class \code{cwm}. 
#' @param wstand Logical; should be the cwm calculated on weighted standardized species attributes? Defaults to the same setting as used in \code{cwm} function (FALSE).
#' @param permutations Number of permutations for the permutation test.
#' @param test Permutation test.
#' @param choices Axes to be plotted.
#' @param display In fitting functions these are ordinary site scores or linear combination scores ("lc") in constrained ordination (cca, rda, capscale). In scores function they are either "vectors" or "factors" (with synonyms "bp" or "cn", resp.). (Description from vegan::envfit help file).
#' @param w Weights used in fitting (concerns mainly \code{\link{cca}} and \code{\link{decorana}} results which have nonconstant weights). (Description from vegan::envfit help file).
#' @param na.rm Remove rows in ordination scores or \code{cwm} variables with missing values - the whole row of data in both ordination and \code{cwm} is removed in case of any missing values and \code{na.rm = TRUE}.
#' @param adjustP Logical; should the P-values be adjusted for multiple testing?
#' @param p.adjust.method Character string; method of P-value adjustement (check \code{\link{p.adjust.methods}} for available option, default \code{'holm'}).
#' @param p.adjust.n Integer; what is the number of comparisons used to adjust the P-values? Default (\code{p.adjust.n = NULL}) means that the number equals to the number of conducted tests (= number of cwm variables).
#' @param ... Parameters passed to \code{\link{scores}}.
#' @details
#' The function is derived from the function \code{\link{envfit}} in \code{library (vegan)}, which fits the environmental vectors or factors onto an ordination and tests their significance using permutation test. The modification is in the permutation schema used to test the fit significance: while the original \code{vegan::\link{envfit}} function builds the null distribution of test statistic (r2) by permuting the values of environmental vector (analogy of \code{standard} or \code{rowbased} test), the function \code{weimea::envfit_cwm} additionally offers also \code{modified} or \code{rowbased} test, based on permuting species attributes prior to calculation of \code{cwm}, and \code{max} test combining results of \code{rowbased} and \code{colbased} as \code{P_max = max (P_row, P_col)}. While the null hypothesis of \code{standard} (\code{rowbased}) permutation test that there is "no statistical relationship between mean indicator values and sample ordination scores" is easy to reject in case that CWM's are numerically derived from the same species composition matrix as sample ordination scores (this applies even if some kind of transformation, e.g. log or presence-absence, was applied on the species data before calculating either ordination or CWM), since such relationship is spurious. The null hypothesis of the alternative \code{modified} (\code{rowbased}) permutation test is that there is "no statistical relationship between information in species attributes and sample ordination scores", which is not spurious if species attributes (\code{traits}) are extrinsic to species composition (i.e. not numerically derived from it). The third, \code{max} test, combines these two results together, and is probably not meaningful in the context of relating CWM of species attributes and sample ordination scores, because it includes the \code{standard} (\code{rowbased}) test which has no meaningful null hypothesis (as explained above).
#'
#'The permutation test requires information about the original vegetation matrix used to calculate the weighted mean of species indicator values, and also original species attributes, which is stored in the object of \code{cwm} class as \code{'com'} and \code{'traits'} attributes. If the argument \code{cwm} is not the object of \code{cwm} class, the function will calculate standard \code{vegan::envfit} and return warning message.
#'
#' The generic functions \code{print}, \code{plot} and \code{scores} for the class envfit can be used to print the results, plot the vectors of CWM onto an ordination diagram and to return the fitting results. Check \code{?envfit} for details which arguments can be used in these functions.
#' @importFrom stats complete.cases cor cor.test p.adjust p.adjust.methods symnum weights
#' @importFrom vegan envfit scores
#' @export
envfit_cwm <- function (ord, cwm, wstand = attr (cwm, 'wstand'), permutations = 999, choices = c(1, 2), display = "sites", w = weights(ord, display), na.rm = FALSE, test = 'modified', adjustP = FALSE, p.adjust.method = 'holm', p.adjust.n = NULL, ...) 
{
  if (!is.cwm (cwm)){
    warning ("The argument cwm is not of 'cwm' class - vegan::envfit function applied instead.")
    return (envfit (ord = ord, env = cwm, permutations = permutations, choices = choices, display = display, w = w, na.rm = na.rm, ...))
  } else {
    TEST <- c('none', 'parametric', 'standard', 'rowbased', 'modified', 'colbased', 'max', 'all')
    test <- match.arg (test, TEST, several.ok = F)
    weights.default <- function(object, ...) NULL
    w <- eval(w)  # originally w < eval(w)
    wstand <- eval (wstand)
    vectors <- NULL
    X <- scores(ord, display = display, choices = choices, ...)
    keep <- complete.cases(X) & complete.cases(cwm)
    if (any(!keep)) {
      if (!na.rm) 
        stop("missing values in data: consider na.rm = TRUE")
      X <- X[keep, , drop = FALSE]
      cwm <- cwm[keep, , drop = FALSE]
      w <- w[keep]
      na.action <- structure(seq_along(keep)[!keep], class = "omit")
    }
    nr <- nrow (X)
    if (missing(w) || is.null(w)) 
      w <- 1
    if (length(w) == 1) 
      w <- rep(w, nrow(X))
    if (nrow(cwm) != nrow(X)) 
      stop("input data have non-matching numbers of observations")
    out.vectorfit <- vectorfit_cwmCpp(X = as.matrix (X), L = as.matrix (extract.cwm (cwm, 'com')), t = as.matrix (extract.cwm (cwm, 'traits')), wstand = wstand, perm = permutations, w = w, test = test)
    
    heads <- t(out.vectorfit$heads)
    if (is.null(colnames(X))) 
      colnames(heads) <- paste("Dim", 1:ncol (X), sep = "")
    else colnames(heads) <- colnames(X)
    rownames (heads) <- colnames (cwm)
    pvals <- as.vector (out.vectorfit$P)
    if (adjustP) pvals <- p.adjust (pvals, method = p.adjust.method, n = if (is.null (p.adjust.n)) length (pvals) else p.adjust.n)
    vectors <- list (arrows = heads, r = as.vector (out.vectorfit$r), permutations = permutations, pvals = pvals)
    class (vectors) = "vectorfit"
    sol <- list(vectors = vectors, factors = NULL)
    if (!is.null(na.action)) 
      sol$na.action <- na.action
    class(sol) <- "envfit"
    return (sol)
  }
}

