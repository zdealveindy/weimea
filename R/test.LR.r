#' Test of the link between sample attributes and matrix of species composition (i.e. whether \bold{L} is linked to \bold{R})
#' 
#' The function testing the link between sample attributes and species composition of the matrix, from which weighted means are calculated. The test is based either on db-RDA (distance-based redundancy analysis), or Moran's I. Significant relationship is considered as an argument to use modified permutation test instead of the standard permutation test for testing the relationship between weighted mean of species attributes and sample attributes.
#' @param M Object of the 'wm' class. Matrix with weighted means of species attributes.
#' @param env Matrix with environmental variables.
#' @param type Currently not implemented. In the future, other types of the test (apart to the one based on db-RDA) should be available.
#' @param cor.coef If \code{type = "cor"}: which correlation coefficient should be used? Partial match to \code{c("pearson", "spearman", "kendal")}.
#' @param exact If \code{type = "cor"} and \code{cor.coef = "spearman"} or \code{cor.coef = 'kendal'} - a logical indicating whether and exact p-value should be computed. 
#' @param alpha Target Type I error rate for Monte Carlo permutation tests (influences number of permutations).
#' @param sqrt Logical value, default FALSE. Should the distance matrix based on Whittaker's index of association be square-rooted to become Euclidean? See Details.
#' @param perm Number of permutations for \code{type = "cor"}. Default = 199.
#' @param x,object Object of the class "testLR"
#' @param digits Number of digits reported for parameters in summary output.
#' @param ... Other arguments passed into \code{print}, \code{summary} or \code{coef} functions (not implemented yet).
#' 
#' @details
#' In case of dbRDA, the matrix of intersample distances is calculated using Whittaker's index of association (\code{\link{ia}}) and significance of the variation explained by sample attributes (R2) is tested by Monte Carlo permutation test. In case of Moran's I, the test is examining wheather the sample attributes variable is compositionally autocorrelated, i.e. whether the Moran's I calculated on this variable using as weighted inverted dissimilarities between sample's species composition is significant.
#' 
#' Whittaker's index of association (calculated as Manhattan type distance on species profiles) is metric, but not Euclidean, and in PCoA (on which dbRDA is based) it can produce negative eigenvalues. After square root transformation, the index becomes both metric and Euclidean.
#' 
#' Variation explained by given environmental variable (R2) can differ between individual weighted means calculated from the same species composition matrix. This happens when different species attributes have different values missing; explained variation is calculated only from those columns (species) of compositional matrix L, which have assigned value for given species attribute. 
#'
#' @return
#' Object of the class "testLR", list of lists, each node containing two parts: results of dbRDA analysis (calculated by \code{\link{capscale}} function from \code{vegan}) and results of Monte Carlo permutation test (calculated by \code{\link{anova.cca}} function, also from \code{vegan}). 
#' @name test.LR
#' @examples
#' data (vltava)
#' test.LR (M = wm (vltava$spe, vltava$ell), vltava$env[,'pH', drop = FALSE], type = 'cor')

#' @export
test.LR <- function (M, env, type = "cor", cor.coef = c("pearson", "kendall", "spearman"), exact = FALSE, alpha = 0.001, sqrt = F, perm = 499)
{
  cor.coef <- match.arg (cor.coef)
  if (!is.wm (M)) stop ("Argument M must be an object of class 'wm'!")
  env <- as.matrix (env)
  if (is.null(colnames (env))) colnames (env) <- 'env'
  res <- list ()
  for (co.M in colnames (M))
    for (co.env in colnames (env))
      res [[co.M]][[co.env]] <- test.LR.0 (M = M[,co.M], env = env[,co.env], type = type, cor.coef = cor.coef, alpha = alpha, sqrt = sqrt, perm = perm)
  class (res) <- 'testLR'
  return (res)
}
#' @rdname test.LR
#' @export
test.LR.0 <- function (M, env, type = 'cor', cor.coef = c('pearson'), exact = FALSE, alpha = 0.001, sqrt = F, perm = 199)
{
  if (type == 'dbRDA')
  {
    sitspe <- attr (M, 'sitspe')[, !is.na (attr (M, 'speatt'))]
    if (sqrt) pcoa.temp <- vegan::capscale (sqrt (ia (sitspe)) ~ env) else pcoa.temp <- vegan::capscale (ia (sitspe) ~ env)
    anova.temp <- anova (pcoa.temp, alpha = alpha)
    res <- list (type = 'dbRDA', pcoa = pcoa.temp, anova = anova.temp)
  }
  if (type == 'moran')
  {
    sitspe <- attr (M, 'sitspe')[, !is.na (attr (M, 'speatt'))]
    residuals <- resid (lm (M ~ as.matrix (env)))
    ia.dist.inv <- as.matrix (1-ia (sitspe))
    diag (ia.dist.inv) <- 0
    moran <- ape::Moran.I (residuals, weight = ia.dist.inv)
    res <- list (type = 'moran', moran = moran)
  }
  if (type == 'corR')
  {
    sitspe <- attr (M, 'sitspe')[, !is.na (attr (M, 'speatt'))]
    env <- as.matrix (env)
    M.art <- wmR (sitspe = sitspe, speatt = wmR (sitspe = t (sitspe), speatt = env))
    cor.res <- cor.test (M.art, env, method = cor.coef)
    obs.stat <- cor.res$statistic
    exp.stat <- replicate (perm, expr = {
      env <- sample (env)
      cor.test (wmR (sitspe = sitspe, speatt = wmR (sitspe = t (sitspe), speatt = env)), env, method = cor.coef)$statistic
    })
    if (cor.coef == 'spearman') P <- sum (abs (c(obs.stat, exp.stat)) <= abs (obs.stat))/(perm+1) else
      P <- sum (abs (c(obs.stat, exp.stat)) >= abs (obs.stat))/(perm+1)  # Spearman's S is sum of variance, increases with decreasing tightness of the correlation
    res <- list (type = 'cor', statistic = obs.stat, estimate = cor.res$estimate, cor = P)
  }
  if (type == 'cor') res <- test_LR_cor (sitspe = as.matrix (attr (M, 'sitspe')), speatt = as.matrix (attr (M, 'speatt')), env = as.matrix (env), cor_coef = cor.coef, perm = perm)
  if (type == 'lm') res <- test_LR_lm (sitspe = as.matrix (attr (M, 'sitspe')), speatt = as.matrix (attr (M, 'speatt')), env = as.matrix (env), perm = perm)
  return (res)
  }

#' @rdname test.LR
#' @export
print.testLR <- function (x, digits = 3, ...)
{
  symnum.pval <- function (pval) symnum( pval, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
  names.speatt <- names (x)
  names.env <- names (x[[1]])
  res <- lapply (x, FUN = function (sp) lapply (sp, FUN = function (en) 
  {
    if (en$type == 'lm') res.temp <- c (format (en$estimate, digits = digits), format (en$lm, digits= digits), format (symnum.pval (en$lm)))
    if (en$type == 'cor') res.temp <- c(format (en$estimate, digits = digits), format (en$cor, digits = digits), format (symnum.pval (en$cor)))
    if (en$type == 'dbRDA') res.temp <- c(format (vegan::RsquareAdj (en$pcoa)$r.squared, digits = digits), format (en$anova[,'Pr(>F)'][1], digits = digits), symnum.pval (en$anova[,'Pr(>F)'][1]))
    if (en$type == 'moran') res.temp <- c(format (en$moran$observed, digits = digits), format (en$moran$p.value, digits = 3), symnum.pval (en$moran$p.value))
    return (res.temp)
  }))
  
  res.m <- matrix (unlist (res), ncol = 3*length (names.env), nrow = length (names.speatt), dimnames = list (names.speatt, as.vector (rbind (names.env, "P", ""))), byrow = T)
  print.default (res.m, quote = F, right = T)
}

#' @rdname test.LR
#' @export
summary.testLR <- function (object, ...)
  print.default (object)

#' @rdname test.LR
#' @export
coef.testLR <- function (object, ...)
{
  names.speatt <- names (object)
  names.env <- names (object[[1]])
  res <- lapply (object, FUN = function (sp) lapply (sp, FUN = function (en) 
  {
    if (en$type == 'cor') res.temp <- c(en$estimate, en$cor)
    if (en$type == 'dbRDA') res.temp <- c(vegan::RsquareAdj (en$pcoa)$r.squared, en$anova[,'Pr(>F)'][1])
    if (en$type == 'moran') res.temp <- c(en$moran$observed, en$moran$p.value)
    return (res.temp)
  }))
  res.m <- matrix (unlist (res), ncol = 2*length (names.env), nrow = length (names.speatt), dimnames = list (names.speatt, as.vector (rbind (names.env, "P values"))), byrow = T)
  return (res.m)
}

