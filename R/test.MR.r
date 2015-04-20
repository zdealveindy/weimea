#' Testing the relationship between weighted-mean of species attributes with sample attributes
#' 
#' Function performing standard, modified, combined and sequential test to calculate significance of relationship between weighted-mean of species attributes and sample attributes.
#' 
#' @name test.MR
#' @param M An object of the class \code{wm} 
#' @param env Vector or matrix with variables. See details.
#' @param method Statistical method used to analyse the relationship between M (of class \code{wm} and env); partial match to \code{'lm'}, \code{'aov'}, \code{'cor'}, \code{'kruskal'}, \code{'slope'} and \code{'fourthcorner'}. Currently only \code{'cor'} method implemented.
#' @param cor.coef Correlation coefficient in case of \code{method = 'cor'}. Partial match to 'pearson', 'spearman' and 'kendal'.
#' @param dependence Should M be dependent variable and env independent (\code{'M ~ env'}), or opposite? Applicable only for \code{method = 'lm'}. Partial match to \code{'M ~ env'} and \code{'env ~ M'}, so to write \code{dep = 'M'} is enough.
#' @param permutations Number of permutations.
#' @param test Vector of character values. Which test should be conducted? Partial match of vector with \code{'standard'}, \code{'modified'}, \code{'combined'}, \code{'sequential'} and \code{'all'}.
#' @param parallel NULL (default) or integer number. Number of cores for parallel calculation of modified permutation test. Maximum number of cores should correspond to number of available cores on the processor.
#' @param object,x object of the class \code{"mopet"}, generated by function \code{mopet}.
#' @param digits number of digits reported by \code{print} method on object of \code{"mopet"} class (default is 3).
#' @param ... Other arguments for \code{print} or \code{summary} functions (not implemented).

#' @export
#' @examples 
#' data (vltava)
#' M <- wm (sitspe = vltava$herbs$spe, speatt = vltava$herbs$traits)
#' re <- test.MR (M = M, env = vltava$env[,c('pH', 'COVERE32')])
#' @details
#' Currently implemented statistical methods are correlation (\code{'cor'}), linear regression (\code{method = 'lm'}), ANOVA (\code{'aov'}) and Kruskal-Wallis test (\code{'kruskal'}).
#'
#' Argument \code{env} can be vector or matrix with one column. Only in case of linear regression (\code{method = 'lm'}) is possible to use matrix with several variables, which will be all used as independent variables in the model. For ANOVA and Kruskal-Wallis test, make sure that 'env' is \code{factor} (warning will be returned if this is not the case, but the calculation will be conducted). 
#' 
#' Difference between \code{method = 'lm'} and \code{'aov'} is in the format of summary tables, returned by \code{summary.mopet} function. In case of 'aov', this summary is expressed in the traditional language of ANOVA rather than linear models.
#' 
#' Both \code{method = 'lm'} and \code{'slope'} are based on linear regression and calculated by function \code{\link{lm}}, but differ by test statistic: while 'lm' is using F value and is testing the strength of the regression (measured by r2), 'slope' is using the slope of the regression line (b). This statistic is added here for comparison with the fourth corner method. While r2 (and r) is influenced by the issue of compositional autocorrelation, slope of regression is not.
#' 
#' Specific issue related to weighted mean is the case of missing species attributes. In current implementation, species with missing species attributes are removed from sample x species matrix prior to permutation of species attributes among species. 
#' @return  Function \code{mopet} returns list of the class \code{"mopet"} (with \code{print} and \code{summary} methods), which contains the following items:
#' \itemize{
#'  \item \code{real.summaries} summary of the method
#'  \item \code{coefs} model coefficients
#'  \item \code{stat} test statistic
#'  \item \code{orig.P} P-values from the original (parametric) test
#'  \item \code{perm.P} P-values from the not-modified permutation test (permuted are the whole rows of matrix M)
#'  \item \code{modif.P} P-values from the modified permutation test (permuted are species attributes in object M)
#'  \item \code{seq.P} P-values from sequential test (\code{max (perm.P, modif.P)})
#'  \item \code{permutations} number of permutations
#'  \item \code{method} method of calculation
#'  \item \code{dependence} dependence (in case of \code{method = 'lm'})
#'  }
#' @seealso \code{\link{wm}}

#' @export
test.MR <- function (M, env, method = c('cor'), cor.coef = c('pearson'), dependence = "M ~ env", permutations = 199, test = "combined", parallel = NULL, P.testLR = 0.05)
{
  METHOD <- c('lm', 'aov', 'cor', 'kruskal', 'slope', 'fourthcorner')
  COR.COEF <- c('pearson', 'spearman', 'kendall')
  TEST <- c('none', 'standard', 'modified', 'combined', 'sequential', 'all')
  DEPENDENCE <- c("M ~ env", "env ~ M")
  method <- match.arg (method, METHOD)
  cor.coef <- match.arg (cor.coef, COR.COEF)
  test <- match.arg (test, TEST, several.ok = T)
  dependence <- match.arg (dependence, DEPENDENCE)
  if (!is.wm (M) & ("modified" %in% test || "combined" %in% test || "sequential" %in% test || "all" %in% test)) stop ("Object M must be of 'wm' class")
  env <- as.data.frame (as.matrix (env))
  sitspe <- attr (M, 'sitspe')
  speatt <- attr (M, 'speatt')

# correlation - function 'cor'
  if (method == 'cor')
  {
    vars.names <- expand.grid (env = colnames (env), M = colnames (M))[,2:1]
    vars.order <- expand.grid (env = 1:ncol (env), M = 1:ncol (M))[,2:1]
    res <- apply (vars.order, 1, FUN = function (var12)
      {
      obs <- cor.test (M[,var12[[1]]], env[,var12[[2]]], method = cor.coef)
      obs.coef <- obs$estimate
      obs.stat <- obs$statistic
      P.param <- obs$p.value
      P.stand <- NULL
      P.modif <- NULL
      P.comb <- NULL
      P.seq <- NULL
      if ('standard' %in% test) {
        exp.stat <- replicate (permutations, cor.test (M[,var12[[1]]], sample (env[,var12[[2]]]), method = cor.coef)$statistic)
        P.stand <- sum (abs (c(obs.stat, exp.stat)) >= abs (obs.stat))/(permutations+1)
      }
      if ('modified' %in% test) {
        exp.stat <- replicate (permutations, cor.test (randomize (M[,var12[[1]]]), env[,var12[[2]]], method = cor.coef)$statistic)
        P.modif <- sum (abs (c(obs.stat, exp.stat)) >= abs (obs.stat))/(permutations+1)
      }
      if ('combined' %in% test) {
        if (test.LR (M[,var12[[1]]], env[,var12[[2]]])[[1]][[1]][[2]] < P.testLR) exp.stat <- replicate (permutations, cor.test (randomize (M[,var12[[1]]]), env[,var12[[2]]], method = cor.coef)$statistic) else exp.stat <- replicate (permutations, cor.test (M[,var12[[1]]], sample (env[,var12[[2]]]), method = cor.coef)$statistic)
        P.comb <- sum (abs (c(obs.stat, exp.stat)) >= abs (obs.stat))/(permutations+1)
      }
      if ('sequential' %in% test) {
        exp.stat.stand <- replicate (permutations, cor.test (M[,var12[[1]]], sample (env[,var12[[2]]]), method = cor.coef)$statistic)
        P.stand.seq <- sum (abs (c(obs.stat, exp.stat.stand)) >= abs (obs.stat))/(permutations+1)
        exp.stat.modif <- replicate (permutations, cor.test (randomize (M[,var12[[1]]]), env[,var12[[2]]], method = cor.coef)$statistic)
        P.modif.seq <- sum (abs (c(obs.stat, exp.stat.modif)) >= abs (obs.stat))/(permutations+1)
        P.seq <- max (P.stand.seq, P.modif.seq)
      }
      list (obs.coef = obs.coef, obs.stat = obs.stat, P.param = P.param, P.stand = P.stand, P.modif = P.modif, P.comb = P.comb, P.seq = P.seq)
    })
  }

names (res) <- apply (vars.names, 1, FUN = function (x) paste (x[1], x[2], sep = ' / '))
res.temp <- matrix (NA, ncol = length (res[[1]]), nrow = length (res))
rownames (res.temp) <- names (res)
colnames (res.temp) <- names (res[[1]])
for (i in seq (1, length (res[[1]])))
{
  repl <- unlist (lapply (res, FUN = function (x) x[[i]]))
  if (!is.null (repl)) res.temp[,i] <- repl  
}
result <- list (summary = res.temp, param = list (method = method, cor.coef = cor.coef, dependence = dependence, permutations = permutations, test = test))
class (result) <- 'testMR'
return (result)
}

#' @export
print.testMR <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
  symnum.pval <- function (pval) symnum( pval, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
  summary.temp <- x$summary[,!is.na (colSums (x$summary))]
  print (summary.temp)
}