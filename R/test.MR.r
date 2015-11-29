#' Testing the relationship between weighted-mean of species attributes with sample attributes
#' 
#' Function performing standard, modified, two-step and sequential test to calculate significance of relationship between weighted-mean of species attributes and sample attributes.
#' 
#' @name test.MR
#' @param M An object of the class \code{wm} 
#' @param env Vector or matrix with variables. See details.
#' @param method Statistical method used to analyse the relationship between M (of class \code{wm}) and env (sample attributes); partial match to \code{'lm'}, \code{'aov'}, \code{'cor'}, \code{'kruskal'}, \code{'slope'} and \code{'fourthcorner'}. Currently only \code{'cor'} method implemented.
#' @param cor.coef Correlation coefficient for \code{method = 'cor'}. Partial match to \code{'pearson'}, \code{'spearman'} and \code{'kendal'}. Currently only \code{'pearson'} implemented.
#' @param dependence Should \code{M} be dependent variable and \code{env} independent (\code{'M ~ env'}), or opposite? Applicable only for \code{method = 'lm'}. Partial match to \code{'M ~ env'} and \code{'env ~ M'}.
#' @param perm Number of permutations.
#' @param test Vector of character values. Which test should be conducted? Partial match of vector with \code{'standard'}, \code{'modified'}, \code{'twostep'}, \code{'sequential'} and \code{'all'}.
#' @param fc.test Test for the fourthcorner method. \code{fc.test = c(2,4,6)}.
#' @param testLR.perm Relevant if \code{test = "twostep"} (default): number of permutations for test of realtionship between \bold{L} and \bold{R} matrices.
#' @param testLR.P Relevant if \code{test = "twostep"} (default): significance value at whitch the relationship between \bold{L} and \bold{R} matrices are significant.
#' @param parallel NULL (default) or integer number. Number of cores for parallel calculation of modified permutation test. Maximum number of cores should correspond to number of available cores on the processor.
#' @param x,object object of the class \code{"wm"}, generated by function \code{wm}.
#' @param digits number of digits reported by \code{print} method on object of \code{"wm"} class (default is 3).
#' @param ... Other arguments for \code{print}, \code{summary} or \code{coef} functions (not implemented yet).

#' @export
#' @examples 
#' data (vltava)
#' M <- wm (sitspe = vltava$herbs$spe, speatt = vltava$herbs$traits)
#' re <- test.MR (M = M, env = vltava$env[,c('pH', 'COVERE32')])
#' @details
#' Currently implemented statistical method is correlation (\code{'cor'}).
#'
#' Argument \code{env} can be vector or matrix with one column. Only in case of linear regression (\code{method = 'lm'}) is possible to use matrix with several variables, which will be all used as independent variables in the model. For ANOVA and Kruskal-Wallis test, make sure that 'env' is \code{factor} (warning will be returned if this is not the case, but the calculation will be conducted). 
#' 
#' Difference between \code{method = 'lm'} and \code{'aov'} is in the format of summary tables, returned by \code{summary.wm} function. In case of 'aov', this summary is expressed in the traditional language of ANOVA rather than linear models.
#' 
#' Both \code{method = 'lm'} and \code{'slope'} are based on linear regression and calculated by function \code{\link{lm}}, but differ by test statistic: while 'lm' is using F value and is testing the strength of the regression (measured by r2), 'slope' is using the slope of the regression line (b). This statistic is added here for comparison with the fourth corner method. While r2 (and r) is influenced by the issue of compositional autocorrelation, slope of regression is not.
#' 
#' Specific issue related to weighted mean is the case of missing species attributes. In current implementation, species with missing species attributes are removed from sample x species matrix prior to permutation of species attributes among species. 
#' @return  Function \code{wm} returns list of the class \code{"wm"} (with \code{print} and \code{summary} methods), which contains the following items:
#' \itemize{
#'  \item \code{real.summaries} summary of the method
#'  \item \code{coefs} model coefficients
#'  \item \code{stat} test statistic
#'  \item \code{orig.P} P-values from the original (parametric) test
#'  \item \code{perm.P} P-values from the not-modified permutation test (permuted are the whole rows of matrix M)
#'  \item \code{modif.P} P-values from the modified permutation test (permuted are species attributes in object M)
#'  \item \code{seq.P} P-values from sequential test (\code{max (perm.P, modif.P)})
#'  \item \code{perm} number of permutations
#'  \item \code{method} method of calculation
#'  \item \code{dependence} dependence (in case of \code{method = 'lm'})
#'  }
#' @seealso \code{\link{wm}}

#' @export
test.MR <- function (M, env, method = c('cor'), cor.coef = c('pearson'), dependence = "M ~ env", perm = 499, testLR.perm = NULL, test = "twostep", fc.test = 6, parallel = NULL, testLR.P = 0.05)
{
  METHOD <- c('lm', 'aov', 'cor', 'kruskal', 'slope', 'fourthcorner')
  COR.COEF <- c('pearson', 'spearman', 'kendall')
  TEST <- c('none', 'standard', 'modified', 'twostep', 'sequential', 'all')
  DEPENDENCE <- c("M ~ env", "env ~ M")
  FC.TEST <- c(2,4,6)
  method <- match.arg (method, METHOD)
  cor.coef <- match.arg (cor.coef, COR.COEF)
  test <- match.arg (test, TEST, several.ok = T)
  dependence <- match.arg (dependence, DEPENDENCE)
  if (is.null (testLR.perm)) testLR.perm <- perm
  if (!is.wm (M) & ("modified" %in% test || "twostep" %in% test || "sequential" %in% test || "all" %in% test)) stop ("Object M must be of 'wm' class")
  env <- as.data.frame (as.matrix (env))
  if (method == 'lm' & test == 'twostep' & dependence == 'M ~ env' & ncol (env) > 1) stop ("Two-step permutation method is not available for multiple linear regression with more than one predictor variable (env)")
  sitspe <- attr (M, 'sitspe')
  speatt <- attr (M, 'speatt')

# correlation - function 'cor'
  if (method == 'cor')
  {
    vars.names <- expand.grid (env = colnames (env), speatt = colnames (speatt), stringsAsFactors = FALSE)[,2:1]

    if (is.null (parallel)) 
      res <- apply (vars.names, 1, FUN = function (var12) test_MR_cor (sitspe = as.matrix (sitspe), speatt = as.matrix (speatt[, var12['speatt'], drop = F]), env = as.matrix (env[, var12['env'], drop = F]), test = test, cor_coef = cor.coef, perm = perm, testLR_P = testLR.P, testLR_perm = testLR.perm)) else 
      {
        cl <- parallel::makeCluster (parallel)
        parallel::clusterExport (cl, varlist = c('vars.names', 'sitspe', 'speatt', 'env', 'test', 'cor.coef', 'perm', 'testLR.P', 'testLR.perm'), envir = environment ())
        parallel::clusterEvalQ (cl, eval (call ('library', 'weimea')))
        res <- parallel::parApply (cl, vars.names, 1, FUN = function (var12) test_MR_cor (sitspe = as.matrix (sitspe), speatt = as.matrix (speatt[, var12['speatt'], drop = F]), env = as.matrix (env[, var12['env'], drop = F]), test = test, cor_coef = cor.coef, perm = perm, testLR_P = testLR.P, testLR_perm = testLR.perm))
      }
    names (res) <- apply (vars.names, 1, FUN = function (x) paste (x[1], x[2], sep = ' / '))
      }
    
    # linear regression - function fastLM_wm
    if (method == 'lm')
    {
      vars.names <- if (dependence == 'M ~ env') colnames (M) else colnames (env)
      
      if (is.null (parallel)) {
        if (dependence == 'M ~ env') 
          res <- lapply (vars.names, FUN = function (var1) test_MR_lm (sitspe = as.matrix (sitspe), speatt = as.matrix (speatt[, var1, drop = F]), env = as.matrix (env), test = test, dependence = dependence, perm = perm, testLR_P = testLR.P, testLR_perm = testLR.perm)) else
            res <- lapply (vars.names, FUN = function (var1) test_MR_lm (sitspe = as.matrix (sitspe), speatt = as.matrix (speatt), env = as.matrix (env[, var1, drop = F]), test = test, dependence = dependence, perm = perm, testLR_P = testLR.P, testLR_perm = testLR.perm))
      } else
        {
          cl <- parallel::makeCluster (parallel)
          parallel::clusterExport (cl, varlist = c('vars.names', 'sitspe', 'speatt', 'env', 'test', 'dependence', 'perm', 'testLR.P', 'testLR.perm'), envir = environment ())
          parallel::clusterEvalQ (cl, eval (call ('library', 'weimea')))
          if (dependence == 'M ~ env') 
            res <- parallel::parLapply (cl, vars.names, FUN = function (var1) test_MR_lm (sitspe = as.matrix (sitspe), speatt = as.matrix (speatt[, var1, drop = F]), env = as.matrix (env), test = test, dependence = dependence, perm = perm, testLR_P = testLR.P, testLR_perm = testLR.perm)) else
              res <- parallel::parLapply (cl, vars.names, FUN = function (var1) test_MR_lm (sitspe = as.matrix (sitspe), speatt = as.matrix (speatt), env = as.matrix (env[, var1, drop = F]), test = test, dependence = dependence, perm = perm, testLR_P = testLR.P, testLR_perm = testLR.perm))
        }
    lapply (res, FUN = function (x) x$coef)  
    names (res) <- vars.names
    for (i in seq (1, length (res))) rownames (res[[i]]$coef) <- c ('intercept', names (env))
    }

  if (method == 'fourthcorner')   # implementation of fourt.corner.ade doesn't allow to go parallel
  {
      res <- fourth.corner.ade (sitspe = sitspe, speatt = speatt, env = env, fc.test = fc.test, perm = perm)
  }
  result <- list (out = res, param = list (method = method, cor.coef = cor.coef, dependence = dependence, perm = perm, test = test))
  if ('twostep' %in%  test) result [c('testLR.perm', 'testLR.P')] <- c (testLR.perm, testLR.P)
  class (result) <- 'testMR'
  return (result)
}

#' @rdname test.MR
#' @export
print.testMR <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
  symnum.pval <- function (pval) symnum( pval, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
  if (x$param$method == 'cor')  out.temp0 <- do.call (rbind.data.frame, x$out)
  if (x$param$method == 'lm') out.temp0 <- do.call (rbind.data.frame, lapply (x$out, FUN = function (y) as.data.frame (t(do.call (rbind.data.frame, y)))))
  if (x$param$method == 'fourthcorner') out.temp0 <- x$out
  out.temp0 <- out.temp0[,!is.na (colSums (out.temp0)), drop = F]
  P.columns <- which (substr (colnames (out.temp0), start = 1, stop = 2) == 'P.')
  out.temp <- out.temp0[,1:(P.columns[1]-1), drop = F]
  for (co in P.columns)
    {
    out.temp <- cbind (out.temp, out.temp0[,co], symnum.pval (out.temp0[,co]))
    colnames (out.temp)[c(ncol (out.temp)-1, ncol (out.temp))] <- c(colnames (out.temp0)[co], ' ')
  }
  print (out.temp)
}

#' @rdname test.MR
#' @export
coef.testMR <- function (object, ...)
  return (do.call (rbind.data.frame, object$out))
  
