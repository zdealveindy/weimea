% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_wreg.R
\name{test_wreg}
\alias{test_wreg}
\alias{print.testWREG}
\alias{coef.testWREG}
\alias{rescale_axis}
\alias{plot.testWREG}
\title{Testing the relationship species attributes and sample attributes by fourth corner method}
\usage{
test_wreg(
  env,
  com = NULL,
  traits = NULL,
  test = "all",
  parallel = NULL,
  p.adjust.method = p.adjust.methods,
  adjustP = FALSE
)

\method{print}{testWREG}(
  x,
  digits = max(3, getOption("digits") - 3),
  missing.summary = FALSE,
  eps.Pvalue = 0.001,
  signif.stars = getOption("show.signif.stars"),
  ...
)

\method{coef}{testWREG}(object, ...)

rescale_axis(x, var_re, var_st)

\method{plot}{testWREG}(
  x,
  alpha = 0.05,
  line = NA,
  cex.lab = 1.5,
  par.mar = c(0.5, 0.5, 0.5, 0.5),
  box.col = c("blue", "red"),
  box.lwd = 2,
  ...
)
}
\arguments{
\item{env}{Vector, matrix or data frame with sample attributes (environemntal variables).}

\item{com}{Matrix of data frame with species composition data.}

\item{traits}{Vector, matrix or data frame with species attributes (traits)}

\item{test}{NEEDS FIX Vector of character values. Which test should be conducted? Partial match to \code{'rowbased'} row-based permutation test (called Model 2 in the Legendre et al. 1997 and Dray & Legendre 2008), \code{'colbased'} for column-based permutation test (Model 4), and \code{'max'} for max test (selecting the higher from \code{rowbased} and \code{colbased} result, Model 6). \code{'all'} includes all three tests. See \code{Details}.}

\item{parallel}{NULL (default) or integer number. Number of cores for parallel calculation. Maximum number of cores should correspond to number of available cores on the processor.}

\item{p.adjust.method}{A string indicating the method of P-value adjustement, see \code{\link{p.adjust.methods}} for possible choices.}

\item{adjustP}{Logical, default FALSE. Should be the P-values adjusted? If \code{adjustP = TRUE}, the last column in the results is adjusted by method selected in \code{p.adjust.method}.}

\item{x, object}{NEEDS FIX object of the class \code{"cwm"}, generated by function \code{cwm}.}

\item{digits}{number of digits reported by \code{print} method on object of \code{"cwm"} class (default is 3).}

\item{missing.summary}{Logical; should be the summary of values missing in \code{env}, \code{cwm} and \code{traits} be printed along to the result of \code{test_cwm} analysis? Default is \code{TRUE}.}

\item{eps.Pvalue}{Values of P below this threshold will be printed as \code{< [eps]} in the output.}

\item{signif.stars}{Logical; if TRUE, P-values are additionally encoded visually as 'significance stars' in order to help scanning of long coefficient tables. It defaults to the show.signif.stars slot of \code{\link{options}}.}

\item{...}{Other arguments for \code{print}, \code{summary} or \code{coef} functions (not implemented yet).}

\item{var_re, var_st}{Parameters for scaling function \code{rescale_axis}}

\item{alpha, line, cex.lab, par.mar, box.col, box.lwd}{Graphical parameters for \code{plot} function.}
}
\value{
Function \code{test_wreg} returns the list of the class \code{"testWREG"} (with \code{print} and \code{summary} methods), which contains the following items: NEEDS TO FIX
\itemize{
 \item \code{call} Function call.
 \item \code{out} Data frame with results of calculation; the first column contains the original fourth corner metric (\code{r_fc}); following columns contain some of the following results of parametric tests: \code{P.row}, \code{P.col} and \code{P.max}. If \code{adjustP = TRUE}, the last column contains P values adjusted by \code{p.adjust} function.
 \item \code{miss} List with missing values of \code{env}, \code{traits} or \code{cwm} (the last provided for compatibility with \code{test_cwm} function and refers to samples which, after removing missing \code{traits} values, do not contain any species.)
 \item \code{param} Setting of input parameters (\code{perm} and \code{test}).
 }
}
\description{
Function calculating the weightet linear regressions between CWM of species attributes (traits) and sample attributes (environmental variables) and between SNC of sample attributes and species attributes (sensu ter Braak et al. 2018), and testing them by parametric tests and their combination into max test.
}
\details{
In the current implementation, species with missing species attributes are removed from sample x species matrix prior to permutation of species attributes among species.
}
\examples{
data (vltava)
vltava.env <- vltava$env[,c('pH', 'COVERE32')]
vltava.com <- vltava$herbs$spe
vltava.traits <- vltava$herbs$traits
wreg <- test_wreg (env = vltava.env, com = vltava.com, traits = vltava.traits)
}
\seealso{
\code{\link{cwm}} \code{\link{test_cwm}}
}
