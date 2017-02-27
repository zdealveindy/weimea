#' Simulation model sensu Dray & Legendre (2008), creating \bold{R}, \bold{L} and \bold{Q} matrices
#' 
#' Dray & Legendre (2008) suggested an algorithm for generating three matrices (\bold{R}, \bold{L} and \bold{Q}) to test the performance of different permutation scenarios in the fourht corner analysis. This function reproduces their algorithm according to the description given in their article, and implements some additional functionality, such as adding proportional noise to species attributes and/or sample attributes matrix. 
#' 
#' @param n Integer. Number of samples.
#' @param p Integer. Number of species.
#' @param mi.tol Mean species tolerance.
#' @param scenario Simulation scenario. Default 1. See details.
#' @param add.noise Logical value. Should noise be added to all matrices? (defaults to FALSE). Refers to scenario 1N in Dray & Legendre (2008)
#' @param pa Logical value, default FALSE. Should the species composition matrix (\bold{L}) be converted into presence-absence values?
#' @param prop.noise.speatt Real value between 0 and 1. Proportion of noise added to matrix of species attributes by permuting given proportion of rows in matrix \bold{R} randomly among each other.
#' @param prop.noise.env Real value between 0 and 1. Proportion of noise added to matrix of sample attributes by permuting given proportion of rows in matrix \bold{Q} randomly among each other.
#' @param set.seed Integer number to set seed (for reproducibility of results).
#' 
#' @details
#' The artificial data are constructed in the following way (simplified): matrix \bold{R} contains one randomly generated variable with uniform distribution, representing sample attributes, while matrix \bold{Q} contains one randomly generated variable with uniform distribution representing species attributes. To make species and sample attributes linked via species composition, matrix \bold{L} is constructed from both sample and species attributes in a way that response of individual species abundances to the gradient is modelled as symmetric Gaussian curve with optima equivalent to species attributes generated in matrix \bold{Q} and species tolerance generated as random value with mean \code{mi.tol}. Abundance of given species in particular sample is then derived as the probability of species occurrence at particular value of environmental gradient (given by sample attribute in \bold{R}) based on constructed species response curve. For technical details, see the original description in Dray & Legendre (2008), pages 3405-3406.
#' 
#' Six scenarios were suggested by Dray & Legendre (2008) (1, 1N, 2, 3, 4, 5), from which only the scenarios 1, 2, 3 and 4 are currently implemented in this function. Additionally, scenario 0 is added.
#' The definition of scenarios is as follows:
#' \describe{
#' \item{\code{scenario = 1}}{All three matrices (\bold{R}, \bold{L} and \bold{Q}) are linked together by the mechanism of simulation model described above.}
#' \item{\code{scenario = 1N}}{Alternative to \code{scenario = 1} with added random noise. All three matrices are linked together; normal random noise is added to tables \bold{R} and \bold{Q} (mean = 5, sd = 1), and also to matrix \bold{L} (mean = 0, sd = 2).}
#' \item{\code{scenario = 2}}{Species composition (\bold{L}) is linked to sample attributes (\bold{R}), but not to species attributes (\bold{Q}). Matrices are created as in scenario 1, and afterwards the rows with species values in matrix \bold{Q} are permuted (cancelling the link between \bold{L} and \bold{Q}).}
#' \item{\code{scenario = 3}}{Species composition (\bold{L}) is linked to species attributes (\bold{Q}), but not to sample attributes (\bold{R}). Matrices are created as in scenario 1, and afterwards rows in matrix \bold{R} are permuted (this cancels link between \bold{L} and \bold{R}).}
#' \item{\code{scenario = 4}}{There is no link between \bold{L} and \bold{Q}, neither between \bold{L} and \bold{R}. Matrices are created as in scenario 1, and afterwards the rows in both matrices \bold{R} and \bold{Q} are permuted, cancelling all links between matrices.}
#' \item{\code{scenario = 0}}{This scenario represents the continuous transition between all above scenarios (1 to 4). Modifiying the proportion of noise in species attributes (argument \code{prop.noise.specatt}) and in sample attributes (argument \code{prop.noise.sampatt}) enables to create intermediary scenarios with varying degree to which the matrices are connected to each other. The following settings of arguments is analogous to the scenarios mentioned above:
#' \itemize{
#' \item \code{prop.noise.speatt = 0} and \code{prop.noise.env = 0} is analogous to \code{scenario = 1}
#' \item \code{prop.noise.speatt = 1} and \code{prop.noise.env = 0} is analogous to \code{scenario = 2}
#' \item \code{prop.noise.speatt = 0} and \code{prop.noise.env = 1} is analogous to \code{scenario = 3}
#' \item \code{prop.noise.speatt = 1} and \code{prop.noise.env = 1} is analogous to \code{scenario = 4}
#' }}}
#' 
#' 
#' @return
#' The function returns the list of three items:
#' \tabular{lll}{
#' \tab \code{sitspe} \tab Matrix with species composition data.\cr
#' \tab \code{env}    \tab Matrix with sample attributes (environmnetal data).\cr
#' \tab \code{speatt} \tab Matrix with species attributes.
#' }
#' @author David Zeleny (zeleny.david@@gmail.com), created mostly according to the description in Dray & Legendre (2008).
#' 
#' @references
#' Dray S. & Legendre P. 2008. Testing the species traits-environment relationships: the fourth-corner problem revisited. Ecology, 89:3400-3412.
#' @export
#' 

# the function to simulate the RLQ matrices (sensu Dray & Legendre 2008)
simul.RLQ <- function (n = 100, p = 100, mi.tol = 10, scenario = 1, add.noise = F, pa = F, prop.noise.speatt = 0, prop.noise.env = 0, set.seed = NULL)
{
  if (!is.null (set.seed)) set.seed (set.seed)
  noise.var <- function (var, prop = 0.5)
  {
    sel <- sample (1:length (var), round (length (var)*prop))
    var.r <- var
    var.r[sel] <- sample (var[sel])
    var.r
  }


  # step 1
  x <- runif (n, 0, 100)  # position along environmental gradient
  
  # step 2
  mi <- runif (p, -1, 101)  # species optima
  
  # step 3
  h <- runif (p, 0.5, 1)  # height of species j in it's optimum
  
  # step 4
  sigma <- rnorm (p, mean = mi.tol, sd = 10)
  
  # step 5
  y <- matrix (nrow = n, ncol = p)
  
  for (j in seq (1, p))
    y[,j] <- h[j]*exp (-(x - mi[j])^2/(2*((sigma[j])^2)))
  
  if (pa) y <- apply (y, 2, FUN = function (x) rbinom (length (x), 1, prob = x)) 
  
  if (scenario == 2 | scenario == 4) mi <- sample (mi)
  
  if (scenario == 3 | scenario == 4) x <- sample (x)
  
  if (scenario == 0) 
  {
    mi <- noise.var (mi, prop.noise.speatt)
    x <- noise.var (x, prop.noise.env)
  }

  if (add.noise)
  {
    y.n <- apply (y, 2, FUN = function (x) x + rnorm (n = length (x), mean = 0, sd = 2))
    y.n[y.n < 0] <- 0
    x.n <- x + rnorm (length (x), mean = 5, sd = 1)
    mi.n <- mi + rnorm (length (mi), mean = 5, sd = 1)
  }
    
  if (add.noise) res <- list (sitspe = as.data.frame (y.n), env = as.data.frame (x.n), speatt = as.data.frame (mi.n)) else res <- list (sitspe = as.data.frame (y), env = as.data.frame (x), speatt = as.data.frame (mi))
  return (res)
}