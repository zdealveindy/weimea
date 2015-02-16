#' Simulation of ecological data structured by two ecological gradients
#' @author David Zeleny (zeleny.david@@gmail.com) - based on the script of Jason Fridley (Fridley et al. 2007, Appendix S2), modified for two dimensions


#' @param totS total number of species in simulation
#' @param gr1.length,gr2.length length of first and second gradient in units
#' @param niche.type shape of species response curves ('random', 'normal', 'skewed')
#' @param max.niche.breath vector of the length = 2, with maximum niche breath along first and second gradient (default = \code{c(gr1.length, gr2.length)})
#' @param min.niche.breath vector of the length = 2, with minimum niche breath along the first and second gradient (default = \code{c(10, 10)}). Cannot be more than \code{max.niche.breath}.
#' @param seed set random seed for reproducing always the same result
#' @param plotting should the species response curves along simulated gradient be drawned? Default = \code{FALSE}.
#' @param highlight.species vector of species numbers which should be highlighted by color in ploted diagram (if \code{plotting = TRUE}). Default = \code{NULL}, which means that if plot is drawned, two species will be highlighted - the most generalized and the most specialized ones.
#' @param simul.comm result of simul.comm.2 function with parameters of individual species response curves
#' @param Np number of samples
#' @param sample.x1,sample.x2 positions of sampling along the first (or second, respectively) gradient (default = NULL, means that random locations are generated)
#' @param no.ind mean number of individuals to be drawn from the species pool (if \code{based.on = 'individuals'})
#' @param k mean proportion of species from species pool to be drawn into local community (if \code{based.on = 'species'})
#' @param pa result table will be generated in presence-absence form (default \code{pa = TRUE})
#' @param based.on should the sampling of species from species pool be based on \code{'individuals'} or \code{'species'}? (default = \code{'individuals'})
#' @examples
#' sc <- sample.comm.2 (simul.comm.2 (gr1.length = 5000, gr2.length = 3000, totS = 300),
#'  Np = 200, pa = FALSE, seed = 123)
#' library (vegan)
#' nmds <- metaMDS (sc$a.mat)
#' par (mfrow = c(1,2))
#' ordiplot (nmds, display = 'si', main = 'Environmental gradient 1')
#' points (nmds, display = 'si', cex = scale (sc$sample.x1, center = FALSE)*2, pch = 21, bg = 'white')
#' ordiplot (nmds, display = 'si', main = 'Environmental gradient 2')
#' points (nmds, display = 'si', cex = scale (sc$sample.x2, center = FALSE)*2, pch = 21, bg = 'white')
#' @return The function \code{simul.comm.2} returns \code{list} with 14 items, describing the set of parameters used to simulate community of species response curves:
#' \itemize{
#' \item \code{totS} Total number of species in simulation.
#' \item \code{gr1.length, gr2.length} Length of the first and second gradient, respectively.
#' \item \code{niche.type} Shape of species response curves.
#' \item \code{Ao1, Ao2} Vector of species amplitudes for the first and second gradient, respectively (heights of species response curves, corresponding to maximum probability of species to be selected to community in species optima).
#' \item \code{m1, m2} Vector of species optima along the first and second gradient, respectively.
#' \item \code{r1, r2} Vector of generated species ranges along the first and second gradient, respectively (generated niche breaths).
#' \item \code{range1, range2} Vector of realised species ranges along the first and second gradient, respectively, considering the truncation of species response curves by gradient margins (these differ from \code{r1} and \code{r2} especially at the gradient margins, where the generated species niche may be wide, but since the margin cuts the species occurrences, realised species niche is narrower.)
#' \item \code{a1, a2, g1, g2} Vectors of shape parameters for curves (used in beta function to generate the shape of the species response curve)
#' \item \code{A1.all, A2.all} Matrix (dim = gradient length x number of species) with simulated probabilites of individual specie at individual location along the gradient. 
#' }
#' The function \code{sample.comm.2} returns \code{list} of 8 items with parameters of generated community data; the last item contains also all items returned by \code{simul.comm.2} function:
#' \itemize{
#' \item \code{a.mat} Matrix (sample x species) of species abundances in samples.
#' \item \code{p.mat} Matrix (sample x species) of species occurrence probabilities in samples.
#' \item \code{sample.x1, sample.x2} Vector with positions of samples along the first and second simulated gradient, respectively (environmental variable).
#' \item \code{sample.comm} List of 6 items storing initial settings of arguments in \code{sample.comm.2} (namely arguments \code{Np, based.on, no.ind, k, seed} and \code{pa}).
#' \item \code{simul.comm} List of 14 items returned by function \code{simul.comm.2}.
#' }




#' @rdname simul.comm.2
#' @export
simul.comm.2 <- function (totS = 300, gr1.length = 5000, gr2.length = 5000, niche.type = 'random', max.niche.breath = c(gr1.length, gr2.length), min.niche.breath = c(10, 10), prop.random.species = c(0,0), seed = NULL, plotting = F, highlight.species = NULL)
{
  if (!is.null (seed)) set.seed (seed)
  
  #This is beta function for generating niches
  curve <- function(Ao,m,r,a,g,x)  {
    (Ao*((((x-m)/r)+(a/(a+g)))^a)*((1-(((x-m)/r)+(a/(a+g))))^g))/(((a/(a+g))^a)*((1-(a/(a+g)))^g))
  }
  
  x1 <- seq(1, gr1.length, by = 1) # gradient 1
  x2 <- seq(1, gr2.length, by = 1) # gradient 2
  S1 <- totS #number of species
  S2 <- totS #number of species
  # species values for random niches
  if(niche.type=="random") {
    Ao1<-rlnorm(S1,2,1)  		#amplitude vector (lognormal distribution)
    m1<-sample(seq(5:max(x1)-5),S1, replace = T)		#location of optima
    r1<-runif(S1,min=min.niche.breath[1], max=max.niche.breath[1])		#range along gradient (niche breadth)
    a1<-(runif(S1,min=.1,max=4))		#shape parameter (alpha)
    g1<-(runif(S1,min=.1,max=4))		#shape parameter (gamma)
    
    Ao2<-rlnorm(S2,2,1)			#amplitude vector (lognormal distribution)
    m2<-sample(seq(5:max(x2)-5),S2, replace = T)		#location of optima
    r2<-runif(S2,min=min.niche.breath[2], max=max.niche.breath[2])		#range along gradient (niche breadth)
    a2<-(runif(S2,min=.1,max=4))		#shape parameter (alpha)
    g2<-(runif(S2,min=.1,max=4))		#shape parameter (gamma)
  }
  
  # species values for normal niches
  if(niche.type=="normal") {
    Ao1<-rlnorm(S1,2,1)			#amplitude vector (lognormal distribution)
    m1<-sample(seq(5:max(x1)-5),S1, replace = T)		#location of optima
    r1<-runif(S1,min=min.niche.breath[1], max=max.niche.breath[1])		#range along gradient (niche breadth)
    a1<-rep(1.99,S1)				#shape parameter (alpha)
    g1<-rep(1.99,S1)				#shape parameter (gamma)
    
    Ao2<-rlnorm(S2,2,1)			#amplitude vector (lognormal distribution)
    m2<-sample(seq(5:max(x2)-5),S2, replace = T)		#location of optima
    r2<-runif(S2,min=min.niche.breath[2], max=max.niche.breath[2])		#range along gradient (niche breadth)
    a2<-rep(1.99,S2)				#shape parameter (alpha)
    g2<-rep(1.99,S2)				#shape parameter (gamma)
  }
  
  # species values for skewed niches
  if(niche.type=="skewed") {
    Ao1<-rlnorm(S1,2,1)			#amplitude vector (lognormal distribution)
    m1<-sample(seq(5:max(x1)-5),S1, replace = T)		#location of optima
    r1<-runif(S1,min=min.niche.breath[1],max=max.niche.breath[2])		#range along gradient (niche breadth)
    #produce skews in either direction, randomly
    a.1<-rep(1.99,S1)
    g.1<-rep(.25,S1)
    samp<-sample(c(1:S1),S1/2,replace=FALSE)
    a.1[samp]<-.25
    g.1[samp]<-1.99
    a1<-a.1				#shape parameter (alpha)
    g1<-g.1				#shape parameter (gamma)
    
    Ao2<-rlnorm(S2,2,1)			#amplitude vector (lognormal distribution)
    m2<-sample(seq(5:max(x2)-5),S2, replace = T)		#location of optima
    r2<-runif(S2,min=min.niche.breath[2],max=max.niche.breath[2])		#range along gradient (niche breadth)
    #produce skews in either direction, randomly
    a.1<-rep(1.99,S2)
    g.1<-rep(.25,S2)
    samp<-sample(c(1:S2),S2/2,replace=FALSE)
    a.1[samp]<-.25
    g.1[samp]<-1.99
    a2<-a.1				#shape parameter (alpha)
    g2<-g.1				#shape parameter (gamma)
  }
  
  
  A1.all <- matrix(0, nrow = gr1.length, ncol = S1) #response abundances for all points along gradient 1
  for(L in 1:S1){
    A1.all[,L] <- curve(Ao1[L],m1[L],r1[L],a1[L],g1[L], x1)
  }
  
  A2.all <- matrix(0, nrow = gr2.length, ncol = S2) #response abundances for all points along gradient 2
  for(L in 1:S2){
    A2.all [,L] <- curve(Ao2[L],m2[L],r2[L],a2[L],g2[L], x2)
  }
  A1.all[is.na (A1.all)] <- 0
  A2.all[is.na (A2.all)] <- 0

  if (prop.random.species[1] > 0)
    for (co in sample (1:ncol (A1.all), round (prop.random.species[1]*ncol (A1.all))))
      A1.all[,co] <- sample (A1.all[,co])
  if (prop.random.species[2] > 0) A2.all
  for (co in sample (1:ncol (A2.all), round (prop.random.species[2]*ncol (A2.all))))
    A2.all[,co] <- sample (A2.all[,co])
    
  range1 <- colSums (A1.all > 0)
  range2 <- colSums (A2.all > 0)
    
  if (plotting)
  {
    # Plot species distributions along gradient
    par(mfrow=c(2,1))
    
    plot(x1,A1.all[,1],xlim=c(0,max(x1)),ylim=c(0,max(Ao1)),type="l",xlab="Gradient 1",ylab="Abundance",cex.lab=1.7,cex.axis=1.5,lwd=2)
    for(L in 2:S1) {
      lines(x1,A1.all[,L])
    }
    if (is.null (highlight.species))
    {
      lines(x1,A1.all[,which (range1 == max (range1))[1]],lwd=3,col=2)
      lines(x1,A1.all[,which (range1 == min (range1))[1]],lwd=3,col=3)
    } else 
      for (co in seq (1, length (highlight.species)))
        lines(x1,A1.all[,highlight.species[co]],lwd=3,col=co+1)
    
    plot(x2,A2.all[,1],xlim=c(0,max(x2)),ylim=c(0,max(Ao2)),type="l",xlab="Gradient 2",ylab="Abundance",cex.lab=1.7,cex.axis=1.5,lwd=2)
    for(L in 2:S2) {
      lines(x2,A2.all[,L])
    }
    if (is.null (highlight.species))
    {
      lines(x2,A2.all[,which (range2 == max (range2))[1]],lwd=3,col=2)
      lines(x2,A2.all[,which (range2 == min (range2))[1]],lwd=3,col=3)  
    } else
      for (co in seq (1, length (highlight.species)))
        lines(x2,A2.all[,highlight.species[co]],lwd=3,col=co+1)
  }  

  
  result <- list (totS = totS, gr1.length = gr1.length, gr2.length = gr2.length, niche.type = niche.type, Ao1 = Ao1, Ao2 = Ao2, m1 = m1, m2 = m2, r1 = r1, r2 = r2, range1 = range1, range2 = range2, a1 = a1, a2 = a2, g1 = g1, g2 = g2, A1.all = A1.all, A2.all = A2.all)
  result
}

#' @rdname simul.comm.2
#' @export
sample.comm.2 <- function (simul.comm = simul.comm.2 (), Np = 300, sample.x1 = NULL, sample.x2 = NULL, no.ind = 100, k = 0.2, seed = NULL, pa = F, based.on = 'individuals')
{
  if (!is.null (seed)) set.seed (seed)
  sc <- simul.comm
  BASED.ON <- c('individuals', 'species')
  based.on <- match.arg (based.on, BASED.ON)
  
  #  #This is beta function for generating niches
  #  curve <- function(Ao,m,r,a,g,x)  {
  #    (Ao*((((x-m)/r)+(a/(a+g)))^a)*((1-(((x-m)/r)+(a/(a+g))))^g))/(((a/(a+g))^a)*((1-(a/(a+g)))^g))
  #  }
  
  #Random sample intervals along gradient
  if (is.null (sample.x1)) sample.x1 <- trunc(sample(c(2:sc$gr1.length)-1,Np, replace = T))
  if (is.null (sample.x2)) sample.x2 <- trunc(sample(c(2:sc$gr2.length)-1,Np, replace = T))
  
  #  A1.all <- matrix(0,nrow=simul.comm$gr1.length, ncol=sc$totS) #response abundances
  #  for(L in 1:sc$totS){
  #    A1.all[,L] <- curve(sc$Ao1[L],sc$m1[L],sc$r1[L],sc$a1[L],sc$g1[L], 1:simul.comm$gr1.length)
  #  }
  #  A2.all <- matrix(0,nrow=simul.comm$gr2.length, ncol=sc$totS) #response abundances
  #  for(L in 1:sc$totS){
  #    A2.all[,L] <- curve(sc$Ao2[L],sc$m2[L],sc$r2[L],sc$a2[L],sc$g2[L], 1:simul.comm$gr2.length)
  #  } 
  
  A1 <- simul.comm$A1.all[sample.x1, ]
  A2 <- simul.comm$A2.all[sample.x2, ]
  
  p.mat <- A1*A2  # probability of occurrence matrix
  p.mat[is.na (p.mat)] <- 0
  
  a.mat <- p.mat*0 # prepared abundance matrix
  draws.rand <- round(rnorm(Np,mean=no.ind,sd=1))
  
  #output data frames
  samp.out.rand <- matrix(0,nrow=Np,ncol=sc$totS)
  
  #Sampling for random-sample-interval based on individuals
  if (based.on == 'individuals')
    for(i in 1:Np) {
      samp.prob<-p.mat[i,]		#probabilities of sampling each species in given location (based on rel abundance)
      tab.samp <- table(sample(c(1:sc$totS),size=draws.rand[i],prob=samp.prob,replace=T))	#tabulated vector of spp identities after choosing "draws" no. of individuals
      a.mat[i,][as.numeric(names(tab.samp))] <- tab.samp
    } 
  #Sampling for random-sample-interval based on no of species
  if (based.on == 'species')
    for (i in 1:Np) {
      samp.prob <- p.mat [i,]
      spec.pool.size <- sum (samp.prob > 0)
      no.spec <- rnorm (1, k*spec.pool.size)
      tab.samp <- table(sample(c(1:sc$totS),size=no.spec,prob=samp.prob,replace=F))  #tabulated vector of spp identities after choosing no of species
      a.mat[i,][as.numeric(names(tab.samp))] <- tab.samp
    }
  
  colnames (a.mat) <- paste ('spec_', 1:dim (a.mat)[2], sep = '')
  if (pa == T) a.mat[a.mat>0] <- 1		#presence-absence version
  
  result <- list (a.mat = a.mat, p.mat = p.mat, sample.x1 = sample.x1, sample.x2 = sample.x2, sample.comm = list (Np = Np, based.on = based.on, no.ind = if (based.on == 'individuals') no.ind else NULL, k = if (based.on == 'species') k else NULL, seed = seed, pa = pa), simul.comm = simul.comm)
  result
}
