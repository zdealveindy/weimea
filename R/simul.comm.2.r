#' Simulation of ecological data structured by two ecological gradients
#' @author David Zeleny (zeleny.david@@gmail.com) - based on the script of Jason Fridley (Fridley et al. 2007, Appendix S2), modified for two dimensions


#' @param totS total number of species in simulation
#' @param gr1.length,gr2.length length of first and second gradient in units
#' @param niche.type shape of species response curves ('random', 'normal', 'skewed')
#' @param max.niche.breath vector of the length = 2, with maximum niche breath along first and second gradient (default = \code{c(gr1.length, gr2.length)})
#' @param seed set random seed for reproducing always the same result
#' @param plotting should the species response curves along simulated gradient be drawned?
#' @param simul.comm result of simul.comm.2 function with parameters of individual species response curves
#' @param Np number of samples
#' @param sample.x1,sample.x2 positions of sampling along the first (or second, respectively) gradient (default = NULL, means that random locations are generated)
#' @param no.ind mean number of individuals to be drawn from the species pool (if \code{based.on = 'individuals'})
#' @param no.spec number of species to be drawn from species pool (if \code{based.on = 'species'})
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

#' @rdname simul.comm.2
#' @export
simul.comm.2 <- function (totS = 300, gr1.length = 5000, gr2.length = 5000, niche.type = 'random', max.niche.breath = c(gr1.length, gr2.length), seed = NULL, plotting = F)
{
  if (!is.null (seed)) set.seed (seed)
  
  #This is beta function for generating niches
  curve <- function(Ao,m,r,a,g,x)  {
    (Ao*((((x-m)/r)+(a/(a+g)))^a)*((1-(((x-m)/r)+(a/(a+g))))^g))/(((a/(a+g))^a)*((1-(a/(a+g)))^g))
  }
  
  x1 <- seq(1, gr1.length, by = 1) # gradient 1
  x2 <- seq(1, gr2.length, by = 1) # gradient 2
  # species values for random niches
  if(niche.type=="random") {
    S1<-totS					#number of species
    Ao1<-rlnorm(S1,2,1)			#amplitude vector (lognormal distribution)
    m1<-sample(seq(5:max(x1)-5),S1, replace = T)		#location of optima
    r1<-runif(S1,min=10,max=max.niche.breath[1])		#range along gradient (niche breadth)
    a1<-(runif(S1,min=.1,max=4))		#shape parameter (alpha)
    g1<-(runif(S1,min=.1,max=4))		#shape parameter (gamma)
    
    S2<-totS					#number of species
    Ao2<-rlnorm(S2,2,1)			#amplitude vector (lognormal distribution)
    m2<-sample(seq(5:max(x2)-5),S2, replace = T)		#location of optima
    r2<-runif(S2,min=10,max=max.niche.breath[2])		#range along gradient (niche breadth)
    a2<-(runif(S2,min=.1,max=4))		#shape parameter (alpha)
    g2<-(runif(S2,min=.1,max=4))		#shape parameter (gamma)
  }
  
  # species values for normal niches
  if(niche.type=="normal") {
    S1<-totS					#number of species
    Ao1<-rlnorm(S1,2,1)			#amplitude vector (lognormal distribution)
    m1<-sample(seq(5:max(x1)-5),S1, replace = T)		#location of optima
    r1<-runif(S1,min=10,max=max(c(x1,x2)))		#range along gradient (niche breadth)
    a1<-rep(1.99,S1)				#shape parameter (alpha)
    g1<-rep(1.99,S1)				#shape parameter (gamma)
    
    S2<-totS					#number of species
    Ao2<-rlnorm(S2,2,1)			#amplitude vector (lognormal distribution)
    m2<-sample(seq(5:max(x2)-5),S2, replace = T)		#location of optima
    r2<-runif(S2,min=10,max=max(c(x1,x2)))		#range along gradient (niche breadth)
    a2<-rep(1.99,S2)				#shape parameter (alpha)
    g2<-rep(1.99,S2)				#shape parameter (gamma)
  }
  
  # species values for skewed niches
  if(niche.type=="skewed") {
    S1<-totS					#number of species
    Ao1<-rlnorm(S1,2,1)			#amplitude vector (lognormal distribution)
    m1<-sample(seq(5:max(x1)-5),S1, replace = T)		#location of optima
    r1<-runif(S1,min=10,max=max(c(x1,x2)))		#range along gradient (niche breadth)
    #produce skews in either direction, randomly
    a.1<-rep(1.99,S1)
    g.1<-rep(.25,S1)
    samp<-sample(c(1:S1),S1/2,replace=FALSE)
    a.1[samp]<-.25
    g.1[samp]<-1.99
    a1<-a.1				#shape parameter (alpha)
    g1<-g.1				#shape parameter (gamma)
    
    S2<-totS					#number of species
    Ao2<-rlnorm(S2,2,1)			#amplitude vector (lognormal distribution)
    m2<-sample(seq(5:max(x2)-5),S2, replace = T)		#location of optima
    r2<-runif(S2,min=10,max=max(c(x1,x2)))		#range along gradient (niche breadth)
    #produce skews in either direction, randomly
    a.1<-rep(1.99,S2)
    g.1<-rep(.25,S2)
    samp<-sample(c(1:S2),S2/2,replace=FALSE)
    a.1[samp]<-.25
    g.1[samp]<-1.99
    a2<-a.1				#shape parameter (alpha)
    g2<-g.1				#shape parameter (gamma)
  }
  
  if (plotting)
  {
    # Summary plotting:
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
    
    # Plot species distributions along gradient
    par(mfrow=c(2,1))
    
    plot(x1,A1.all[,1],xlim=c(0,max(x1)),ylim=c(0,max(Ao1)),type="l",xlab="Gradient",ylab="Abundance",cex.lab=1.7,cex.axis=1.5,lwd=2)
    for(L in 2:S1) {
      lines(x1,A1.all[,L])
    }
    lines(x1,A1.all[,51],lwd=3,col=2)
    lines(x1,A1.all[,52],lwd=3,col=3)
    
    plot(x2,A2.all[,1],xlim=c(0,max(x2)),ylim=c(0,max(Ao2)),type="l",xlab="Gradient",ylab="Abundance",cex.lab=1.7,cex.axis=1.5,lwd=2)
    for(L in 2:S2) {
      lines(x2,A2.all[,L])
    }
    lines(x2,A2.all[,51],lwd=3,col=2)
    lines(x2,A2.all[,52],lwd=3,col=3)
  }  
  
  result <- list (totS = totS, gr1.length = gr1.length, gr2.length = gr2.length, niche.type = niche.type, Ao1 = Ao1, Ao2 = Ao2, m1 = m1, m2 = m2, r1 = r1, r2 = r2, a1 = a1, a2 = a2, g1 = g1, g2 = g2)
  result
}

#' @rdname simul.comm.2
#' @export
sample.comm.2 <- function (simul.comm = simul.comm.2 (), Np = 300, sample.x1 = NULL, sample.x2 = NULL, no.ind = 100, no.spec = 20, seed = NULL, pa = F, based.on = 'individuals')
{
  if (!is.null (seed)) set.seed (seed)
  sc <- simul.comm
  
  #This is beta function for generating niches
  curve <- function(Ao,m,r,a,g,x)  {
    (Ao*((((x-m)/r)+(a/(a+g)))^a)*((1-(((x-m)/r)+(a/(a+g))))^g))/(((a/(a+g))^a)*((1-(a/(a+g)))^g))
  }
  
  #Random sample intervals along gradient
  if (is.null (sample.x1)) sample.x1 <- trunc(sample(c(2:sc$gr1.length)-1,Np, replace = T))
  if (is.null (sample.x2)) sample.x2 <- trunc(sample(c(2:sc$gr2.length)-1,Np, replace = T))
  
  
  A1 <- matrix(0,nrow=length (sample.x1),ncol=sc$totS) #response abundances
  for(L in 1:sc$totS){
    A1[,L] <- curve(sc$Ao1[L],sc$m1[L],sc$r1[L],sc$a1[L],sc$g1[L], sample.x1)
  }
  A2 <- matrix(0,nrow=length (sample.x2),ncol=sc$totS) #response abundances
  for(L in 1:sc$totS){
    A2[,L] <- curve(sc$Ao2[L],sc$m2[L],sc$r2[L],sc$a2[L],sc$g2[L], sample.x2)
  }
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
      tab.samp <- table(sample(c(1:sc$totS),size=no.spec,prob=samp.prob,replace=F))  #tabulated vector of spp identities after choosing no of species
      a.mat[i,][as.numeric(names(tab.samp))] <- tab.samp
    }
  
  colnames (a.mat) <- paste ('spec_', 1:dim (a.mat)[2], sep = '')
  if (pa == T) a.mat[a.mat>0] <- 1		#presence-absence version
  
  result <- list (a.mat = a.mat, p.mat = p.mat, sample.x1 = sample.x1, sample.x2 = sample.x2, simul.comm = simul.comm)
  result
}
