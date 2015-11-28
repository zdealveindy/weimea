#' calculates the fourth-corner statistic (according to Legendre et al. 1997, compiled according to Peres-Neto et al. 2012)
#' @param R Matrix of sample attributes.
#' @param L Matrix of species composition.
#' @param Q Matrix of species attributes
#' @param fc.test Test to be chosen for fourthcorner analysis.
#' @param perm Number of permutations.
fourth.corner <- function (R, L, Q)
{
  R <- as.matrix (R)
  L <- as.matrix (L)
  Q <- as.matrix (Q)
  
  missing.Q <- ifelse (any (is.na (Q)), TRUE, FALSE)
  missing.R <- ifelse (any (is.na (R)), TRUE, FALSE)
  
  if (!missing.Q & !missing.R) r.h <- fourth.corner.0 (R, L, Q)
  if (missing.Q & !missing.R)
    r.h <- apply (Q, 2, FUN = function (q) {L <- L[,!is.na (q)]; q <- q[!is.na (q)]; fourth.corner.0 (R, L, q)})
  if (!missing.Q & missing.R)
    r.h <- t(apply (R, 2, FUN = function (r) {L <- L[!is.na (r),]; r <- r[!is.na (r)]; fourth.corner.0 (r, L, Q)}))
  if (missing.Q & missing.R)
    r.h <- apply (Q, 2, FUN = function (q) apply (R, 2, FUN = function (r) {L <- L[!is.na (r), !is.na (q)]; q <- q[!is.na (q)]; r <- r[!is.na (r)]; fourth.corner.0 (r, L, q)}))
  if (is.null (colnames (r.h))) colnames (r.h) <- colnames (Q)
  if (is.null (rownames (r.h))) rownames (r.h) <- colnames (R)
  class (r.h) <- 'fc'
  return (r.h)
}  


fourth.corner.0 <- function (R, L, Q)  # this can be used only in case of no missing values in Q
{
  W.k <- diag (colSums (L))  
  Q <- as.matrix (Q)
  R <- as.matrix (R)
  L <- as.matrix (L)
  W.n <- diag (rowSums (L))
  I.k <- matrix (rep (1, ncol (L)))
  I.n <- matrix (rep (1, nrow (L)))
  I.g <- matrix (rep (1, ncol (Q)))
  trace.m <- function (M) sum (diag (M))
  
  Q.std.1 <- (Q - I.k %*% t(I.k) %*% W.k %*% Q * (1/trace.m (W.k)))
  Q.std.2 <- (I.k %*% ((t(I.k) %*% W.k %*% (Q * Q) - ((t(I.k) %*% W.k %*% Q) * (t(I.k) %*% W.k %*% Q))/trace.m (W.k)) / trace.m (W.k)) ^ 0.5)
  Q.std <- Q.std.1 / Q.std.2
  
  R.std.1 <- (R - I.n %*% t(I.n) %*% W.n %*% R * (1/trace.m (W.n)))
  R.std.2 <- (I.n %*% ((t(I.n) %*% W.n %*% (R * R) - ((t(I.n) %*% W.n %*% R) * (t(I.n) %*% W.n %*% R))/trace.m (W.n)) / trace.m (W.n)) ^ 0.5)
  R.std <- R.std.1 / R.std.2
  
  #Q.L <- apply (Q.std, 2, FUN = function (Q.std.1) I.n %*% t (Q.std.1) * L)
  X.P <- L %*% Q.std / W.n %*% I.n %*% t(I.g)
  
  r.h <- t(apply (R.std, 2, FUN = function (R.h.std) (t(R.h.std) %*% W.n %*% R.h.std) ^-1 %*% t(R.h.std) %*% W.n %*% X.P))
  colnames (r.h) <- colnames (Q)
  return (r.h)
}

#' @export
print.fc <- function (x, ...) print.default (as.matrix (unclass (x)))

#' @export
summary.fc <- function (object, ...) print (object)

#' @export
fourth.corner.ade <- function (sitspe, speatt, env, fc.test, perm)
{
  R <- as.data.frame (env)
  L <- as.data.frame (sitspe)
  Q <- as.data.frame (speatt)
  
  missing.Q <- ifelse (any (is.na (Q)), TRUE, FALSE)
  missing.R <- ifelse (any (is.na (R)), TRUE, FALSE)
  
  if (!missing.Q & !missing.R) r.h <- fourth.corner.ade0 (R = R, L = L, Q = Q, fc.test = fc.test, perm = perm)
  if (missing.Q & !missing.R){
    r.h <- sapply (names (Q), FUN = function (q) {L1 <- L[,!is.na (Q[,q]), drop = FALSE]; Q1 <- Q[!is.na (Q[,q]), q, drop = FALSE]; fourth.corner.ade0 (R = R, L = L1, Q = Q1, fc.test = fc.test, perm = perm)}, simplify = FALSE)
    r.h <- do.call (rbind.data.frame, r.h)
    rownames (r.h) <- r.h$names
    r.h <- r.h[,-1]
  }
  
  if (!missing.Q & missing.R){
    r.h <- sapply (names (R), FUN = function (r) {L1 <- L[!is.na (R[,r]),, drop = FALSE]; R1 <- R[!is.na (R[,r]), r, drop = FALSE]; fourth.corner.ade0 (R = R1, L = L1, Q = Q, fc.test = fc.test, perm = perm)}, simplify = FALSE)
    r.h <- do.call (rbind.data.frame, r.h)
    rownames (r.h) <- r.h$names
    r.h <- r.h[,-1]
  }
  
  if (missing.Q & missing.R){
    r.h <- sapply (names (Q), FUN = function (q) sapply (names (R), FUN = function (r) {L1 <- L[!is.na (R[,r]), !is.na (Q[,q])]; Q1 <- Q[!is.na (Q[,q]), q, drop = FALSE]; R1 <- R[!is.na (R[,r]), r, drop = FALSE]; fourth.corner.ade0 (R = R1, L = L1, Q = Q1, fc.test = fc.test, perm = perm)}, simplify = FALSE), simplify = FALSE)
    r.h <- do.call (rbind.data.frame, lapply (r.h, FUN = function (y) do.call (rbind.data.frame, y)))
    rownames (r.h) <- r.h$names
    r.h <- r.h[,-1]
  }
  return (r.h)
}


#' @export
fourth.corner.ade0 <- function (R, L, Q, fc.test, perm)
{
  temp <- ade4::fourthcorner (tabR = R, tabL = L, tabQ = Q, modeltype = fc.test, nrepet = perm)
  temp <- do.call (cbind.data.frame, temp$tabD[c('names', 'obs', 'pvalue')])
  names (temp) <-  c('names', 'fourthcorner', 'P.value')
  return (temp)
}
