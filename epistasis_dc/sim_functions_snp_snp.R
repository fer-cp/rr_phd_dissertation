# Here we define a simple hypothesis testing function for epistasis with dcov

require(coga)
require(momentchi2)

eigvec <- function(X,b) {
  # eigenvalues for a matrix of SNPs
  # result: two vectors (one for first eigenvalues, one for second eigenvalues)
  n <- nrow(X)
  q <- ncol(X)
  p0 <- sapply(1:q, function(u) length(which(X[,u]==0)))/n
  p1 <- sapply(1:q, function(u) length(which(X[,u]==1)))/n
  p2 <- 1-p1-p0
  a11 <- b/2*(p0+p2-(p0-p2)^2)   
  a12sq  <- b/2*(4-b)/2 * ((p2-p0)*p1) * ((p2-p0)*p1)
  a22 <- (4-b)/2*(p1-p1^2)
  Tr <- a11+a22
  Dt <- a11*a22-a12sq
  Term1 <- Tr/2
  Term2 <- sqrt(Term1^2 - Dt)
  eig1 <- Term1+Term2
  eig2 <- Term1-Term2
  return(list("eig1vec"=eig1, "eig2vec" = eig2))
}

## Aux function for the testing:
stat1 <- function(X) {
  n <- nrow(X)
  het <- 1*(X==1)
  mix <- cov(X, het)^2
  dcov <- 1/4*cov(X)^2 + 9/4 * cov(het)^2 + 3/4* (mix + t(mix))
  return(n * dcov)
}

## Test for the discrete metric:
test1 <- function(X) {
  eigs <- eigvec(X, b=1)
  eig1 <- eigs$eig1vec
  eig2 <- eigs$eig2vec
  q <- ncol(X)
  st <- stat1(X)
  pv <- matrix(ncol=q,nrow=q)
  for (j in 1:(q-1)) {
    for (k in (j+1):q) {
      eigs <- c(eig1[j]*eig1[k], eig1[j]*eig2[k], eig2[j]*eig1[k], eig2[j]*eig2[k])
      eigs <- eigs[(eigs>1e-8)]
      pv[j,k] <- pv[k,j] <-  1-hbe(eigs,st[j,k])
      #pv[j,k] <- pv[k,j] <- 1-pcoga(st[j,k], rep(1/2,4), 1/eigs/2)
      }
  }
  return(pv)
}