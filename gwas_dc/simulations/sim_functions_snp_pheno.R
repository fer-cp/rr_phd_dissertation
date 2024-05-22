## Eigenvalue
require(coga)
require(reticulate)
# py_install("mpmath")
file <- "pvalue_python.py"
source_python(file)

dyad <- function(number, epsilon = -12) {
  j <- epsilon
  num <- 2^epsilon
  while (num < number) {
    j <- j + 1
    num <- num * 2
  }
  maxj <- j - 1

  vec <- rep(NA, maxj - epsilon + 1)
  dy <- 2^((maxj:epsilon))
  rest <- number
  k <- 1
  for (j in maxj:epsilon) {
    vec[k] <- floor(rest / dy[k])
    rest <- rest %% dy[k]
    k <- k + 1
  }
  return(sum(vec * dy))
}


dwF <- function(x, lambda1, lambda2, v) {
  A <- sqrt(lambda2 / lambda1)
  a <- v * lambda2 / 2
  B01 <- exp((v / 2) * log(a) - ((v + 2) / 2) * log(a + x)) * v / 2
  A * B01 * as.double(hypergeo((v + 2) / 2, 1 / 2, 1, (1 - lambda2 / lambda1) * x / (a + x)))
}


pwF <- function(x, lambda1, lambda2, v, hyp2f1 = "forrey", userflag = 16) {
  # v <- ceiling(v/2)*2
  # v <- 50
  #  lambda2 <- lambda2/v
  #  lambda1 <- lambda1/v
  # x <- x/v
  A <- sqrt(lambda2 / lambda1)
  a <- v * lambda2 / 2
  B01 <- exp((v / 2) * log(a) - ((v + 2) / 2) * log(a + x)) * v / 2 * x
  A * B01 * as.double(appellf1((v + 2) / 2, 1 / 2, 1, 2, (1 - lambda2 / lambda1) * x / (a + x), x / (a + x), hyp2f1 = hyp2f1, userflag = userflag))[3]
}


pvalue_appell <- function(x, lambda1, lambda2, v, digits) {
  a <- v * lambda2 / 2
  xax <- x / (a + x)
  l21 <- lambda2 / lambda1
  pv <- pvalue_python(xax, l21, v, digits)
  if (!is.numeric(pv)) {
    pv <- NA
  }
  return(pv)
}




eigvec <- function(X, b, covariates = NULL) {
  n <- nrow(X)
  q <- ncol(X)
  if (is.null(covariates)) {
    p0 <- sapply(1:q, function(u) length(which(X[, u] == 0))) / n
    p1 <- sapply(1:q, function(u) length(which(X[, u] == 1))) / n
    p2 <- 1 - p1 - p0
    a11 <- b / 2 * (p0 + p2 - (p0 - p2)^2)
    a12sq <- b / 2 * (4 - b) / 2 * ((p2 - p0) * p1) * ((p2 - p0) * p1)
    a22 <- (4 - b) / 2 * (p1 - p1^2)
    Tr <- a11 + a22
    Dt <- a11 * a22 - a12sq
    Term1 <- Tr / 2
    Term2 <- sqrt(Term1^2 - Dt)
    eig1 <- Term1 + Term2
    eig2 <- pmax(Term1 - Term2, 0)
  } else {
    Z <- cbind(rep(1, n), covariates)
    IH <- diag(1, n) - Z %*% solve(t(Z) %*% Z) %*% t(Z)

    U <- sqrt(b / 2) * X
    V <- sqrt((4 - b) / 2) * 1 * (X == 1)

    a11 <- 1 / n * sapply(1:q, function(u) t(U[, u]) %*% IH %*% U[, u])
    a12 <- 1 / n * sapply(1:q, function(u) t(V[, u]) %*% IH %*% U[, u])
    a22 <- 1 / n * sapply(1:q, function(u) t(V[, u]) %*% IH %*% V[, u])
    a12sq <- a12^2

    Tr <- a11 + a22
    Dt <- a11 * a22 - a12sq
    Term1 <- Tr / 2
    Term2 <- sqrt(Term1^2 - Dt)
    eig1 <- Term1 + Term2
    eig2 <- pmax(Term1 - Term2, 0)
  }

  return(list("eig1vec" = eig1, "eig2vec" = eig2))
}

stat <- function(X, Y, b, covariates = NULL) {
  n <- nrow(X)
  het <- 1 * (X == 1)

  if (is.null(covariates)) {
    stat <- 2 * n * (b / 4 * cov(X, Y)^2 + (4 - b) / 4 * cov(het, Y)^2) / var(Y) * (n - 1) / n
  } else {
    Z <- cbind(rep(1, n), covariates)
    IH <- diag(1, n) - Z %*% solve(t(Z) %*% Z) %*% t(Z)
    Yadj <- as.vector(IH %*% Y)
    stat <- 2 * n * (b / 4 * cov(X, Yadj)^2 + (4 - b) / 4 * cov(het, Yadj)^2) / var(Yadj) * (n - 1) / n
  }
  return(stat)
}






screentest <- function(eig, stat, n, covariates = NULL) {
  k <- stat
  lam1 <- eig$eig1vec
  lam2 <- eig$eig2vec

  qcov <- 0
  if (!is.null(covariates)) {
    qcov <- ncol(as.matrix(covariates))
  }

  nadj <- n - qcov


  q <- length(k)
  pv <- rep(NA, q)

  pos <- which(lam2 * n - k > 1e-6)

  neg <- which(lam2 * n - k <= 1e-6)

  kpos <- k[pos]
  kneg <- k[neg]

  lam1pos <- lam1[pos]
  lam2pos <- lam2[pos]
  lam1neg <- lam1[neg]

  consneg <- pf(kneg * (nadj - 3) / (lam1neg * n - kneg), 1, nadj - 3, lower.tail = FALSE)
  acneg <- pf(kneg * (nadj - 2) / (lam1neg * n - kneg), 1, nadj - 2, lower.tail = FALSE)


  # conspos <- 5* (pf(kpos * (nadj-3) / (lam1pos * n + lam2pos * n - 2 * kpos), 1, nadj-3, lower.tail=FALSE))
  conspos <- (pf(kpos * (nadj - 3) / (lam1pos * n + lam2pos * n - 2 * kpos), 1, nadj - 3, lower.tail = FALSE))
  acpos1 <- pf(kpos * (nadj - 3) / (lam1pos * n - kpos), 1, nadj - 3, lower.tail = FALSE)

  acpos2 <- pf(kpos * (nadj - 3) / sqrt((lam1pos * n - kpos) * (lam2pos * n - kpos)), 2, nadj - 3, lower.tail = FALSE)

  qpos <- length(pos)
  acpos3 <- rep(NA, qpos)

  if (qpos > 0) {
    for (j in 1:qpos) {
      acpos3[j] <- 1 - pcoga2dim(kpos[j] * (nadj - 3) / n, 1 / 2, 1 / 2, 1 / (lam1pos[j] - kpos[j] / n) / 2, 1 / (lam2pos[j] - kpos[j] / n) / 2)

      if (acpos3[j] < 1e-14) {
        acpos3[j] <- 0
      }
      #   if (scale2 > 1e-2) {
      #     approx1 <- 1 - pcoga2dim(statcog, 1/2, 1/2 , 1/scale1/2, 1/scale2/2)
      #   } else {
      #     approx1 <- 0
      #   }



      #  if (approx2[j] < 1e-15) {
      #   pv[j] <- approx2[j]
      #  } else {
      #    pv[j] <-max(approx1, approx2[j])
      #  }
    }
  }




  acpos <- pmax(acpos1, acpos2, acpos3)

  p.ac <- rep(NA, q)
  p.co <- rep(NA, q)

  p.ac[pos] <- acpos
  p.ac[neg] <- acneg

  p.co[pos] <- conspos
  p.co[neg] <- consneg

  return(list("p.ac" = p.ac, "p.co" = p.co))
}



.getPold <- function(stat, eig1, eig2, n) {
  pv <- .getP(c(eig1 - stat / n, eig2 - stat / n, rep(-stat / n, n - 3)))
  return(pv)
}

# .getPnew <- function(stat,eig1, eig2 , n) {
#  distobj <- AbscontDistribution(d = function(x) dwF(x, lambda1 = 2*(eig1-stat/n), lambda2 = 2*(eig2 - stat/n), v =n-3),  low1 = 0)
#  pv <- 1-distr::p(distobj)(stat*(n-3)/n)
#  return(pv)
# }


.getPnew <- function(stat, eig1, eig2, n, co, ac, alg, qcov, digits) {
  # distobj <- AbscontDistribution(d = function(x) dwF(x, lambda1 = 2*(eig1-stat/n), lambda2 = 2*(eig2 - stat/n), v =n-3),  low1 = 0)
  lambda1 <- 2 * (eig1 - stat / n)
  lambda2 <- 2 * (eig2 - stat / n)
  if (alg == "appell" & (lambda2 / lambda1) > 1e-2) {
    pv <- try(pvalue_appell(x = stat * (n - 3 - qcov) / n, lambda1 = lambda1, lambda2 = lambda2, v = n - 3 - qcov, digits = digits))
    if (class(pv) == "try-error") {
      pv <- co
    }
    # bad <- class(pv) == "try-error" | pv>co | pv < ac

    #  if (bad) {
    #   pv <-.getP(c(eig1-stat/n, eig2 - stat/n, rep(-stat/n, n-3-qcov)))
    # }
  } else {
    # if (lambda2 < 0) {
    pv <- co
    # } else {
    #  pv <-.getP(c(eig1-stat/n, eig2 - stat/n, rep(-stat/n, n-3-qcov)))
    # }
  }

  return(pv)
}




snptest <- function(X, Y, b = 2, M = 1e-4, m = 1e-64,
                    nchisq = 3e4, alg = "appell",
                    covariates = NULL) {
  # X: genotype matrix
  # Y: pheno vector
  # N>nchisq => No python, but just a Chisq-test
  # alg: we don't change
  # covariates: an n x k matrix (k:=#covariates)
  # Categorical variables: 0's and 1's with dummies
  n <- nrow(X)
  q <- ncol(X)

  qcov <- 0
  if (!is.null(covariates)) {
    qcov <- ncol(as.matrix(covariates))
  }

  stats <- stat(X, Y, b, covariates)
  eigs <- eigvec(X, b, covariates)
  screen <- screentest(eigs, stats, n, covariates)

  p.recalc <- p.recalc2 <- rep(NA, q)

  small <- which(screen$p.co < m)
  p.recalc[small] <- m

  big <- which(screen$p.ac > M)
  p.recalc[big] <- screen$p.ac[big]

  sel <- setdiff(1:q, union(small, big))


  eig1 <- eigs$eig1vec
  eig2 <- eigs$eig2vec

  const <- which(eig1 < 1e-3)
  sel <- setdiff(sel, const)

  if (n < nchisq & b > 0 & b < 4) {
    for (j in sel) {
      p.recalc[j] <- .getPnew(stats[j], eig1[j], eig2[j], n, screen$p.co[j], screen$p.ac[j], alg, qcov, digits = abs(floor(log10(screen$p.ac[j]))))
    }
  } else {
    p.recalc <- screen$p.ac
  }


  return(list("eigen1" = eig1, "eigen2" = eig2, "stats" = stats, "p.co" = screen$p.co, "p.ac" = screen$p.ac, "p.recalc" = p.recalc))
}
