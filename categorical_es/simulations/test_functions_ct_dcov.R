require(CompQuadForm)

statis.dcov <- function(ct) {
  n <- sum(ct)
  q <- rowSums(ct) / n
  r <- colSums(ct) / n
  e <- n * (q %o% r)
  ts <- sum((ct - e)^2) / n # Important to divide by "n"
  return(ts)
}

coef.dcov <- function(ct) {
  I <- nrow(ct)
  J <- ncol(ct)
  n <- sum(ct)
  q <- rowSums(ct) / n
  r <- colSums(ct) / n
  A <- diag(q) - q %*% t(q)
  B <- diag(r) - r %*% t(r)
  lambda <- eigen(A, T, T)$values
  mu <- eigen(B, T, T)$values
  coef <- as.numeric(lambda %o% mu)
  ac <- abs(coef)
  if (any(ac > 1E-7)) {
    coef <- coef[ac > 1E-9]
  } else {
    coef <- coef[which.max(ac)]
    if (all.equal(0, coef)) coef <- rep(1, (I - 1) * (J - 1))
    # The previous line prevents farebrothers>1
  }
  return(coef)
}

test.dcov <- function(ct) {
  ts <- statis.dcov(ct)
  if (ts < 1E-7) ts <- 1E-7 # Avoids farebrothers > 1
  coef <- coef.dcov(ct)
  pval <- farebrother(ts, coef)$Qq
  if (pval < 0 | pval > 1) cat("Check Farebrother's behaviour.\n")
  return(pval)
}

test.dcov.perm <- function(data, B = 1E4) {
  rowTotals <- rowSums(data)
  colTotals <- colSums(data)
  obs.value <- statis.dcov(ct)
  ts.perm <- vector(mode = "numeric", length = B)
  ts.perm <- sapply(r2dtable(B, rowTotals, colTotals), statis.dcov)
  pval <- mean(ts.perm > obs.value)
  return(pval)
}

statis.chi <- function(ct) {
  n <- sum(ct)
  q <- rowSums(ct) / n
  r <- colSums(ct) / n
  e <- n * (q %o% r)
  if (any(e == 0)) {
    return(NA)
  }
  ts <- sum((ct - e)^2 / e)
  return(ts)
}

test.chi <- function(ct) {
  ts <- statis.chi(ct)
  if (is.na(ts)) {
    return(NA)
  }
  I <- nrow(ct)
  J <- ncol(ct)
  pval <- farebrother(ts, rep(1, (I - 1) * (J - 1)))$Qq
  return(pval)
}

test.chisq <- function(ct) {
  ts <- statis.chi(ct)
  if (is.na(ts)) {
    return(NA)
  }
  I <- nrow(ct)
  J <- ncol(ct)
  pval <- 1 - pchisq(ts, (I - 1) * (J - 1))
  return(pval)
}

statis.g <- function(ct) {
  n <- sum(ct)
  q <- rowSums(ct) / n
  r <- colSums(ct) / n
  e <- n * (q %o% r)
  if (any(e == 0 | ct == 0)) {
    ii <- which(e == 0 | ct == 0)
    # This is also done by DescTools
    ct <- ct[-ii]
    e <- e[-ii]
  }
  ts <- 2 * sum(ct * log(ct / e))
  # ts<-unname( DescTools::GTest(ct)$statistic ) # Slower
  return(ts)
}

test.g <- function(ct) {
  ts <- statis.g(ct)
  if (is.na(ts)) {
    return(NA)
  }
  I <- nrow(ct)
  J <- ncol(ct)
  pval <- 1 - pchisq(ts, (I - 1) * (J - 1))
  # Same result as:
  # pval<-unname( DescTools::GTest(ct,correct="none")$p.value )
  return(pval)
}

test.g.perm <- function(data, B = 1E4) {
  rowTotals <- rowSums(data)
  colTotals <- colSums(data)
  obs.value <- statis.g(ct)
  ts.perm <- vector(mode = "numeric", length = B)
  ts.perm <- sapply(r2dtable(B, rowTotals, colTotals), statis.g)
  pval <- mean(ts.perm > obs.value)
  return(pval)
}

statis.e <- function(obs, p0) {
  n <- sum(obs)
  e <- n * p0
  ts <- sum((obs - e)^2) / n # Important to divide by "n"
  return(ts)
}

coef.e <- function(obs, p0) {
  I <- length(obs)
  n <- sum(obs)
  A <- diag(p0) - p0 %*% t(p0)
  lambda <- eigen(A, T, T)$values
  coef <- lambda
  ac <- abs(coef)
  if (any(ac > 1E-7)) {
    coef <- coef[ac > 1E-9]
  } else {
    coef <- coef[which.max(ac)]
    if (all.equal(0, coef)) coef <- rep(1, I - 1)
    # The previous line prevents farebrothers>1
  }
  return(coef)
}

test.e <- function(obs, p0) {
  ts <- statis.e(obs, p0)
  if (ts < 1E-7) ts <- 1E-7 # Avoids farebrothers > 1
  coef <- coef.e(obs, p0)
  pval <- farebrother(ts, coef)$Qq
  if (pval < 0 | pval > 1) cat("Check Farebrother's behaviour.\n")
  return(pval)
}

statis.chi.gof <- function(obs, p0) {
  n <- sum(obs)
  e <- n * p0
  if (any(e == 0)) {
    return(NA)
  }
  ts <- sum((obs - e)^2 / e)
  return(ts)
}

test.chi.gof <- function(obs, p0, df = NA) {
  ts <- statis.chi.gof(obs, p0)
  if (is.na(ts)) {
    return(NA)
  }
  I <- length(obs)
  if (is.na(df)) df <- I - 1
  pval <- farebrother(ts, rep(1, df))$Qq
  return(pval)
}

test.chisq.gof <- function(obs, p0, df = NA) {
  ts <- statis.chi.gof(obs, p0)
  if (is.na(ts)) {
    return(NA)
  }
  I <- length(obs)
  if (is.na(df)) df <- I - 1
  pval <- 1 - pchisq(ts, df)
  return(pval)
}

statis.usp <- function(ct) {
  n <- sum(ct)
  q <- rowSums(ct) / n
  r <- colSums(ct) / n
  e <- n * (q %o% r)
  ts <- sum((ct - e)^2 / (n * (n - 3)) - 4 * ct * e / (n * (n - 2) * (n - 3)))
  return(ts)
}

test.uspa <- function(ct) {
  ts <- statis.usp(ct)
  coef <- coef.dcov(ct)
  n <- sum(ct)
  kk <- n * ts + sum(coef) # Important to multiply by "n"
  if (kk < 1E-7) kk <- 1E-7 # Avoids farebrothers > 1
  pval <- farebrother(kk, coef)$Qq
  if (pval < 0 | pval > 1) cat("Check Farebrother's behaviour.\n")
  return(pval)
}

test.usp <- function(ct) USPDiscrete(ct, 999)$p.value
