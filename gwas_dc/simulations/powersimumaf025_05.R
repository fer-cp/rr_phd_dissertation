library(parallel)

source("sim_functions_snp_pheno.R")
nsim <- 100000
n <- 300
maf <- 0.25
# h <- 0.5
diff02 <- 1
siglev <- 5e-2
RNGkind("L'Ecuyer-CMRG")
set.seed(1)

pv0 <- pv1 <- pv2 <- pv3 <- pv4 <- pvcat <- pvmax <- pvrec <- rep(NA, nsim)

pvalues <- as.list(rep(NA, 11))

for (k in 1:11) {
  h <- 0.1 * (k - 1)
  pvalues[[k]] <- mclapply(1:nsim, function(j) {
    X <- rbinom(n, 2, maf)
    Y <- diff02 * (1 * (X == 2) + h * 1 * (X == 1)) + 5 * rnorm(n)
    # Y <- rnorm(n)
    # U <- floor(runif(1, 0, 3))
    # Y <- diff02 * (X==U) + 2*rnorm(n)
    Xmat <- cbind(X, rbinom(n, 2, 0.5))

    pv0 <- snptest(Xmat, Y, b = 0, M = 1, m = 0)$p.recalc[1]
    pv1 <- snptest(Xmat, Y, b = 1, M = 1, m = 0)$p.recalc[1]
    pv2 <- snptest(Xmat, Y, b = 2, M = 1, m = 0)$p.recalc[1]
    pv3 <- snptest(Xmat, Y, b = 3, M = 1, m = 0)$p.recalc[1]
    pv4 <- snptest(Xmat, Y, b = 4, M = 1, m = 0)$p.recalc[1]

    facX <- as.factor(X)
    fit <- lm(Y ~ facX)
    pvcat <- anova(fit)[[5]][1]


    pvmax <- try(nmax3(Y, X)$p.value)
    if (class(pvmax) == "try-error") {
      pvmax <- 1
    }

    # rec <- 1*(X==2)
    # fit <- lm(Y ~ rec)
    # pvrec <- anova(fit)[[5]][1]

    return(list(pv0, pv1, pv2, pv3, pv4, pvmax, pvcat))
  }, mc.cores = 50)
}



pow0 <- pow1 <- pow2 <- pow3 <- pow4 <- powcat <- powmax <- powrec <- rep(NA, 11)
siglev <- 5e-2
for (k in 1:11) {
  pvalues_mat <- matrix(unlist(pvalues[[k]]), byrow = TRUE, ncol = 7)
  pow0[k] <- length(which(pvalues_mat[, 1] < siglev)) / nsim
  pow1[k] <- length(which(pvalues_mat[, 2] < siglev)) / nsim
  pow2[k] <- length(which(pvalues_mat[, 3] < siglev)) / nsim
  pow3[k] <- length(which(pvalues_mat[, 4] < siglev)) / nsim
  pow4[k] <- length(which(pvalues_mat[, 5] < siglev)) / nsim
  powmax[k] <- length(which(pvalues_mat[, 6] < siglev)) / nsim
  powcat[k] <- length(which(pvalues_mat[, 7] < siglev)) / nsim
  # powrec[k] <- length(which(pvalues_mat[,8] < siglev ))/ nsim
}
xax <- 0.1 * (0:10)
plot(xax, pow3, ylim = c(0, 0.5), type = "l", col = "blue")
lines(xax, pow4, ylim = c(0, 1), type = "l", col = "black")
lines(xax, pow2, ylim = c(0, 1), type = "l", col = "green")
lines(xax, powcat, ylim = c(0, 1), type = "l", col = "goldenrod")
lines(xax, powmax, ylim = c(0, 1), type = "l", col = "purple")

save.image("maf05_025_300.RData")
