library(parallel)
library(AssocTests)
source("sim_functions_snp_pheno.R")
nsim <- 100000
n <- 300
maf <- 0.5
# h <- 0.5
diff02 <- 0
siglev <- 5e-2
RNGkind("L'Ecuyer-CMRG")
set.seed(1)


pv0 <- pv1 <- pv2 <- pv3 <- pv4 <- pvcat <- pvmax <- pvrec <- rep(NA, nsim)

pvalues <- as.list(rep(NA, 11))

for (k in 1:5) {
  pvalues[[k]] <- mclapply(1:nsim, function(j) {
    maf <- 0.1 * k
    X <- rbinom(n, 2, maf)
    Y <- rnorm(n)
    Xmat <- cbind(X, rbinom(n, 2, 0.5))

    pv0 <- snptest(Xmat, Y, b = 0, M = 1, m = 0)$p.recalc[1]
    pv1 <- snptest(Xmat, Y, b = 1, M = 1, m = 0)$p.recalc[1]
    pv2 <- snptest(Xmat, Y, b = 2, M = 1, m = 0)$p.recalc[1]
    pv3 <- snptest(Xmat, Y, b = 3, M = 1, m = 0)$p.recalc[1]
    pv4 <- snptest(Xmat, Y, b = 4, M = 1, m = 0)$p.recalc[1]

    facX <- as.factor(X)
    fit <- lm(Y ~ facX)
    pvcat <- try(anova(fit)[[5]][1])
    if (class(pvcat) == "try-error") pvcat <- 1

    pvmax <- try(nmax3(Y, X)$p)
    if (class(pvmax) == "try-error") pvmax <- 1

    # rec <- 1*(X==2)
    # fit <- lm(Y ~ rec)
    # pvrec <- anova(fit)[[5]][1]
    pvrec <- 1

    return(list(pv0, pv1, pv2, pv3, pv4, pvmax, pvcat, pvrec))
  }, mc.cores = 50)
}



pow0 <- pow1 <- pow2 <- pow3 <- pow4 <- powcat <- powmax <- powrec <- rep(NA, 5)
siglev <- 5e-2
for (k in 1:5) {
  pvalues_mat <- matrix(unlist(pvalues[[k]]), byrow = TRUE, ncol = 8)
  pow0[k] <- length(which(pvalues_mat[, 1] < siglev)) / nsim
  pow1[k] <- length(which(pvalues_mat[, 2] < siglev)) / nsim
  pow2[k] <- length(which(pvalues_mat[, 3] < siglev)) / nsim
  pow3[k] <- length(which(pvalues_mat[, 4] < siglev)) / nsim
  pow4[k] <- length(which(pvalues_mat[, 5] < siglev)) / nsim
  powmax[k] <- length(which(pvalues_mat[, 6] < siglev)) / nsim
  powcat[k] <- length(which(pvalues_mat[, 7] < siglev)) / nsim
  powrec[k] <- length(which(pvalues_mat[, 8] < siglev)) / nsim
}
xax <- 0.1 * (1:5)
plot(xax, pow3, ylim = c(0, 0.01), type = "l", col = "blue")
lines(xax, pow4, ylim = c(0, 1), type = "l", col = "black")
lines(xax, pow2, ylim = c(0, 1), type = "l", col = "green")
lines(xax, powcat, ylim = c(0, 1), type = "l", col = "goldenrod")
lines(xax, powmax, ylim = c(0, 1), type = "l", col = "purple")
abline(a = 5e-3, b = 0, lty = 2)
save.image("typeI_005.RData")
# df <- data.frame("maf" = rep(xax,5), "power" = c(pow2,pow3,pow4,powcat,powmax), "method" = c(rep("b2",5),rep("b3",5),rep("add",5), rep("anova",5), rep("nmax3",5)))
