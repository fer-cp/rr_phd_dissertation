library(parallel)
library(AssocTests)
library(microbenchmark)
source("sim_functions_snp_pheno.R")
nsim <- 100
# n <- 2000
maf <- 0.5
p <- 100e3
Xmat <- matrix(rbinom(n * p, 2, maf), ncol = p)
Y <- rnorm(n)
RNGkind("L'Ecuyer-CMRG")
set.seed(1)


pv0 <- pv1 <- pv2 <- pv3 <- pv4 <- pvcat <- pvmax <- pvrec <- rep(NA, nsim)

comptime <- as.list(rep(NA, 5))

for (k in 1:5) {
  comptime[[k]] <- mclapply(1:nsim, function(u) {
    n <- 2000 * k
    Xmat <- matrix(rbinom(n * p, 2, maf), ncol = p)
    Y <- rnorm(n)

    pv2naiv <- microbenchmark(snptest(Xmat, Y, b = 2, M = 1, m = 0), times = 1)$time
    pv2fast <- microbenchmark(snptest(Xmat, Y, b = 2, M = 1e-3, m = 0), times = 1)$time

    pv3naiv <- microbenchmark(snptest(Xmat, Y, b = 3, M = 1, m = 0), times = 1)$time
    pv3fast <- microbenchmark(snptest(Xmat, Y, b = 3, M = 1e-3, m = 0), times = 1)$time

    pv4 <- microbenchmark(snptest(Xmat, Y, b = 4, M = 1, m = 0), times = 1)$time

    pvmax <- 100 * microbenchmark(for (i in 1:1000) {
      nmax3(Y, Xmat[, i])
    }, times = 1)$time
    return(list(pv2naiv, pv2fast, pv3naiv, pv3fast, pv4, pvmax))
  }, mc.cores = 50)
}

save.image("comptime.RData")



pv2naiv <- pv2fast <- pv3naiv <- pv3fast <- pv4 <- pvmax <- matrix(ncol = 5, nrow = 100)

for (k in 1:5) {
  for (j in 1:100) {
    if (is.numeric(comptime[[k]][[j]][[1]])) {
      pv2naiv[j, k] <- comptime[[k]][[j]][[1]]
      pv2fast[j, k] <- comptime[[k]][[j]][[2]]
      pv3naiv[j, k] <- comptime[[k]][[j]][[3]]
      pv3fast[j, k] <- comptime[[k]][[j]][[4]]
      pv4[j, k] <- comptime[[k]][[j]][[5]]
      pvmax[j, k] <- comptime[[k]][[j]][[6]]
    }
  }
}

mean_na <- function(x) {
  mean(x, na.rm = TRUE)
}

timeb2naiv <- apply(pv2naiv, 2, mean_na) / 60e9
timeb2fast <- apply(pv2fast, 2, mean_na) / 60e9
timeb3naiv <- apply(pv3naiv, 2, mean_na) / 60e9
timeb3fast <- apply(pv3fast, 2, mean_na) / 60e9
timeb4 <- apply(pv4, 2, mean_na) / 60e9
timemax <- apply(pvmax, 2, mean_na) / 60e9

save.image("comptime.RData")



png("comptime.png", width = 2400 / 75 * 240, height = 1200 / 75 * 240, res = 300)
# setwd("/users/edelmand/documents")
# load("sim_lowcens.RData")
par(mfrow = c(1, 2))
par(mar = c(5, 6, 4, 1) + .1)
plot(ss, timeb2fast, type = "b", ylim = c(0, 3), ylab = "Time in Minutes", xlab = "Sample Size", col = "blue", cex = 2, pch = 20, cex.lab = 2, cex.axis = 2, cex.main = 2, cex.sub = 2, lwd = 2)
points(ss, timeb3fast, type = "b", pch = 15, col = "red", lwd = 2, cex = 2, lty = 2)
points(ss, timeb4, type = "b", pch = 17, col = "black", lwd = 2, cex = 2, lty = 3)
legend("bottomright",
  legend = c("b=2", "b=3", "add. (b=4)"),
  pch = c(20, 15, 17), col = c("blue", "red", "black"), bty = "n", cex = 2
)

plot(ss, log10(timeb2fast), type = "b", ylim = c(-1, 4), ylab = "Time in Minutes", yaxt = "n", xlab = "Sample Size", col = "blue", cex = 2, pch = 20, cex.lab = 2, cex.axis = 2, cex.main = 2, cex.sub = 2, lwd = 2)
points(ss, log10(timeb3fast), type = "b", pch = 15, col = "red", lwd = 2, cex = 2, lty = 2)
points(ss, log10(timeb4), type = "b", pch = 17, col = "black", lwd = 2, cex = 2, lty = 3)
points(ss, log10(timeb2naiv), type = "b", pch = 12, col = "darkgreen", lwd = 2, cex = 2, lty = 3)
points(ss, log10(timeb3naiv), type = "b", pch = 10, col = "goldenrod", lwd = 2, cex = 2, lty = 3)
points(ss, log10(timemax), type = "b", pch = 8, col = "purple", lwd = 2, cex = 2, lty = 3)
axis(2, at = -1:4, labels = cumprod(rep(10, 6)) * 0.01, cex.axis = 2)
legend("bottomright",
  legend = c("b=2", "b=3", "add. (b=4)", "b=2 (w/o screen)", "b=3 (w/o screen)", "nmax3"),
  pch = c(20, 15, 17, 12, 10, 8), col = c("blue", "red", "black", "darkgreen", "goldenrod", "purple"), bty = "n", cex = 2
)

dev.off()
