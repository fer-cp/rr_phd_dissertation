## Computation Time

load("comptime.RData")

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


png("comptime.png", width = 2400 / 75 * 240, height = 1200 / 75 * 240, res = 300)
# setwd("/users/edelmand/documents")
# load("sim_lowcens.RData")
par(mfrow = c(1, 2))
par(mar = c(5, 6, 4, 1) + .1)
plot(ss, timeb2fast, type = "b", ylim = c(0, 3), ylab = "Time in Minutes", xlab = "Sample Size", col = "blue", cex = 2, pch = 20, cex.lab = 2, cex.axis = 2, cex.main = 2, cex.sub = 2, lwd = 2)
points(ss, timeb3fast, type = "b", pch = 15, col = "red", lwd = 2, cex = 2, lty = 2)
points(ss, timeb4, type = "b", pch = 17, col = "black", lwd = 2, cex = 2, lty = 3)
legend("bottomright",
  legend = c("b=2 (fast)", "b=3 (fast)", "add. (b=4)"),
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
  legend = c("b=2 (fast)", "b=3 (fast)", "add. (b=4)", "b=2 (naiv)", "b=3 (naiv)", "nmax3"),
  pch = c(20, 15, 17, 12, 10, 8), col = c("blue", "red", "black", "darkgreen", "goldenrod", "purple"), bty = "n", cex = 2
)

dev.off()

## Type I Error

nsim <- 1e5
load("typeI_005.RData")
t.05.2 <- t.05.3 <- t.05.4 <- t.05.max <- t.05.cat <- rep(NA, 5)
siglev <- 5e-2
for (k in 1:5) {
  pvalues_mat <- matrix(unlist(pvalues[[k]]), byrow = TRUE, ncol = 8)
  t.05.2[k] <- length(which(pvalues_mat[, 3] < siglev)) / nsim
  t.05.3[k] <- length(which(pvalues_mat[, 4] < siglev)) / nsim
  t.05.4[k] <- length(which(pvalues_mat[, 5] < siglev)) / nsim
  t.05.max[k] <- length(which(pvalues_mat[, 6] < siglev)) / nsim
  t.05.cat[k] <- length(which(pvalues_mat[, 7] < siglev)) / nsim
}


load("typeI_5e5.RData")
nsim <- 1e8
t.5e5.2 <- t.5e5.3 <- t.5e5.4 <- rep(NA, 5)
siglev <- 5e-5
for (k in 1:5) {
  pvalues_vec <- unlist(pvalues[[k]])
  pvalues0 <- unlist(lapply(0:9999, function(u) {
    pvalues_vec[(1:10000) + 5e4 * u]
  }))
  pvalues1 <- unlist(lapply(0:9999, function(u) {
    pvalues_vec[1e4 + (1:10000) + 5e4 * u]
  }))
  pvalues2 <- unlist(lapply(0:9999, function(u) {
    pvalues_vec[2e4 + (1:10000) + 5e4 * u]
  }))
  pvalues3 <- unlist(lapply(0:9999, function(u) {
    pvalues_vec[3e4 + (1:10000) + 5e4 * u]
  }))
  pvalues4 <- unlist(lapply(0:9999, function(u) {
    pvalues_vec[4e4 + (1:10000) + 5e4 * u]
  }))

  t.5e5.2[k] <- length(which(pvalues2 < siglev)) / nsim
  t.5e5.3[k] <- length(which(pvalues3 < siglev)) / nsim
  t.5e5.4[k] <- length(which(pvalues4 < siglev)) / nsim
}

xax <- 0.1 * (1:5)
png("typeI.png", width = 2400 / 75 * 240, height = 1200 / 75 * 240, res = 300)
par(mfrow = c(1, 2))
par(mar = c(5, 6, 4, 1) + .1)
plot(xax, t.05.2, type = "b", ylim = c(0, 0.1), xlab = "Minor Allele Frequency (MAF)", ylab = "Type I Error", col = "blue", cex = 2, pch = 20, cex.lab = 2, cex.axis = 2, cex.main = 2, cex.sub = 2, lwd = 2)
points(xax, t.05.3, type = "b", pch = 15, col = "red", lwd = 2, cex = 2, lty = 2)
points(xax, t.05.4, type = "b", pch = 17, col = "black", lwd = 2, cex = 2, lty = 3)
points(xax, t.05.cat, type = "b", pch = 10, col = "goldenrod", lwd = 2, cex = 2, lty = 3)
points(xax, t.05.max, type = "b", pch = 8, col = "purple", lwd = 2, cex = 2, lty = 3)
legend("bottomright",
  legend = c("b=2", "b=3", "add. (b=4)", "categorical (ANOVA)", "nmax3"),
  pch = c(20, 15, 17, 10, 8), col = c("blue", "red", "black", "goldenrod", "purple"), bty = "n", cex = 2
)
abline(a = 5e-2, b = 0, lty = 2)
plot(xax, t.5e5.2, type = "b", ylim = c(0, 1e-4), xlab = "Minor Allele Frequency (MAF)", ylab = "Type I Error", col = "blue", cex = 2, pch = 20, cex.lab = 2, cex.axis = 2, cex.main = 2, cex.sub = 2, lwd = 2)
points(xax, t.5e5.3, type = "b", pch = 15, col = "red", lwd = 2, cex = 2, lty = 2)
points(xax, t.5e5.4, type = "b", pch = 17, col = "black", lwd = 2, cex = 2, lty = 3)
legend("bottomright",
  legend = c("b=2", "b=3", "add. (b=4)"),
  pch = c(20, 15, 17), col = c("blue", "red", "black"), bty = "n", cex = 2
)
abline(a = 5e-5, b = 0, lty = 2)

dev.off()



## Power fuer alpha 0.05

load("maf05_01_300.RData")
pow2.01 <- pow2
pow3.01 <- pow3
pow4.01 <- pow4
powcat.01 <- powcat
powmax.01 <- powmax

load("maf05_03_300.RData")
pow2.03 <- pow2
pow3.03 <- pow3
pow4.03 <- pow4
powcat.03 <- powcat
powmax.03 <- powmax


load("maf05_05_300.RData")
xax <- 0.1 * (0:10)
pow2.05 <- pow2
pow3.05 <- pow3
pow4.05 <- pow4
powcat.05 <- powcat
powmax.05 <- powmax


xax <- 0.1 * (0:10)
png("power05.png", width = 2400 / 75 * 240, height = 1200 / 75 * 240, res = 300)
par(mfrow = c(1, 3))
par(mar = c(12, 6, 4, 1) + .1)
plot(xax, pow2.01, type = "b", ylim = c(0, 0.5), xlab = "Heterozygous effect (h)", ylab = "Power", main = "MAF 0.1", col = "blue", cex = 3, pch = 20, cex.lab = 3, cex.axis = 3, cex.main = 3, cex.sub = 3, lwd = 3)
points(xax, pow3.01, type = "b", pch = 15, col = "red", lwd = 3, cex = 3, lty = 2)
points(xax, pow4.01, type = "b", pch = 17, col = "black", lwd = 3, cex = 3, lty = 3)
points(xax, powcat.01, type = "b", pch = 10, col = "goldenrod", lwd = 3, cex = 3, lty = 3)
points(xax, powmax.01, type = "b", pch = 8, col = "purple", lwd = 3, cex = 3, lty = 3)
legend("topleft",
  legend = c("b=2", "b=3", "add. (b=4)", "categorical (ANOVA)", "nmax3"),
  pch = c(20, 15, 17, 10, 8), col = c("blue", "red", "black", "goldenrod", "purple"), bty = "n", cex = 4
)
abline(a = 5e-2, b = 0, lty = 2)

plot(xax, pow2.03, type = "b", ylim = c(0, 0.5), xlab = "Heterozygous effect (h)", ylab = "Power", main = "MAF 0.3", col = "blue", cex = 3, pch = 20, cex.lab = 3, cex.axis = 3, cex.main = 3, cex.sub = 3, lwd = 3)
points(xax, pow3.03, type = "b", pch = 15, col = "red", lwd = 3, cex = 3, lty = 2)
points(xax, pow4.03, type = "b", pch = 17, col = "black", lwd = 3, cex = 3, lty = 3)
points(xax, powcat.03, type = "b", pch = 10, col = "goldenrod", lwd = 3, cex = 3, lty = 3)
points(xax, powmax.03, type = "b", pch = 8, col = "purple", lwd = 3, cex = 3, lty = 3)
legend("topleft",
  legend = c("b=2", "b=3", "add. (b=4)", "categorical (ANOVA)", "nmax3"),
  pch = c(20, 15, 17, 10, 8), col = c("blue", "red", "black", "goldenrod", "purple"), bty = "n", cex = 4
)
abline(a = 5e-2, b = 0, lty = 2)

plot(xax, pow2.05, type = "b", ylim = c(0, 0.5), xlab = "Heterozygous effect (h)", ylab = "Power", main = "MAF 0.5", col = "blue", cex = 3, pch = 20, cex.lab = 3, cex.axis = 3, cex.main = 3, cex.sub = 3, lwd = 3)
points(xax, pow3.05, type = "b", pch = 15, col = "red", lwd = 3, cex = 3, lty = 2)
points(xax, pow4.05, type = "b", pch = 17, col = "black", lwd = 3, cex = 3, lty = 3)
points(xax, powcat.05, type = "b", pch = 10, col = "goldenrod", lwd = 3, cex = 3, lty = 3)
points(xax, powmax.05, type = "b", pch = 8, col = "purple", lwd = 3, cex = 3, lty = 3)
legend("topleft",
  legend = c("b=2", "b=3", "add. (b=4)", "categorical (ANOVA)", "nmax3"),
  pch = c(20, 15, 17, 10, 8), col = c("blue", "red", "black", "goldenrod", "purple"), bty = "n", cex = 4
)
abline(a = 5e-2, b = 0, lty = 2)

dev.off()





## Power fuer alpha 5e-8 - Version 1

load("maf01_5e8.RData")
pow0 <- pow1 <- pow2 <- pow3 <- pow4 <- powcat <- powmax <- powrec <- rep(NA, 11)
siglev <- 5e-8
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
pow2.001 <- pow2
pow3.001 <- pow3
pow4.001 <- pow4
powcat.001 <- powcat
powmax.001 <- powmax

load("maf025_5e8.RData")
pow0 <- pow1 <- pow2 <- pow3 <- pow4 <- powcat <- powmax <- powrec <- rep(NA, 11)
siglev <- 5e-8
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
pow2.01 <- pow2
pow3.01 <- pow3
pow4.01 <- pow4
powcat.01 <- powcat
powmax.01 <- powmax


load("maf04_5e8.RData")
pow0 <- pow1 <- pow2 <- pow3 <- pow4 <- powcat <- powmax <- powrec <- rep(NA, 11)
siglev <- 5e-8
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
pow2.03 <- pow2
pow3.03 <- pow3
pow4.03 <- pow4
powcat.03 <- powcat
powmax.03 <- powmax



load("maf05_5e8.RData")
pow0 <- pow1 <- pow2 <- pow3 <- pow4 <- powcat <- powmax <- powrec <- rep(NA, 11)
siglev <- 5e-8
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
pow2.05 <- pow2
pow3.05 <- pow3
pow4.05 <- pow4
powcat.05 <- powcat
powmax.05 <- powmax




xax <- 0.1 * (0:10)
png("power5e8.png", width = 2400 / 75 * 240, height = 1200 / 75 * 240 * 2.2, res = 300)
par(mfrow = c(2, 2))
par(mar = c(5, 6, 4, 1) + .1)

plot(xax, pow2.001, type = "b", ylim = c(0, 0.5), xlab = "Heterozygous effect (h)", ylab = "Power", main = "MAF 0.1", col = "blue", cex = 3, pch = 20, cex.lab = 3, cex.axis = 3, cex.main = 3, cex.sub = 3, lwd = 3)
points(xax, pow3.001, type = "b", pch = 15, col = "red", lwd = 3, cex = 3, lty = 2)
points(xax, pow4.001, type = "b", pch = 17, col = "black", lwd = 3, cex = 3, lty = 3)
points(xax, powcat.001, type = "b", pch = 10, col = "goldenrod", lwd = 3, cex = 3, lty = 3)
points(xax, powmax.001, type = "b", pch = 8, col = "purple", lwd = 3, cex = 3, lty = 3)
legend("topleft",
  legend = c("b=2", "b=3", "add. (b=4)", "categorical (ANOVA)", "nmax3"),
  pch = c(20, 15, 17, 10, 8), col = c("blue", "red", "black", "goldenrod", "purple"), bty = "n", cex = 4
)
abline(a = 5e-8, b = 0, lty = 2)


plot(xax, pow2.01, type = "b", ylim = c(0, 0.5), xlab = "Heterozygous effect (h)", ylab = "Power", main = "MAF 0.25", col = "blue", cex = 3, pch = 20, cex.lab = 3, cex.axis = 3, cex.main = 3, cex.sub = 3, lwd = 3)
points(xax, pow3.01, type = "b", pch = 15, col = "red", lwd = 3, cex = 3, lty = 2)
points(xax, pow4.01, type = "b", pch = 17, col = "black", lwd = 3, cex = 3, lty = 3)
points(xax, powcat.01, type = "b", pch = 10, col = "goldenrod", lwd = 3, cex = 3, lty = 3)
points(xax, powmax.01, type = "b", pch = 8, col = "purple", lwd = 3, cex = 3, lty = 3)
legend("topleft",
  legend = c("b=2", "b=3", "add. (b=4)", "categorical (ANOVA)", "nmax3"),
  pch = c(20, 15, 17, 10, 8), col = c("blue", "red", "black", "goldenrod", "purple"), bty = "n", cex = 4
)
abline(a = 5e-8, b = 0, lty = 2)

plot(xax, pow2.03, type = "b", ylim = c(0, 0.5), xlab = "Heterozygous effect (h)", ylab = "Power", main = "MAF 0.4", col = "blue", pch = 20, cex = 3, cex.lab = 3, cex.axis = 3, cex.main = 3, cex.sub = 3, lwd = 3)
points(xax, pow3.03, type = "b", pch = 15, col = "red", lwd = 3, cex = 3, lty = 2)
points(xax, pow4.03, type = "b", pch = 17, col = "black", lwd = 3, cex = 3, lty = 3)
points(xax, powcat.03, type = "b", pch = 10, col = "goldenrod", lwd = 3, cex = 3, lty = 3)
points(xax, powmax.03, type = "b", pch = 8, col = "purple", lwd = 3, cex = 3, lty = 3)
legend("topleft",
  legend = c("b=2", "b=3", "add. (b=4)", "categorical (ANOVA)", "nmax3"),
  pch = c(20, 15, 17, 10, 8), col = c("blue", "red", "black", "goldenrod", "purple"), bty = "n", cex = 4
)
abline(a = 5e-8, b = 0, lty = 2)

plot(xax, pow2.05, type = "b", ylim = c(0, 0.5), xlab = "Heterozygous effect (h)", ylab = "Power", main = "MAF 0.5", col = "blue", pch = 20, cex = 3, cex.lab = 3, cex.axis = 3, cex.main = 3, cex.sub = 3, lwd = 3)
points(xax, pow3.05, type = "b", pch = 15, col = "red", lwd = 3, cex = 3, lty = 2)
points(xax, pow4.05, type = "b", pch = 17, col = "black", lwd = 3, cex = 3, lty = 3)
points(xax, powcat.05, type = "b", pch = 10, col = "goldenrod", lwd = 2, cex = 3, lty = 3)
points(xax, powmax.05, type = "b", pch = 8, col = "purple", lwd = 3, cex = 3, lty = 3)
legend("topleft",
  legend = c("b=2", "b=3", "add. (b=4)", "categorical (ANOVA)", "nmax3"),
  pch = c(20, 15, 17, 10, 8), col = c("blue", "red", "black", "goldenrod", "purple"), bty = "n", cex = 4
)
abline(a = 5e-8, b = 0, lty = 2)

dev.off()







## Power fuer alpha 0.05, Version 2

load("maf05_01_300.RData")
pow2.001 <- pow2
pow3.001 <- pow3
pow4.001 <- pow4
powcat.001 <- powcat
powmax.001 <- powmax

load("maf05_025_300.RData")
pow2.01 <- pow2
pow3.01 <- pow3
pow4.01 <- pow4
powcat.01 <- powcat
powmax.01 <- powmax

load("maf05_04_300.RData")
pow2.03 <- pow2
pow3.03 <- pow3
pow4.03 <- pow4
powcat.03 <- powcat
powmax.03 <- powmax


load("maf05_05_300.RData")
xax <- 0.1 * (0:10)
pow2.05 <- pow2
pow3.05 <- pow3
pow4.05 <- pow4
powcat.05 <- powcat
powmax.05 <- powmax




xax <- 0.1 * (0:10)
png("power005.png", width = 2400 / 75 * 240, height = 1200 / 75 * 240 * 2.2, res = 300)
par(mfrow = c(2, 2))
par(mar = c(5, 6, 4, 1) + .1)

plot(xax, pow2.001, type = "b", ylim = c(0, 0.5), xlab = "Heterozygous effect (h)", ylab = "Power", main = "MAF 0.1", col = "blue", cex = 3, pch = 20, cex.lab = 3, cex.axis = 3, cex.main = 3, cex.sub = 3, lwd = 3)
points(xax, pow3.001, type = "b", pch = 15, col = "red", lwd = 3, cex = 3, lty = 2)
points(xax, pow4.001, type = "b", pch = 17, col = "black", lwd = 3, cex = 3, lty = 3)
points(xax, powcat.001, type = "b", pch = 10, col = "goldenrod", lwd = 3, cex = 3, lty = 3)
points(xax, powmax.001, type = "b", pch = 8, col = "purple", lwd = 3, cex = 3, lty = 3)
legend("topleft",
  legend = c("b=2", "b=3", "add. (b=4)", "categorical (ANOVA)", "nmax3"),
  pch = c(20, 15, 17, 10, 8), col = c("blue", "red", "black", "goldenrod", "purple"), bty = "n", cex = 4
)
abline(a = 5e-2, b = 0, lty = 2)


plot(xax, pow2.01, type = "b", ylim = c(0, 0.5), xlab = "Heterozygous effect (h)", ylab = "Power", main = "MAF 0.25", col = "blue", cex = 3, pch = 20, cex.lab = 3, cex.axis = 3, cex.main = 3, cex.sub = 3, lwd = 3)
points(xax, pow3.01, type = "b", pch = 15, col = "red", lwd = 3, cex = 3, lty = 2)
points(xax, pow4.01, type = "b", pch = 17, col = "black", lwd = 3, cex = 3, lty = 3)
points(xax, powcat.01, type = "b", pch = 10, col = "goldenrod", lwd = 3, cex = 3, lty = 3)
points(xax, powmax.01, type = "b", pch = 8, col = "purple", lwd = 3, cex = 3, lty = 3)
legend("topleft",
  legend = c("b=2", "b=3", "add. (b=4)", "categorical (ANOVA)", "nmax3"),
  pch = c(20, 15, 17, 10, 8), col = c("blue", "red", "black", "goldenrod", "purple"), bty = "n", cex = 4
)
abline(a = 5e-2, b = 0, lty = 2)

plot(xax, pow2.03, type = "b", ylim = c(0, 0.5), xlab = "Heterozygous effect (h)", ylab = "Power", main = "MAF 0.4", col = "blue", pch = 20, cex = 3, cex.lab = 3, cex.axis = 3, cex.main = 3, cex.sub = 3, lwd = 3)
points(xax, pow3.03, type = "b", pch = 15, col = "red", lwd = 3, cex = 3, lty = 2)
points(xax, pow4.03, type = "b", pch = 17, col = "black", lwd = 3, cex = 3, lty = 3)
points(xax, powcat.03, type = "b", pch = 10, col = "goldenrod", lwd = 3, cex = 3, lty = 3)
points(xax, powmax.03, type = "b", pch = 8, col = "purple", lwd = 3, cex = 3, lty = 3)
legend("topleft",
  legend = c("b=2", "b=3", "add. (b=4)", "categorical (ANOVA)", "nmax3"),
  pch = c(20, 15, 17, 10, 8), col = c("blue", "red", "black", "goldenrod", "purple"), bty = "n", cex = 4
)
abline(a = 5e-2, b = 0, lty = 2)

plot(xax, pow2.05, type = "b", ylim = c(0, 0.5), xlab = "Heterozygous effect (h)", ylab = "Power", main = "MAF 0.5", col = "blue", pch = 20, cex = 3, cex.lab = 3, cex.axis = 3, cex.main = 3, cex.sub = 3, lwd = 3)
points(xax, pow3.05, type = "b", pch = 15, col = "red", lwd = 3, cex = 3, lty = 2)
points(xax, pow4.05, type = "b", pch = 17, col = "black", lwd = 3, cex = 3, lty = 3)
points(xax, powcat.05, type = "b", pch = 10, col = "goldenrod", lwd = 2, cex = 3, lty = 3)
points(xax, powmax.05, type = "b", pch = 8, col = "purple", lwd = 3, cex = 3, lty = 3)
legend("topleft",
  legend = c("b=2", "b=3", "add. (b=4)", "categorical (ANOVA)", "nmax3"),
  pch = c(20, 15, 17, 10, 8), col = c("blue", "red", "black", "goldenrod", "purple"), bty = "n", cex = 4
)
abline(a = 5e-2, b = 0, lty = 2)

dev.off()



## Rare SNPs





load("maf001_5e8_30k.RData")
pow2.03 <- pow2
pow3.03 <- pow3
pow4.03 <- pow4
powcat.03 <- powcat
powmax.03 <- powmax


load("maf005_5e8_6k.RData")
xax <- 0.1 * (0:10)
pow2.05 <- pow2
pow3.05 <- pow3
pow4.05 <- pow4
powcat.05 <- powcat
powmax.05 <- powmax




xax <- 0.1 * (0:10)
png("rare.png", width = 2400 / 75 * 240, height = 1200 / 75 * 240, res = 300)
par(mfrow = c(1, 2))
par(mar = c(5, 6, 4, 1) + .1)

plot(xax, pow2.05, type = "b", ylim = c(0, 0.5), xlab = "Heterozygous effect (h)", ylab = "Power", main = "n= 6,000, MAF = 0.05", col = "blue", pch = 20, cex = 3, cex.lab = 3, cex.axis = 3, cex.main = 3, cex.sub = 3, lwd = 3)
points(xax, pow3.05, type = "b", pch = 15, col = "red", lwd = 3, cex = 3, lty = 2)
points(xax, pow4.05, type = "b", pch = 17, col = "black", lwd = 3, cex = 3, lty = 3)
points(xax, powcat.05, type = "b", pch = 10, col = "goldenrod", lwd = 2, cex = 3, lty = 3)
points(xax, powmax.05, type = "b", pch = 8, col = "purple", lwd = 3, cex = 3, lty = 3)
legend("topleft",
  legend = c("b=2", "b=3", "add. (b=4)", "categorical (ANOVA)", "nmax3"),
  pch = c(20, 15, 17, 10, 8), col = c("blue", "red", "black", "goldenrod", "purple"), bty = "n", cex = 4
)
abline(a = 5e-2, b = 0, lty = 2)

plot(xax, pow2.03, type = "b", ylim = c(0, 0.5), xlab = "Heterozygous effect (h)", ylab = "Power", main = "n = 30,000, MAF = 0.01", col = "blue", pch = 20, cex = 3, cex.lab = 3, cex.axis = 3, cex.main = 3, cex.sub = 3, lwd = 3)
points(xax, pow3.03, type = "b", pch = 15, col = "red", lwd = 3, cex = 3, lty = 2)
points(xax, pow4.03, type = "b", pch = 17, col = "black", lwd = 3, cex = 3, lty = 3)
points(xax, powcat.03, type = "b", pch = 10, col = "goldenrod", lwd = 3, cex = 3, lty = 3)
points(xax, powmax.03, type = "b", pch = 8, col = "purple", lwd = 3, cex = 3, lty = 3)
legend("topleft",
  legend = c("b=2", "b=3", "add. (b=4)", "categorical (ANOVA)", "nmax3"),
  pch = c(20, 15, 17, 10, 8), col = c("blue", "red", "black", "goldenrod", "purple"), bty = "n", cex = 4
)
abline(a = 5e-2, b = 0, lty = 2)



dev.off()
