M <- 1E4

#### Plot for model 2S ####
finame <- "power_CS_2allele_modelS"
CS <- read.table(paste("pv_gof/", finame, ".txt", sep = ""))
(CS <- as.matrix(CS))
x <- CS[, 1]
y <- CS[, 2]

finame <- "power_ED_2allele_modelS"
ED <- read.table(paste("pv_gof/", finame, ".txt", sep = ""))
(ED <- as.matrix(ED))
xx <- ED[, 1]
z <- ED[, 2]

all.equal(x, xx) # TRUE

x <- x[1:11]
x <- c(x, (3:6) / 20)
y <- y[1:11]
y <- c(y, rep(1, 4))
z <- z[1:11]
z <- c(z, rep(1, 4))

SEy <- sqrt(y * (1 - y) / M)
SEz <- sqrt(z * (1 - z) / M)

pdf(paste("power_2s.pdf", sep = ""), width = 7, height = 7 / 1.3)
plot(x, y,
  type = "l", ylim = c(0, 1), xlab = "s",
  ylab = "Power", las = 1, cex.lab = 1.4, col = 4
)
lines(x, z, col = 2)

for (i in 1:length(y)) {
  points(x[i], y[i], cex = 0.75, col = 4, pch = 19)
  segments(x[i], y[i] - 3 * SEy[i],
    x[i], y[i] + 3 * SEy[i],
    col = 4
  )

  points(x[i], z[i], cex = 0.75, col = 2, pch = 19)
  segments(x[i], z[i] - 3 * SEz[i],
    x[i], z[i] + 3 * SEz[i],
    col = 2
  )
}
dev.off()



#### Plot for model 2K ####
finame <- "power_CS_2allele_modelK"
CS <- read.table(paste("pv_gof/", finame, ".txt", sep = ""))
(CS <- as.matrix(CS))
x <- CS[, 1]
y <- CS[, 2]

finame <- "power_ED_2allele_modelK"
ED <- read.table(paste("pv_gof/", finame, ".txt", sep = ""))
(ED <- as.matrix(ED))
xx <- ED[, 1]
z <- ED[, 2]

all.equal(x, xx) # TRUE

x <- x[3:47]
y <- y[3:47]
z <- z[3:47]

SEy <- sqrt(y * (1 - y) / M)
SEz <- sqrt(z * (1 - z) / M)

pdf(paste("power_2k.pdf", sep = ""), width = 7, height = 7 / 1.3)
plot(x, y,
  type = "l", ylim = c(0, 1), xlab = "k",
  ylab = "Power", las = 1, cex.lab = 1.4, col = 4
)
lines(x, z, col = 2)

for (i in 1:length(y)) {
  points(x[i], y[i], cex = 0.75, col = 4, pch = 19)
  segments(x[i], y[i] - 3 * SEy[i],
    x[i], y[i] + 3 * SEy[i],
    col = 4
  )

  points(x[i], z[i], cex = 0.75, col = 2, pch = 19)
  segments(x[i], z[i] - 3 * SEz[i],
    x[i], z[i] + 3 * SEz[i],
    col = 2
  )
}
dev.off()



#### Plot for model 3S ####
finame <- "power_CS_3allele_modelS"
CS <- read.table(paste("pv_gof/", finame, ".txt", sep = ""))
(CS <- as.matrix(CS))
x <- CS[, 1]
y <- CS[, 2]

finame <- "power_ED_3allele_modelS"
ED <- read.table(paste("pv_gof/", finame, ".txt", sep = ""))
(ED <- as.matrix(ED))
xx <- ED[, 1]
z <- ED[, 2]

all.equal(x, xx) # TRUE

x <- x[1:11]
x <- c(x, (11:20) / 100)
y <- y[1:11]
y <- c(y, rep(1, 10))
z <- z[1:11]
z <- c(z, rep(1, 10))

SEy <- sqrt(y * (1 - y) / M)
SEz <- sqrt(z * (1 - z) / M)

pdf(paste("power_3s.pdf", sep = ""), width = 7, height = 7 / 1.3)
plot(x, y,
  type = "l", ylim = c(0, 1), xlab = "s",
  ylab = "Power", las = 1, cex.lab = 1.4, col = 4
)
lines(x, z, col = 2)

for (i in 1:length(y)) {
  points(x[i], y[i], cex = 0.75, col = 4, pch = 19)
  segments(x[i], y[i] - 3 * SEy[i],
    x[i], y[i] + 3 * SEy[i],
    col = 4
  )

  points(x[i], z[i], cex = 0.75, col = 2, pch = 19)
  segments(x[i], z[i] - 3 * SEz[i],
    x[i], z[i] + 3 * SEz[i],
    col = 2
  )
}
dev.off()




#### Plot for model 3K ####
finame <- "power_CS_3allele_modelK"
CS <- read.table(paste("pv_gof/", finame, ".txt", sep = ""))
(CS <- as.matrix(CS))
x <- CS[, 1]
y <- CS[, 2]

finame <- "power_ED_3allele_modelK"
ED <- read.table(paste("pv_gof/", finame, ".txt", sep = ""))
(ED <- as.matrix(ED))
xx <- ED[, 1]
z <- ED[, 2]

all.equal(x, xx) # TRUE

x <- x[1:21]
x <- c(x, .3, .4)
y <- y[1:21]
y <- c(y, rep(1, 2))
z <- z[1:21]
z <- c(z, rep(1, 2))


SEy <- sqrt(y * (1 - y) / M)
SEz <- sqrt(z * (1 - z) / M)

pdf(paste("power_3k.pdf", sep = ""), width = 7, height = 7 / 1.3)
plot(x, y,
  type = "l", ylim = c(0, 1), xlab = "k",
  ylab = "Power", las = 1, cex.lab = 1.4, col = 4
)
lines(x, z, col = 2)

for (i in 1:length(y)) {
  points(x[i], y[i], cex = 0.75, col = 4, pch = 19)
  segments(x[i], y[i] - 3 * SEy[i],
    x[i], y[i] + 3 * SEy[i],
    col = 4
  )

  points(x[i], z[i], cex = 0.75, col = 2, pch = 19)
  segments(x[i], z[i] - 3 * SEz[i],
    x[i], z[i] + 3 * SEz[i],
    col = 2
  )
}
dev.off()
