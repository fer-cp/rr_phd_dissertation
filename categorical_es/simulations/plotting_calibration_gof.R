# rm(list=ls())
require(ggplot2)
R <- 1E4
alf <- c(1, 2, 5, 10) / 100

## Just edit the following line and run all:
# finame<-"Calibr_data"; R<-1E4
finame <- "alpha_CS_2allele"
E <- read.table(paste("pv_gof/", finame, ".txt", sep = ""))
E <- as.matrix(E)
colnames(E) <- NULL
E
CSbi03 <- data.frame(Nominal = alf, Estimated = E[2, ])
CSbi05 <- data.frame(Nominal = alf, Estimated = E[3, ])

finame <- "alpha_ED_2allele"
G <- read.table(paste("pv_gof/", finame, ".txt", sep = ""))
G <- as.matrix(G)
colnames(G) <- NULL
G
EDbi03 <- data.frame(Nominal = alf, Estimated = G[2, ])
EDbi05 <- data.frame(Nominal = alf, Estimated = G[3, ])

finame <- "alpha_CS_3allele"
E <- read.table(paste("pv_gof/", finame, ".txt", sep = ""))
E <- as.matrix(E)
colnames(E) <- NULL
E <- drop(E)
(CStri07 <- data.frame(Nominal = alf, Estimated = E))

finame <- "alpha_ED_3allele"
G <- read.table(paste("pv_gof/", finame, ".txt", sep = ""))
G <- as.matrix(G)
colnames(G) <- NULL
G <- drop(G)
(EDtri07 <- data.frame(Nominal = alf, Estimated = G))

finame <- "alpha_CS_3allele_bis"
E <- read.table(paste("pv_gof/", finame, ".txt", sep = ""))
E <- as.matrix(E)
colnames(E) <- NULL
E <- drop(E)
(CStri03 <- data.frame(Nominal = alf, Estimated = E))

finame <- "alpha_ED_3allele_bis"
G <- read.table(paste("pv_gof/", finame, ".txt", sep = ""))
G <- as.matrix(G)
colnames(G) <- NULL
G <- drop(G)
(EDtri03 <- data.frame(Nominal = alf, Estimated = G))


## Interval limits with the "exact" binomial ##
k <- length(alf)
J <- matrix(NA, nrow = k, ncol = 2)
colnames(J) <- c("L", "U")
b <- .05
for (s in 1:k) {
  h <- alf[s]
  J[s, ] <- qbinom(c(b / 2, 1 - b / 2), R, h) / R
}
# all.equal(apply(J,1,mean),D[,1]) # Not always centered


#### Plots with the "exact" intervals (J) ####
DD <- rbind(rep(.005, 2), EDbi03, rep(.105, 2))
alff <- c(.005, alf, .105)
kk <- k + 2
JJ <- matrix(NA, nrow = kk, ncol = 2)
colnames(JJ) <- c("L", "U")
b <- .05
for (s in 1:kk) {
  h <- alff[s]
  JJ[s, ] <- qbinom(c(b / 2, 1 - b / 2), R, h) / R
}
JJ <- as.data.frame(JJ)

# With ggplot(2)
pdf(paste("alpha_2allele_th03.pdf", sep = ""), width = 7, height = 7 / 1.11)
ggplot(DD, aes(x = Nominal, y = Estimated)) +
  # coord_cartesian(ylim = c(0, .13)) +
  geom_ribbon(aes(ymin = JJ[, 1], ymax = JJ[, 2]), alpha = .2) + # alpha=transp
  geom_segment(
    mapping = aes(
      x = .005, y = .005, xend = .105,
      yend = .105, color = "black"
    ),
    linewidth = 1, color = "black"
  ) +
  geom_point(color = 4, size = 4, data = CSbi03, shape = 19) +
  geom_point(color = 2, size = 4, data = EDbi03, shape = c(19, 1, 19, 19)) +
  # ED: red; CS: blue
  scale_y_continuous(
    breaks = c(0, .05, .1),
    minor_breaks = c(1, 2, 3, 4, 6, 7, 8, 9, 11, 12) / 100
  ) +
  scale_x_continuous(
    breaks = c(.01, .05, .1),
    minor_breaks = c(1, 2, 3, 4, 6, 7, 8, 9) / 100
  )
dev.off()

# Inside<-J[,1]<=D[,2]&D[,2]<=J[,2]
# cbind(D,J,Inside)


DD <- rbind(rep(.005, 2), EDbi05, rep(.105, 2))
pdf(paste("alpha_2allele_th05.pdf", sep = ""), width = 7, height = 7 / 1.11)
ggplot(DD, aes(x = Nominal, y = Estimated)) +
  # coord_cartesian(ylim = c(0, .13)) +
  geom_ribbon(aes(ymin = JJ[, 1], ymax = JJ[, 2]), alpha = .2) + # alpha=transp
  geom_segment(
    mapping = aes(
      x = .005, y = .005, xend = .105,
      yend = .105, color = "black"
    ),
    linewidth = 1, color = "black"
  ) +
  geom_point(color = 4, size = 4, data = CSbi05, shape = 19) +
  geom_point(color = 2, size = 4, data = EDbi05, shape = c(1, 19, 19, 19)) +
  # ED: red; CS: blue
  scale_y_continuous(
    breaks = c(0, .05, .1),
    minor_breaks = c(1, 2, 3, 4, 6, 7, 8, 9, 11, 12) / 100
  ) +
  scale_x_continuous(
    breaks = c(.01, .05, .1),
    minor_breaks = c(1, 2, 3, 4, 6, 7, 8, 9) / 100
  )
dev.off()



DD <- rbind(rep(.005, 2), EDtri07, rep(.105, 2))
pdf(paste("alpha_3allele_th07.pdf", sep = ""), width = 7, height = 7 / 1.11)
ggplot(DD, aes(x = Nominal, y = Estimated)) +
  # coord_cartesian(ylim = c(0, .13)) +
  geom_ribbon(aes(ymin = JJ[, 1], ymax = JJ[, 2]), alpha = .2) + # alpha=transp
  geom_segment(
    mapping = aes(
      x = .005, y = .005, xend = .105,
      yend = .105, color = "black"
    ),
    linewidth = 1, color = "black"
  ) +
  geom_point(color = 4, size = 4, data = CStri07, shape = 19) +
  geom_point(color = 2, size = 4, data = EDtri07, shape = c(19, 19, 19, 19)) +
  # ED: red; CS: blue
  scale_y_continuous(
    breaks = c(0, .05, .1),
    minor_breaks = c(1, 2, 3, 4, 6, 7, 8, 9, 11, 12) / 100
  ) +
  scale_x_continuous(
    breaks = c(.01, .05, .1),
    minor_breaks = c(1, 2, 3, 4, 6, 7, 8, 9) / 100
  )
dev.off()




DD <- rbind(rep(.005, 2), EDtri03, rep(.105, 2))
pdf(paste("alpha_3allele_th03.pdf", sep = ""), width = 7, height = 7 / 1.11)
ggplot(DD, aes(x = Nominal, y = Estimated)) +
  # coord_cartesian(ylim = c(0, .13)) +
  geom_ribbon(aes(ymin = JJ[, 1], ymax = JJ[, 2]), alpha = .2) + # alpha=transp
  geom_segment(
    mapping = aes(
      x = .005, y = .005, xend = .105,
      yend = .105, color = "black"
    ),
    linewidth = 1, color = "black"
  ) +
  geom_point(color = 4, size = 4, data = CStri03, shape = 19) +
  geom_point(color = 2, size = 4, data = EDtri03, shape = c(19, 19, 19, 19)) +
  # ED: red; CS: blue
  scale_y_continuous(
    breaks = c(0, .05, .1),
    minor_breaks = c(1, 2, 3, 4, 6, 7, 8, 9, 11, 12) / 100
  ) +
  scale_x_continuous(
    breaks = c(.01, .05, .1),
    minor_breaks = c(1, 2, 3, 4, 6, 7, 8, 9) / 100
  )
dev.off()
