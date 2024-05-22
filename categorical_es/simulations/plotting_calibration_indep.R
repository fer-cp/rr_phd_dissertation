# rm(list=ls())
require(ggplot2)


## Just edit the following line and run all:
# finame<-"Calibr_data"; R<-1E4
finame <- "Calib_M1E4"
R <- 1E4

(E <- read.table(paste(finame, ".txt", sep = "")))
colnames(E) <- c("Nominal", "Estimated", "Method")
ii <- which(E$Method == "DC")
(D <- E[ii, 1:2])
k <- nrow(D)
#
ii <- which(E$Method == "CS")
(CS <- E[ii, 1:2])
k == nrow(CS) # TRUE
#
ii <- which(E$Method == "CSP")
(CSP <- E[ii, 1:2])
k == nrow(CSP) # TRUE
#
ii <- which(E$Method == "USP")
(USP <- E[ii, 1:2])
k == nrow(USP) # TRUE
#
ii <- which(E$Method == "F")
(FET <- E[ii, 1:2])
k == nrow(FET) # TRUE
#
ii <- which(E$Method == "G")
(G <- E[ii, 1:2])
k == nrow(G) # TRUE

## Interval limits with the "exact" binomial ##
J <- matrix(NA, nrow = k, ncol = 2)
colnames(J) <- c("L", "U")
b <- .05
for (s in 1:k) {
  h <- D[s, 1]
  J[s, ] <- qbinom(c(b / 2, 1 - b / 2), R, h) / R
}
# all.equal(apply(J,1,mean),D[,1]) # Not always centered


#### Plots with the "exact" intervals (J) ####
DD <- rbind(rep(.005, 2), D, rep(.105, 2))
kk <- k + 2
JJ <- matrix(NA, nrow = kk, ncol = 2)
colnames(JJ) <- c("L", "U")
b <- .05
for (s in 1:kk) {
  h <- DD[s, 1]
  JJ[s, ] <- qbinom(c(b / 2, 1 - b / 2), R, h) / R
}

# With ggplot(2)
pdf(paste(finame, "_ggplot2.pdf", sep = ""), width = 7, height = 7 / 1.54)
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
  geom_point(color = 4, size = 2, data = CS, shape = 8) + # shape=19
  geom_point(color = 2, size = 2, data = CSP, shape = 1) +
  geom_point(color = "black", size = 2, data = USP, shape = 1) +
  geom_point(color = 3, size = 2, data = FET, shape = 1) +
  geom_point(color = 7, size = 2, data = D) +
  geom_point(color = 6, size = 2, data = G) +
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
