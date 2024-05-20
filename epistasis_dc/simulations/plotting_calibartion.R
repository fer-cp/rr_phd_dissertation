require(ggplot2)

## Just edit the following line (comment/uncomment) and run the script:
finame <- "Calibr_indep_alpha_BOOST" # Figure 2a
# finame<-"Calibr_rexp_e10_alpha_BOOST"  # Figure 2b


R <- 1E3

woboost <- function(x) { # Manager of filenames (removes suffix "_BOOST")
  x <- unlist(strsplit(x, "_"))
  x <- x[x != "BOOST"]
  return(paste(x, collapse = "_"))
}

( E <- read.table(paste(finame, ".dat", sep = "")) )
D <- E[, 5:6]
colnames(D) <- c("Nominal", "Estimated")
D
k <- nrow(D)
#
(E <- read.table(paste(woboost(finame), ".dat", sep = "")))
DC <- E[, 5:6]
colnames(DC) <- c("Nominal", "Estimated")
DC
k == nrow(DC) # TRUE

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
  geom_point(color = "red", size = 2, data = D) + # shape=19
  geom_point(color = "blue", size = 2, data = DC) +
  scale_y_continuous(
    breaks = c(0, .05, .1),
    minor_breaks = c(1, 2, 3, 4, 6, 7, 8, 9, 11, 12) / 100
  ) +
  scale_x_continuous(
    breaks = c(.01, .05, .1),
    minor_breaks = c(1, 2, 3, 4, 6, 7, 8, 9) / 100
  )
dev.off()

Inside <- J[, 1] <= D[, 2] & D[, 2] <= J[, 2]
cbind(D, J, Inside)
