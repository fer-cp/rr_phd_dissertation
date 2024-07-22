rm(list = ls())
source("bs_test_functions.R")
source("test_functions_ct_dcov.R")

GenerateData <- function(n, p) {
  r <- dim(p)[1]
  c <- dim(p)[2]
  data <- rmultinom(n, 1, p)
  loc <- colSums(data * (1:(r * c)))
  x <- loc %% r
  x[x == 0] <- r
  y <- (loc + 1) %/% r
  xydata <- cbind(x, y)

  dataagr <- rowSums(data)
  freq <- matrix(dataagr, ncol = c)

  output <- list()
  output[[1]] <- xydata
  output[[2]] <- freq
  return(output)
}


#### Decaying marginals example ####

# M<-700
M <- 10000

n <- 100
B <- 999
r <- 5
c <- 8 # Sometimes 6 x 6

px <- (1 / 2)^(1:r)
py <- (1 / 2)^(1:c)
Eps <- (0:10) * 0.012 # good for n=100
# Note that eps cannot be higher than any of these two:
(ub1 <- 1 / (8 * (1 - 2^(-r)) * (1 - 2^(-c)))) # 0.1295383
(ub2 <- 1 - 1 / (4 * (1 - 2^(-r)) * (1 - 2^(-c)))) # 0.7409235
min(ub1, ub2)

perturb <- matrix(rep(0, r * c), nrow = r)
perturb[1, 1] <- 1
perturb[1, 2] <- -1
perturb[2, 1] <- -1
perturb[2, 2] <- 1

DCpower <- rep(NA, length(Eps))
chisqpower <- rep(NA, length(Eps))
Gpower <- rep(NA, length(Eps))
USPpower <- rep(NA, length(Eps))
# DCpermpower=rep(NA,length(Eps))
chisqpermpower <- rep(NA, length(Eps))
# Gpermpower=rep(NA,length(Eps))
Fisherpower <- rep(NA, length(Eps))

# DCStats=list()
# USPStats=list()

for (i in 1:length(Eps)) {
  # cat("i=",i,"out of",length(Eps))
  eps <- Eps[i]
  p <- px %*% t(py)
  p <- p / sum(p)
  p <- p + eps * perturb

  DCpval <- rep(NA, M)
  chisqpval <- rep(NA, M)
  Gpval <- rep(NA, M)
  USPpval <- rep(NA, M)
  # DCpermpval=rep(NA,M)
  chisqpermpval <- rep(NA, M)
  # Gpermpval=rep(NA,M)
  Fisherpval <- rep(NA, M)

  # DCstat=rep(NA,M)
  # USPstat=rep(NA,M)

  set.seed(1)
  for (m in 1:M) {
    cat("i=", i, "out of", length(Eps), ". m=", m, "out of", M, "\n")
    data <- GenerateData(n, p)
    freq <- data[[2]]

    # USPtest=USPDiscrete(freq,B)
    USPpval[m] <- USPDiscrete(freq, B)$p.value
    DCpval[m] <- test.dcov(freq)
    # DCpermpval[m]=test.dcov.perm(freq,B=999)

    # DCstat[m]=statis.dcov(freq)
    # freqc=colSums(freq); freqr=rowSums(freq); n=sum(freq)
    # sumr2=sum(freqr^2); sumc2=sum(freqc^2)
    # USPstat[m]=USPtest$TestStat+
    #   (sumr2+sumc2)/(n*(n-1)*(n-3))+
    #   (3*n-2)*sumr2*sumc2/(n^3*(n-1)*(n-2)*(n-3))-n/((n-1)*(n-3))

    chisq <- ChisqStat(freq)
    chisqpval[m] <- pchisq(chisq, df = (r - 1) * (c - 1), lower.tail = F)
    # chisqpval[m]=test.chisq(freq) # Mine generates NA's for empty col/rows

    chisqpermpval[m] <- ChisqPerm(freq, B)

    # gstat=GStat(freq)
    # Gpval[m]=pchisq(gstat, df=(r-1)*(c-1), lower.tail=F)
    Gpval[m] <- test.g(freq)

    # Gpermpval[m]=GPerm(freq,B)

    Fisherpval[m] <- fisher.test(freq, simulate.p.value = T, B = 999)$p.value
  }
  # USPStats[[i]]=USPstat
  # DCStats[[i]]=DCstat

  (DCpower[i] <- sum(DCpval <= 0.05) / M)
  (chisqpower[i] <- sum(chisqpval <= 0.05) / M)
  (Gpower[i] <- sum(Gpval <= 0.05) / M)
  (USPpower[i] <- sum(USPpval <= 0.05) / M)
  # (DCpermpower[i]=sum(DCpermpval<=0.05)/M)
  (chisqpermpower[i] <- sum(chisqpermpval <= 0.05) / M)
  # (Gpermpower[i]=sum(Gpermpval<=0.05)/M)
  (Fisherpower[i] <- sum(Fisherpval <= 0.05) / M)

  if (i == 1) {
    write.table(DCpval, "pvals/DCpval_H0_dec_marg_M1E4.txt", row.names = F, col.names = F)
    write.table(USPpval, "pvals/USPpval_H0_dec_marg_M1E4.txt", row.names = F, col.names = F)
    write.table(chisqpval, "pvals/CSpval_H0_dec_marg_M1E4.txt", row.names = F, col.names = F)
    write.table(chisqpermpval, "pvals/CSPpval_H0_dec_marg_M1E4.txt", row.names = F, col.names = F)
    write.table(Gpval, "pvals/Gpval_H0_dec_marg_M1E4.txt", row.names = F, col.names = F)
    write.table(Fisherpval, "pvals/Fpval_H0_dec_marg_M1E4.txt", row.names = F, col.names = F)
  }
}

SEb <- sqrt(USPpower * (1 - USPpower) / M)
SEc <- sqrt(chisqpower * (1 - chisqpower) / M)
SEcp <- sqrt(chisqpermpower * (1 - chisqpermpower) / M)
SEg <- sqrt(Gpower * (1 - Gpower) / M)
# SEgp=sqrt(Gpermpower*(1-Gpermpower)/M)
SEd <- sqrt(DCpower * (1 - DCpower) / M)
# SEdp=sqrt(DCpermpower*(1-DCpermpower)/M)
SEf <- sqrt(Fisherpower * (1 - Fisherpower) / M)

# Color codes:
# DC: col=7 ~ "goldenrod"
# chisq: col=4 ~ "lightskyblue"
# chisqperm: col=2 ~ "firebrick" ("darkred")
# USP: col=\emptyset="black"
# Fisher: col=3 ~ "limegreen"
# G: col=6 ~ "darkorchid" (purple - magenta)


plot(Eps, USPpower, type = "l", ylim = c(0, 1), xlab = expression(epsilon), ylab = "Power", las = 1, cex.lab = 1.4)
lines(Eps, DCpower, col = 7)
lines(Eps, chisqpower, col = 4)
lines(Eps, Gpower, col = 6)

lines(Eps, chisqpermpower, col = 2)
# lines(Eps,DCpermpower,col=8)
# lines(Eps,Gpermpower,col=3)
lines(Eps, Fisherpower, col = 3)


for (i in 1:length(USPpower)) {
  points(Eps[i], USPpower[i], cex = 0.75)
  segments(Eps[i], USPpower[i] - 3 * SEb[i], Eps[i], USPpower[i] + 3 * SEb[i])

  points(Eps[i], chisqpower[i], cex = 0.75, col = 4)
  segments(Eps[i], chisqpower[i] - 3 * SEc[i], Eps[i], chisqpower[i] + 3 * SEc[i], col = 4)

  points(Eps[i], chisqpermpower[i], cex = 0.75, col = 2)
  segments(Eps[i], chisqpermpower[i] - 3 * SEcp[i], Eps[i], chisqpermpower[i] + 3 * SEcp[i], col = 2)

  points(Eps[i], DCpower[i], cex = 0.75, col = 7)
  segments(Eps[i], DCpower[i] - 3 * SEd[i], Eps[i], DCpower[i] + 3 * SEd[i], col = 7)

  # points(Eps[i],DCpermpower[i],cex=0.75,col=8)
  # segments(Eps[i],DCpermpower[i]-3*SEdp[i],Eps[i],DCpermpower[i]+3*SEdp[i],col=8)

  points(Eps[i], Gpower[i], cex = 0.75, col = 6)
  segments(Eps[i], Gpower[i] - 3 * SEg[i], Eps[i], Gpower[i] + 3 * SEg[i], col = 6)

  # points(Eps[i],Gpermpower[i],cex=0.75,col=3)
  # segments(Eps[i],Gpermpower[i]-3*SEgp[i],Eps[i],Gpermpower[i]+3*SEgp[i],col=3)

  points(Eps[i], Fisherpower[i], cex = 0.75, col = 3)
  segments(Eps[i], Fisherpower[i] - 3 * SEf[i], Eps[i], Fisherpower[i] + 3 * SEf[i], col = 3)
}
# With n=100; M=10000; B=999; r=5; c=8; Eps=0.012*(0:10)

write.table(DCpower, "pvals/DCpower_dec_marg_M1E4.txt", row.names = F, col.names = F)
write.table(USPpower, "pvals/USPpower_dec_marg_M1E4.txt", row.names = F, col.names = F)
write.table(chisqpower, "pvals/CSpower_dec_marg_M1E4.txt", row.names = F, col.names = F)
write.table(chisqpermpower, "pvals/CSPpower_dec_marg_M1E4.txt", row.names = F, col.names = F)
write.table(Gpower, "pvals/Gpower_dec_marg_M1E4.txt", row.names = F, col.names = F)
write.table(Fisherpower, "pvals/Fpower_dec_marg_M1E4.txt", row.names = F, col.names = F)

# Exporting for the calibration plot
Calib <- matrix(NA, 24, 3)
avals <- c(.01, .02, .05, .1)
methods <- c("CS", "CSP", "DC", "F", "G", "USP")
for (i in 1:6) {
  m <- methods[i]
  finame <- paste("pvals/", m, "pval_H0_dec_marg_M1E4.txt", sep = "")
  for (j in 1:4) {
    a <- avals[j]
    pv <- drop(read.table(finame))
    Calib[4 * (i - 1) + j, ] <- c(a, mean(pv < a), m)
  }
}
Calib
write.table(Calib, "Calib_M1E4.txt", row.names = F, col.names = F, quote = F)
