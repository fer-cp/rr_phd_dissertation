rm(list = ls())

source("test_functions_ct_dcov.R")

# Adult:
M <- 1E4 # 10
n <- 500 # 50

# # Toy:
# M<-10
# n<-50


#### Calibration biallelic ####

alf <- c(1, 2, 5, 10) / 100
la <- length(alf)

A_ED <- matrix(NA, 3, 4)
A_CS <- matrix(NA, 3, 4)

tt <- c(1 / 10, 1 / 3, 1 / 2)
lt <- length(tt)

for (i in 1:lt) {
  t <- tt[i]
  cat("Taking theta =", t, "\n")
  p0 <- c(t^2, 2 * t * (1 - t), (1 - t)^2)
  p <- p0 # Calibration
  pv_ED <- rep(NA, M)
  pv_CS <- rep(NA, M)
  set.seed(1) # Determ
  for (m in 1:M) {
    x <- sample(1:3, n, replace = T, prob = p)
    o <- table(c(1:3, x)) - 1
    pv_ED[m] <- test.e(o, p0)
    pv_CS[m] <- test.chisq.gof(o, p0)
  }
  rt <- as.character(100 * round(t, 2))
  # fi<-paste("pv_gof/pv_ED_H0_2allele_maf_0_",rt,".txt")
  # write.table(pv_ED,fi,row.names=F,col.names=F)
  # fi<-paste("pv_gof/pv_CS_H0_2allele_maf_0_",rt,".txt")
  # write.table(pv_CS,fi,row.names=F,col.names=F)
  for (j in 1:la) {
    a <- alf[j]
    A_ED[i, j] <- mean(pv_ED < a)
    A_CS[i, j] <- mean(pv_CS < a)
  }
}

fi <- paste("pv_gof/alpha_ED_2allele.txt")
write.table(A_ED, fi, row.names = F, col.names = F)
fi <- paste("pv_gof/alpha_CS_2allele.txt")
write.table(A_CS, fi, row.names = F, col.names = F)




#### Power biallelic: Model R ####
# r=1 is H0, but I repeat it to check that everything works

rgrid <- c(1:10, 10 * (2:7), 81)
nr <- length(rgrid)

B_ED <- rep(NA, nr)
B_CS <- rep(NA, nr)

p0 <- c(81, 18, 1) / 100

for (i in 1:nr) {
  r <- rgrid[i]
  cat("Taking r =", r, "\n")
  p <- c(
    .81 + (1 - r) / 100,
    .18,
    r / 100
  )
  pv_ED <- rep(NA, M)
  pv_CS <- rep(NA, M)
  set.seed(1) # Determ
  for (m in 1:M) {
    x <- sample(1:3, n, replace = T, prob = p)
    o <- table(c(1:3, x)) - 1
    pv_ED[m] <- test.e(o, p0)
    pv_CS[m] <- test.chisq.gof(o, p0)
  }
  B_ED[i] <- mean(pv_ED < .05)
  B_CS[i] <- mean(pv_CS < .05)
}

fi <- paste("pv_gof/power_ED_2allele_modelR.txt")
write.table(cbind(rgrid, B_ED), fi, row.names = F, col.names = F)
fi <- paste("pv_gof/power_CS_2allele_modelR.txt")
write.table(cbind(rgrid, B_CS), fi, row.names = F, col.names = F)



#### Power biallelic: Model S ####
# s=0 is H0, but I repeat it to check that everything works

sgrid <- c((0:9) / 100, (1:10) / 10)
ns <- length(sgrid)

B_ED <- rep(NA, ns)
B_CS <- rep(NA, ns)

p0 <- c(4, 4, 1) / 9

for (i in 1:ns) {
  s <- sgrid[i]
  cat("Taking s =", s, "\n")
  p <- p0
  p[-3] <- p[-3] * (1 - s)
  p[3] <- (1 + 8 * s) / 9

  pv_ED <- rep(NA, M)
  pv_CS <- rep(NA, M)
  set.seed(1) # Determ
  for (m in 1:M) {
    x <- sample(1:3, n, replace = T, prob = p)
    o <- table(c(1:3, x)) - 1
    pv_ED[m] <- test.e(o, p0)
    pv_CS[m] <- test.chisq.gof(o, p0)
  }
  B_ED[i] <- mean(pv_ED < .05)
  B_CS[i] <- mean(pv_CS < .05)
}

fi <- paste("pv_gof/power_ED_2allele_modelS.txt")
write.table(cbind(sgrid, B_ED), fi, row.names = F, col.names = F)
fi <- paste("pv_gof/power_CS_2allele_modelS.txt")
write.table(cbind(sgrid, B_CS), fi, row.names = F, col.names = F)


#### Power biallelic: Model K ####
# k=0 is H0, but I repeat it to check that everything works

kgrid <- c((-5:-1) / 5, (-19:19) / 100, (1:5) / 5)
if (n == 100) {
  kgrid <- c((-5:-2) / 5, (-39:39) / 100, (2:5) / 5)
}
nk <- length(kgrid)

B_ED <- rep(NA, nk)
B_CS <- rep(NA, nk)

p0 <- c(1, 2, 1) / 4

for (i in 1:nk) {
  k <- kgrid[i]
  cat("Taking k =", k, "\n")
  p <- p0
  p[-2] <- p[-2] - k / 4
  p[2] <- p[2] + k / 2

  pv_ED <- rep(NA, M)
  pv_CS <- rep(NA, M)
  set.seed(1) # Determ
  for (m in 1:M) {
    x <- sample(1:3, n, replace = T, prob = p)
    o <- table(c(1:3, x)) - 1
    pv_ED[m] <- test.e(o, p0)
    pv_CS[m] <- test.chisq.gof(o, p0)
  }
  B_ED[i] <- mean(pv_ED < .05)
  B_CS[i] <- mean(pv_CS < .05)
}

fi <- paste("pv_gof/power_ED_2allele_modelK.txt")
write.table(cbind(kgrid, B_ED), fi, row.names = F, col.names = F)
fi <- paste("pv_gof/power_CS_2allele_modelK.txt")
write.table(cbind(kgrid, B_CS), fi, row.names = F, col.names = F)


#### Calibration triallelic ####

alf <- c(1, 2, 5, 10) / 100
la <- length(alf)

A_ED <- rep(NA, la)
A_CS <- rep(NA, la)

cat("Doing the triallelic calibration\n")
t1 <- .7
t2 <- .25
t3 <- .05
p0 <- c(t1^2, t2^2, t3^2, 2 * t1 * t2, 2 * t1 * t3, 2 * t2 * t3)
p0bis <- c(490, 62.5, 2.5, 350, 70, 25) / 1000
if (!all(abs(p0 - p0bis) < 1e-7)) cat("Problems...\n")
p <- p0 # Calibration
pv_ED <- rep(NA, M)
pv_CS <- rep(NA, M)
set.seed(1) # Determ
for (m in 1:M) {
  x <- sample(1:6, n, replace = T, prob = p)
  o <- table(c(1:6, x)) - 1
  pv_ED[m] <- test.e(o, p0)
  pv_CS[m] <- test.chisq.gof(o, p0)
}
# fi<-paste("pv_gof/pv_ED_H0_3allele.txt")
# write.table(pv_ED,fi,row.names=F,col.names=F)
# fi<-paste("pv_gof/pv_CS_H0_3allele.txt")
# write.table(pv_CS,fi,row.names=F,col.names=F)
for (j in 1:la) {
  a <- alf[j]
  A_ED[j] <- mean(pv_ED < a)
  A_CS[j] <- mean(pv_CS < a)
}


fi <- paste("pv_gof/alpha_ED_3allele.txt")
write.table(A_ED, fi, row.names = F, col.names = F)
fi <- paste("pv_gof/alpha_CS_3allele.txt")
write.table(A_CS, fi, row.names = F, col.names = F)


#### Calibration triallelic (bis) ####

alf <- c(1, 2, 5, 10) / 100
la <- length(alf)

A_ED <- rep(NA, la)
A_CS <- rep(NA, la)

cat("Doing the triallelic calibration\n")
t1 <- 1 / 3
t2 <- 1 / 3
t3 <- 1 / 3
p0 <- c(t1^2, t2^2, t3^2, 2 * t1 * t2, 2 * t1 * t3, 2 * t2 * t3)
p0bis <- c(rep(1, 3), rep(2, 3)) / 9
if (!all(abs(p0 - p0bis) < 1e-7)) cat("Problems...\n")
p <- p0 # Calibration
pv_ED <- rep(NA, M)
pv_CS <- rep(NA, M)
set.seed(1) # Determ
for (m in 1:M) {
  x <- sample(1:6, n, replace = T, prob = p)
  o <- table(c(1:6, x)) - 1
  pv_ED[m] <- test.e(o, p0)
  pv_CS[m] <- test.chisq.gof(o, p0)
}
# fi<-paste("pv_gof/pv_ED_H0_3allele_bis.txt")
# write.table(pv_ED,fi,row.names=F,col.names=F)
# fi<-paste("pv_gof/pv_CS_H0_3allele_bis.txt")
# write.table(pv_CS,fi,row.names=F,col.names=F)
for (j in 1:la) {
  a <- alf[j]
  A_ED[j] <- mean(pv_ED < a)
  A_CS[j] <- mean(pv_CS < a)
}


fi <- paste("pv_gof/alpha_ED_3allele_bis.txt")
write.table(A_ED, fi, row.names = F, col.names = F)
fi <- paste("pv_gof/alpha_CS_3allele_bis.txt")
write.table(A_CS, fi, row.names = F, col.names = F)



#### Power triallelic: Model R ####
# r=1 is H0, but I repeat it to check that everything works

rgrid <- c(1:8, 11, 14, 17.81818)
nr <- length(rgrid)

B_ED <- rep(NA, nr)
B_CS <- rep(NA, nr)

t1 <- .7
t2 <- .25
t3 <- .05
p0 <- c(t1^2, t2^2, t3^2, 2 * t1 * t2, 2 * t1 * t3, 2 * t2 * t3)
p0bis <- c(490, 62.5, 2.5, 350, 70, 25) / 1000
if (!all(abs(p0 - p0bis) < 1e-7)) cat("Problems...\n")

for (i in 1:nr) {
  r <- rgrid[i]
  cat("Taking r =", r, "\n")
  p <- p0
  p[c(3, 6)] <- r * p[c(3, 6)]
  p[1] <- p[1] - .0275 * (1 - r)
  pv_ED <- rep(NA, M)
  pv_CS <- rep(NA, M)
  set.seed(1) # Determ
  for (m in 1:M) {
    x <- sample(1:6, n, replace = T, prob = p)
    o <- table(c(1:6, x)) - 1
    pv_ED[m] <- test.e(o, p0)
    pv_CS[m] <- test.chisq.gof(o, p0)
  }
  B_ED[i] <- mean(pv_ED < .05)
  B_CS[i] <- mean(pv_CS < .05)
}

fi <- paste("pv_gof/power_ED_3allele_modelR.txt")
write.table(cbind(rgrid, B_ED), fi, row.names = F, col.names = F)
fi <- paste("pv_gof/power_CS_3allele_modelR.txt")
write.table(cbind(rgrid, B_CS), fi, row.names = F, col.names = F)


#### Power triallelic: Model S ####
# s=0 is H0, but I repeat it to check that everything works

sgrid <- c((0:10) / 100, .2, .4, .7, 1)
ns <- length(sgrid)

B_ED <- rep(NA, ns)
B_CS <- rep(NA, ns)

t1 <- .7
t2 <- .25
t3 <- .05
p0 <- c(t1^2, t2^2, t3^2, 2 * t1 * t2, 2 * t1 * t3, 2 * t2 * t3)

for (i in 1:ns) {
  s <- sgrid[i]
  cat("Taking s =", s, "\n")
  p <- p0
  p[-2] <- p[-2] * (1 - s)
  p[2] <- (1 + 15 * s) / 16

  pv_ED <- rep(NA, M)
  pv_CS <- rep(NA, M)
  set.seed(1) # Determ
  for (m in 1:M) {
    x <- sample(1:6, n, replace = T, prob = p)
    o <- table(c(1:6, x)) - 1
    pv_ED[m] <- test.e(o, p0)
    pv_CS[m] <- test.chisq.gof(o, p0)
  }
  B_ED[i] <- mean(pv_ED < .05)
  B_CS[i] <- mean(pv_CS < .05)
}

fi <- paste("pv_gof/power_ED_3allele_modelS.txt")
write.table(cbind(sgrid, B_ED), fi, row.names = F, col.names = F)
fi <- paste("pv_gof/power_CS_3allele_modelS.txt")
write.table(cbind(sgrid, B_CS), fi, row.names = F, col.names = F)


#### Power triallelic: Model K ####
# k=0 is H0, but I repeat it to check that everything works

kgrid <- c((0:20) / 100, .4, .7, 1)
nk <- length(kgrid)

B_ED <- rep(NA, nk)
B_CS <- rep(NA, nk)

t1 <- 1 / 3
t2 <- 1 / 3
t3 <- 1 / 3
p0 <- c(t1^2, t2^2, t3^2, 2 * t1 * t2, 2 * t1 * t3, 2 * t2 * t3)

for (i in 1:nk) {
  k <- kgrid[i]
  cat("Taking k =", k, "\n")
  p <- p0
  p[1:3] <- (2 * k + 1) / 9
  p[4:6] <- (2 - 2 * k) / 9
  pv_ED <- rep(NA, M)
  pv_CS <- rep(NA, M)
  set.seed(1) # Determ
  for (m in 1:M) {
    x <- sample(1:6, n, replace = T, prob = p)
    o <- table(c(1:6, x)) - 1
    pv_ED[m] <- test.e(o, p0)
    pv_CS[m] <- test.chisq.gof(o, p0)
  }
  B_ED[i] <- mean(pv_ED < .05)
  B_CS[i] <- mean(pv_CS < .05)
}

fi <- paste("pv_gof/power_ED_3allele_modelK.txt")
write.table(cbind(kgrid, B_ED), fi, row.names = F, col.names = F)
fi <- paste("pv_gof/power_CS_3allele_modelK.txt")
write.table(cbind(kgrid, B_CS), fi, row.names = F, col.names = F)
