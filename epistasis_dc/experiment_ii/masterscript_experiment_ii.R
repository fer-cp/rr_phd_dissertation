# We begin with the 6456 SNPs at the "data_final" database
# which in turn was generated with the "filtering_snps_experiment_ii.R"
# using the whole GWAS data that we are not allowed to share.

# At many points of the script, we generate intermediate result files,
# in order to ease running only parts of it.

rm(list = ls())

#### Auxiliary functions ####
source("../sim_functions_snp_snp.R") # Our test
source("../aux_functions_triu.R") # Working with triangular matrices
source("aux_functions_nearness.R") # Studying physical distance

#### Reading the anonymised data ####
X <- as.matrix(read.table("Matrix_X.dat"))
Y <- as.matrix(read.table("Matrix_Y.dat"))

#### Applying our test procedure ####
tY <- test1(Y) # test2 for b=2 (unnecessary) # matrix
(p <- ncol(Y)) # 6456 x 6456 symmetric, with NA's on the diagonal
rm(Y)
write.table(tY, "pval_Y_brute.dat", col.names = F, row.names = F)
rm(tY)

tX <- test1(X) # 6456 x 6456
(p <- ncol(X))
# p==ncol(X) # TRUE
rm(X)
write.table(tX, "pval_X_brute.dat", col.names = F, row.names = F)
rm(tX)

#### FDR correction via Benjamini-Hochberg ####
tYc <- as.matrix(read.table("pval_Y_brute.dat")) # DF to matrix
tYc[lower.tri(tYc)] <- NA # triu
tYc[1:5, 1:5]
dim(tYc) # L x L
tYv <- as.numeric(tYc) # triu to vec
length(tYv) == p^2 # TRUE
rm(tYc)
tYv <- tYv[!is.na(tYv)]
summary(tYv) # without NA's
length(tYv) == (p^2 - p) / 2 # TRUE
tY.BH <- p.adjust(tYv, method = "BH")
summary(tY.BH)
length(tY.BH) == (p^2 - p) / 2 # TRUE
rm(tYv)
write.table(tY.BH, "tY_BH.dat")
rm(tY.BH)

tXc <- as.matrix(read.table("pval_X_brute.dat")) # DF to matrix
tXc[lower.tri(tXc)] <- NA # triu
tXc[1:5, 1:5]
dim(tXc) # L x L
tXv <- as.numeric(tXc) # triu to vec
length(tXv) == p^2 # TRUE
rm(tXc)
tXv <- tXv[!is.na(tXv)]
summary(tXv) # without NA's
length(tXv) == (p^2 - p) / 2 # TRUE
tX.BH <- p.adjust(tXv, method = "BH")
summary(tX.BH)
length(tX.BH) == (p^2 - p) / 2 # TRUE
rm(tXv)
write.table(tX.BH, "tX_BH.dat")
rm(tX.BH)


#### Counting (un)detected SNP pairs ####
tX.BH <- drop(as.matrix(read.table("tX_BH.dat")))
tY.BH <- drop(as.matrix(read.table("tY_BH.dat")))
# vectors of length (p^2-p)/2
alpha <- .05
ix <- tX.BH < alpha
iy <- tY.BH < alpha
table(ix, iy)
# Rows correspond to tX (cases); cols to tY (controls)
#           FALSE     TRUE
# FALSE  20833623     1326
# TRUE      1465      326

p <- 6456

indA <- which(ix & iy)
mean(ix & iy) # 0.00156 % of pairs convey ancestrality
length(indA) # 326 out of (p^2-p)/2
max(indA) / (p^2 - p) * 2 # 0.9978325

indY <- which(ix & !iy)
mean(ix & !iy) # 0.00703 % of pairs are in presumed epistasis
length(indY) # 1465 out of (p^2-p)/2
max(indY) / (p^2 - p) * 2 # 0.9996335

indX <- which(!ix & iy) # FPs (cf. Biom)
mean(!ix & iy) # 0.00636 % of pairs are in fake epistasis
length(indX) # 1326 out of (p^2-p)/2
max(indX) / (p^2 - p) * 2 # 0.9996721

indN <- which(!ix & !iy)
mean(!ix & !iy) # 99.99 % of pairs present no signal at all
length(indN) # 20 833 623 out of (p^2-p)/2
max(indN) / (p^2 - p) * 2 # 1

length(intersect(indX, indY)) # 0

prop.test( c(length(indY), length(indX)),
           alternative = "greater",  rep((p^2 - p) / 2, 2) ) # P = 0.0045

#### Translating the indices into SNP pairs ####
hitsY <- t(sapply(indY, itoij)) # Pairs by row
write.table(hitsY, "snppairs_hits_Y.dat", row.names = F, col.names = F)
#
hitsX <- t(sapply(indX, itoij)) # Pairs by row
write.table(hitsX, "snppairs_hits_X.dat", row.names = F, col.names = F)
#
hitsA <- t(sapply(indA, itoij)) # Pairs by row
write.table(hitsA, "snppairs_hits_A.dat", row.names = F, col.names = F)
#
hitsN <- t(sapply(indN, itoij)) # Pairs by row (runs slower)
write.table(hitsN, "snppairs_hits_N.dat", row.names = F, col.names = F)

#### Filtering by physical distance ####
hitsY <- as.matrix(read.table("snppairs_hits_Y.dat"))
hitsYc <- hitsY[!apply(hitsY, 1, nearby), ]
(nrhyc <- nrow(hitsYc)) # 1272  (result present in the manuscript)
hitsYc[c(1:4, (nrhyc - 2):nrhyc), ]
snpsYc <- as.numeric(as.matrix(hitsYc))
length(snpsYc) == 2 * nrhyc # TRUE
snpsYcu <- unique(snpsYc)
length(snpsYcu) # 1539  (result present in the manuscript)
(fryc <- table(table(snpsYc)))
# #Interactions involved 1    2   3   4   5   6  7  8   9  10  11
# #SNPs                 960 329 150  61  23   5  7  1   1   1   1
sum(fryc) == length(snpsYcu) # TRUE: 1539

hitsX <- as.matrix(read.table("snppairs_hits_X.dat"))
hitsXc <- hitsX[!apply(hitsX, 1, nearby), ]
(nrhxc <- nrow(hitsXc)) # 1137  (result present in the manuscript)
hitsXc[c(1:4, (nrhxc - 2):nrhxc), ]
snpsXc <- as.numeric(as.matrix(hitsXc))
length(snpsXc) == 2 * nrhxc # TRUE
snpsXcu <- unique(snpsXc)
length(snpsXcu) # 1439  (result present in the manuscript)
(frxc <- table(table(snpsXc)))
# #Interactions involved 1    2   3   4   5   6  7
# #SNPs                 912 324 128  51  19   4  1
sum(frxc) == length(snpsXcu) # TRUE: 1439

hitsA <- as.matrix(read.table("snppairs_hits_A.dat"))
hitsAc <- hitsA[!apply(hitsA, 1, nearby), ]
(nrhac <- nrow(hitsAc)) # 7

hitsN <- as.matrix(read.table("snppairs_hits_N.dat"))
hitsNc <- hitsN[!apply(hitsN, 1, nearby), ]
(nrhnc <- nrow(hitsNc)) # 20813042

rm(hitsN)
rm(hitsNc)

#### Testing for excess in pairs and unique SNPs ####
p <- 6456
(tt <- nrhyc + nrhxc + nrhac + nrhnc) # 20815458
prop.test( c(nrow(hitsYc), nrow(hitsXc)), rep(tt, 2),
          alternative = "greater" ) # 0.0032 (result present in the manuscript)
prop.test( c(length(snpsYcu), length(snpsXcu)), rep(p, 2),
           alternative = "greater" ) # 0.019 (result present in the manuscript)


#### Controlling for common SNPs between *X and *Y ####

length(intersect(snpsXcu, snpsYcu)) # 800 common SNPs

snpsYnc <- snpsYc[!is.element(snpsYc, snpsXc)]
(py <- length(unique(snpsYnc))) # 739 = 1539 - 800
(frync <- table(table(snpsYnc)))
# #Interactions involved 1   2    3   4   5   6   7   9   10
# #SNPs                 558 115  40  15   6   1   2   1   1
sum(frync) == py # TRUE # 739

snpsXnc <- snpsXc[!is.element(snpsXc, snpsYc)]
(px <- length(unique(snpsXnc))) # 639 = 1439 - 800
(frxnc <- table(table(snpsXnc)))
# #Interactions involved    1   2    3   4   5   6
# #SNPs                   478 103  44  10   3   1
sum(frxnc) == px # TRUE # 639

p <- 6456
prop.test( c(py, px), rep(p, 2),  alternative = "greater" )
    # P=0.0024  (result present in the manuscript)


#### Wilcoxon-Mann-Whitney test ####
alpha <- .05
indYc <- apply(hitsYc, 1, ijtoi)
length(indYc) # 1272
all(tX.BH[indYc] < alpha) # TRUE
pYmw <- tY.BH[indYc]
all(pYmw > alpha) # TRUE
rktY.BH <- rank(tY.BH) # 10 seconds
rktX.BH <- rank(tX.BH) # 10 seconds
rdYmw <- rktY.BH[indYc] - rktX.BH[indYc]
length(rdYmw) # 1272
summary(rdYmw)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 932  6131697  6132324  9278838 13267413 18981048

indXc <- apply(hitsXc, 1, ijtoi)
length(indXc) # 1137
all(tY.BH[indXc] < alpha) # TRUE
pXmw <- tX.BH[indXc]
all(pXmw > alpha) # TRUE
rdXmw <- rktX.BH[indXc] - rktY.BH[indXc]
length(rdXmw) # 1137
summary(rdXmw)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 1328  6029136  6029794  9367372 13632450 18978566

wilcox.test(rdYmw, rdXmw) # P < 2.2e-16
wilcox.test(rdYmw, rdXmw, alternative = "less") # P = 1
wilcox.test(rdYmw, rdXmw, alternative = "greater") # P < 2.2e-16
# The significant difference is in the desired direction
# (result present in the manuscript)

rm(tX.BH)
rm(tY.BH) # Freeing up memory
