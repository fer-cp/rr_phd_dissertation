# We begin with the 8030 SNPs at the "miss_final" database
# which in turn was generated with the "filtering_snps_experiment_i.R"
# using the whole GWAS data that we are not allowed to share.

# At many points of the script, we generate intermediate result files,
# in order to ease running only parts of it.

source("../sim_functions_snp_snp.R")
X <- as.matrix(read.table("Matrix_X.dat"))
Y <- as.matrix(read.table("Matrix_Y.dat"))

tY <- test1(Y) # matrix
(p <- ncol(Y)) # 8030
rm(Y)
write.table(tY, "pval_Y_brute.dat", col.names = F, row.names = F)
rm(tY)

tX <- test1(X) # 8030 x 8030 symmetric matrix, with NA's on the diagonal
p == ncol(X) # TRUE
rm(X)
write.table(tX, "pval_X_brute.dat", col.names = F, row.names = F)
rm(tX)


#### FDR correction by BH ####
tYc <- as.matrix(read.table("pval_Y_brute.dat")) # DF to matrix
tYc[lower.tri(tYc)] <- NA # triu
tYc[1:5, 1:5]
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

#### Counting pairs ####
tX.BH <- drop(as.matrix(read.table("tX_BH.dat")))
tY.BH <- drop(as.matrix(read.table("tY_BH.dat")))
alpha <- .05
table(tX.BH < alpha, tY.BH < alpha) # Rows correspond to tX; cols to tY
#            FALSE     TRUE
# FALSE  32234729      501
# TRUE        615      590

indY <- which(tX.BH < alpha & tY.BH > alpha)
mean(tX.BH < alpha & tY.BH > alpha) # 0.0019 % of pairs are in presumed epistasis
                                    # 615 out of (p^2-p)/2
length(indY) # 615
p <- 8030
max(indY) / (p^2 - p) * 2 # 0.9995019


indX <- which(tX.BH > alpha & tY.BH < alpha) # FPs
mean(tX.BH > alpha & tY.BH < alpha) # 0.0016 % of pairs are apparent false positives
                                    # 501 out of (p^2-p)/2
length(indX) # 501
max(indX) / (p^2 - p) * 2 # 0.9997509

prop.test(c(length(indY), length(indX)),
  rep((p^2 - p) / 2, 2),
  alternative = "greater"
) # P = 0.000359

indA <- which(tX.BH < alpha & tY.BH < alpha) # Ancestrality
mean(tX.BH < alpha & tY.BH < alpha) # 0.0018 % of pairs reflect ancestrality
                                    # 590 out of (p^2-p)/2
length(indA) # 590
max(indA) / (p^2 - p) * 2 # 0.9997507

indN <- which(tX.BH > alpha & tY.BH > alpha)
mean(tX.BH > alpha & tY.BH > alpha) # 99.99471 % of pairs show no signal whatsoever
length(indN) # 32234729
max(indN) / (p^2 - p) * 2 # 1

# Warning: The names "indX"/"indY", "hitsX"/"hitsY" can be confusing.
# Check twice what is what before making any changes.

rm(tX.BH)
rm(tY.BH)

#### Auxiliary functions for dealing with triangular matrices ####
source("../aux_functions_triu.R") # Working with triangular matrices

#### Counting SNPs and SNP pairs ####
hitsY <- t(sapply(indY, itoij)) # Pairs by row
nrhy <- nrow(hitsY)
hitsY[c(1:4, seq(nrhy, nrhy - 2, -1)), ]
length(unique(as.numeric(hitsY))) # 1138
write.table(hitsY, "snppairs_hits_Y.dat", row.names = F, col.names = F)

hitsX <- t(sapply(indX, itoij)) # Pairs by row
nrhx <- nrow(hitsX)
hitsX[c(1:4, seq(nrhx, nrhx - 2, -1)), ]
length(unique(as.numeric(hitsX))) # 934
write.table(hitsX, "snppairs_hits_X.dat", row.names = F, col.names = F)

hitsA <- t(sapply(indA, itoij))
nrha <- nrow(hitsA)
hitsA[c(1:4, seq(nrha, nrha - 2, -1)), ]
length(unique(as.numeric(hitsA))) # 1094
write.table(hitsA, "snppairs_hits_A.dat", row.names = F, col.names = F)

hitsN <- t(sapply(indN, itoij)) # Pairs by row
(nrhn <- nrow(hitsN)) # 32234729
hitsN[c(1:4, seq(nrhn, nrhn - 2, -1)), ]
length(unique(as.numeric(hitsN))) # 8030 (every SNP is involved at least one TN)
write.table(hitsN, "snppairs_hits_N.dat", row.names = F, col.names = F)

prop.test(
  c(
    length(unique(as.numeric(hitsY))),
    length(unique(as.numeric(hitsX)))
  ),
  rep(p, 2),
  alternative = "greater"
) # P = 8.8E-7

# Chromosome and position:
chr_pos <- as.matrix(read.table("chr_pos.dat"))
dim(chr_pos) # 8030 x 2
head(chr_pos)
tail(chr_pos)
  # Please note that some SNP IDs have been altered for the purpose of
  # sharing this script for reproducibility, in a way that preserves
  # anonymity of sampled individuals and maintains the results we present.


# Functions for testing "nearness"
samechr <- function(h, cp = chr_pos) {
  if (length(h) != 2) {
    return(NA)
    cat("Trouble...\n")
  }
  return(length(unique(cp[h, 1])) == 1)
}
nearby <- function(h, cp = chr_pos) { # Flags candidates to LD exclusion
  if (length(h) != 2) {
    return(NA)
    cat("Trouble...\n")
  }
  if (samechr(h) & abs(cp[h[1], 2] - cp[h[2], 2]) < 1E6) {
    return(T)
  }
  return(F)
}

#### Filtering by physical distance within the genome ####
# Studying "nearness" on the TPs
hitsY <- as.matrix(read.table("snppairs_hits_Y.dat"))
dim(hitsY) # 615   2
mean(apply(hitsY, 1, samechr)) # 83.4 % our hits lay on the same chromosome
sum(apply(hitsY, 1, samechr)) # 513 out of 615

hitsYsc <- hitsY[apply(hitsY, 1, samechr), ]
mean(apply(hitsYsc, 1, nearby)) # 97.9 % are < 1 Mb away
sum(apply(hitsYsc, 1, nearby)) # 502 out of 513

# Studying "nearness" on the FPs
hitsX <- as.matrix(read.table("snppairs_hits_X.dat"))
dim(hitsX) # 501   2
mean(apply(hitsX, 1, samechr)) # 81.8 % our hits lay on the same chromosome
sum(apply(hitsX, 1, samechr)) # 410 out of 501

hitsXsc <- hitsX[apply(hitsX, 1, samechr), ]
mean(apply(hitsXsc, 1, nearby)) # 99.0 % are < 1 Mb away
sum(apply(hitsXsc, 1, nearby)) # 406 out of 410

# Studying "nearness" on the APs
hitsA <- as.matrix(read.table("snppairs_hits_A.dat"))
dim(hitsA) # 590   2
mean(apply(hitsA, 1, samechr)) # 100.0 % our hits lay on the same chromosome
sum(apply(hitsA, 1, samechr)) # 590 out of 590

hitsAsc <- hitsA[apply(hitsA, 1, samechr), ]
mean(apply(hitsAsc, 1, nearby)) # 99.9 % are < 1 Mb away
sum(apply(hitsAsc, 1, nearby)) # 589 out of 590

# Studying "nearness" on the TNs
hitsN <- as.matrix(read.table("snppairs_hits_N.dat"))
dim(hitsN) # 32234729   2
mean(apply(hitsN, 1, samechr)) # 5.29 % our hits lay on the same chromosome
sum(apply(hitsN, 1, samechr)) # 1 706 096 out of 32 234 729

hitsNsc <- hitsN[apply(hitsN, 1, samechr), ]
mean(apply(hitsNsc, 1, nearby)) # 2.37 % are < 1 Mb away
sum(apply(hitsNsc, 1, nearby)) # 40 431 out of 1 706 096


#### Filtering by physical distance ####
hitsYc <- hitsY[!apply(hitsY, 1, nearby), ]
nrow(hitsYc) #  113 (result present in the manuscript)
snpsY <- unique(as.numeric(as.matrix(hitsYc)))
length(snpsY) # 222 < 2*113 (result present in the manuscript)

hitsXc <- hitsX[!apply(hitsX, 1, nearby), ]
nrow(hitsXc) # 95 (result present in the manuscript)
snpsX <- unique(as.numeric(as.matrix(hitsXc)))
length(snpsX) # 189 < 2*95 (result present in the manuscript)

hitsAc <- matrix(hitsA[!apply(hitsA, 1, nearby), ], ncol = 2, byrow = T)
nrow(hitsAc) # 1
snpsA <- unique(as.numeric(as.matrix(hitsAc)))
length(snpsA) # 2

hitsNc <- hitsN[!apply(hitsN, 1, nearby), ]
nrow(hitsNc) # 32194298
snpsN <- unique(as.numeric(as.matrix(hitsNc)))
length(snpsN) # 8030

(t <- nrow(hitsYc) + nrow(hitsXc) + nrow(hitsAc) + nrow(hitsNc))
prop.test( c(nrow(hitsYc), nrow(hitsXc)), rep(t, 2),
          alternative = "greater" ) # 0.1193  (result present in the manuscript)
prop.test(c(length(snpsY), length(snpsX)), rep(p, 2),
          alternative = "greater" ) # 0.05491  (result present in the manuscript)

#### And we map all into genes ####
chrpos <- apply(chr_pos, 1, paste, collapse = ":")
length(chrpos) # 8030
snpsX.cp <- chrpos[snpsX]
snpsY.cp <- chrpos[snpsY]
write.table(chrpos, "chrpos.dat", row.names = F, col.names = F)
write.table(snpsX, "snpsX.dat", row.names = F, col.names = F)
write.table(snpsY, "snpsY.dat", row.names = F, col.names = F)
write.table(snpsX.cp, "snpsX.cp.dat", row.names = F, col.names = F)
write.table(snpsY.cp, "snpsY.cp.dat", row.names = F, col.names = F)

snpsX.cp <- as.character(as.matrix(read.table("snpsX.cp.dat")))
snpsY.cp <- as.character(as.matrix(read.table("snpsY.cp.dat")))
cp <- as.character(as.matrix(read.table("chrpos.dat")))
chrpos <- cp

#### A few auxiliary functions ####
atomize <- function(x) unlist(strsplit(x, "")) # Smashes character strings
is.only.numbers <- function(x) { # Does the arg contain only digits?
  return(all(is.element(atomize(x), as.character(0:9))))
}
is.gene <- function(x) { # Is the arg an ENSEMBL gene ID?
  if (!is.character(x)) {
    return(NA)
    cat("Trouble...\n")
  }
  t <- unlist(strsplit(x, "G"))
  if (t[1] != "ENS") {
    return(F)
  }
  return(is.only.numbers(t[2]))
}
is.lrg <- function(x) { # Is the arg an ENSEMBL regulatory region ID?
  if (!is.character(x)) {
    return(NA)
    cat("Trouble...\n")
  }
  t <- unlist(strsplit(x, "_"))
  if (t[1] != "LRG") {
    return(F)
  }
  return(is.only.numbers(t[2]))
}
is.snp <- function(x) { # Is the arg a SNP? (Format chr:pos)
  if (!is.character(x)) {
    return(NA)
    cat("Trouble...\n")
  }
  t <- unlist(strsplit(x, ":"))
  q1 <- is.element(t[1], as.character(1:22))
  q2 <- is.only.numbers(t[2])
  return(q1 & q2)
}

snp2gene <- as.matrix(read.table("snp2gene.dat"))
colnames(snp2gene) <- NULL
snp2gene <- snp2gene[!sapply(snp2gene[, 2], is.lrg), ]
all(sapply(snp2gene[, 2], is.gene)) # TRUE
all(sapply(snp2gene[, 1], is.snp)) # TRUE
all(sapply(chrpos, is.snp)) # TRUE

dup2gene <- as.matrix(read.table("dup2gene.dat"))
colnames(dup2gene) <- NULL
all(sapply(dup2gene[, 2], is.gene)) # TRUE
all(sapply(dup2gene[, 1], is.snp)) # TRUE

all(sapply(snpsY.cp, is.snp)) # TRUE
ridY <- is.element(snp2gene[, 1], snpsY.cp)
genesYnu <- snp2gene[ridY, 2]
length(genesYnu) # 223
# Exactly one SNP is associated to two genes:
taux.genesY <- table(snp2gene[ridY, 1])
table(taux.genesY)
# 1   2
# 221 1
taux.genesY[taux.genesY == 2]
# 8:22449181
# 2
iigY <- snp2gene[ridY, 1] == "8:22449181"
ssgY <- snp2gene[ridY, ]
ssgY[iigY, ]
# [1,] "8:22449181" "ENSG00000120913"
# [2,] "8:22449181" "ENSG00000248235"
genesY <- unique(genesYnu)
rm(genesYnu)
(ngY <- length(genesY)) # 220 (result present in the manuscript)
length(unique(snp2gene[which(ridY), 1])) ==
  length(snpsY.cp)
length(snpsY.cp) # TRUE; 222

all(sapply(snpsX.cp, is.snp)) # TRUE
ridX <- is.element(snp2gene[, 1], snpsX.cp)
genesXnu <- snp2gene[ridX, 2]
length(genesXnu) # 192
# Exactly one SNP is associated to two genes:
taux.genesX <- table(snp2gene[ridX, 1])
table(taux.genesX)
# 1   2
# 221 1
taux.genesX[taux.genesX == 2]
# 11:14810762 2:234637853  3:57414071
# 2           2           2
iigX <- snp2gene[ridX, 1] %in% names(taux.genesX[taux.genesX == 2])
ssgX <- snp2gene[ridX, ]
ssgX[iigX, ]
# [1,] "2:234637853" "ENSG00000241635"
# [2,] "2:234637853" "ENSG00000243135"
# [3,] "11:14810762" "ENSG00000261923"
# [4,] "11:14810762" "ENSG00000152270"
# [5,] "3:57414071"  "ENSG00000174844"
# [6,] "3:57414071"  "ENSG00000262840"
genesX <- unique(genesXnu)
(ngX <- length(genesX)) # 191  (result present in the manuscript)
length(unique(snp2gene[which(ridX), 1])) ==
  length(snpsX.cp)
length(snpsX.cp) # TRUE; 189

any(is.element(snpsX.cp, dup2gene[, 1])) # FALSE
any(is.element(snpsY.cp, dup2gene[, 1])) # FALSE
all(is.element(snpsY.cp, snp2gene[, 1])) # TRUE
all(is.element(snpsX.cp, snp2gene[, 1])) # TRUE

all(sapply(genesX, is.gene)) # TRUE
all(sapply(genesY, is.gene)) # TRUE

## For the background, we need to map into genes
## all the $p$ initial SNPs

genes <- snp2gene[is.element(snp2gene[, 1], cp), 2]
length(genes) # 8151 > 8030
genes <- unique(genes)
(ng <- length(genes)) # 6217
length(unique(snp2gene[which(is.element(snp2gene[, 1], cp)), 1])) ==
  length(cp)
length(cp) # TRUE; 8030

prop.test(c(ngY, ngX), c(ng, ng), alternative = "greater") # P=0.080

any(is.element(cp, dup2gene[, 1])) # FALSE
all(is.element(cp, snp2gene[, 1])) # TRUE
all(sapply(genes, is.gene)) # TRUE
all(sapply(cp, is.snp)) # TRUE

write.table(genesX, "genesX.dat", row.names = F, col.names = F)
write.table(genesY, "genesY.dat", row.names = F, col.names = F)
write.table(genes, "genes.dat", row.names = F, col.names = F)
write.table(genesX, "genesX_noquote.dat",
  row.names = F, col.names = F, quote = F
)
write.table(genesY, "genesY_noquote.dat",
  row.names = F, col.names = F, quote = F
)

## Processing the HUGO names
genesYH <- as.matrix(read.table("genesY_HUGO.txt",
  header = T, sep = "\t"
))
nrow(genesYH) == ngY
genesYH <- genesYH[order(genesYH[, 2]), ]
write.table(genesYH, "genesY_HUGO_ordered.txt",
  row.names = F, col.names = F, quote = F
)
genesXH <- as.matrix(read.table("genesX_HUGO.txt",
  header = T, sep = "\t"
))
nrow(genesXH) == ngX
genesXH <- genesXH[order(genesXH[, 2]), ]
write.table(genesXH, "genesX_HUGO_ordered.txt",
  row.names = F, col.names = F, quote = F
)

## How many genes are common
genesX <- as.character(as.matrix(read.table("genesX.dat")))
genesY <- as.character(as.matrix(read.table("genesY.dat")))
length(genesX) # 191
length(genesY) # 220
length(intersect(genesX, genesY)) # 13
ncgX <- setdiff(genesX, genesY)
length(ncgX) # 191-13=178  (result present in the manuscript)
ncgY <- setdiff(genesY, genesX)
length(ncgY) # 220-13=207  (result present in the manuscript)

prop.test(c(length(ncgY), length(ncgX)),
          c(ng, ng),  alternative = "greater" )
          # P=0.074  (result present in the manuscript)

#### GSEA of synapse genes for the non-common genes ####
syngo <- drop(as.matrix(read.table("synapse_genes.txt",  quote = "")))
is.character(syngo) # TRUE
all(sapply(syngo, is.gene)) # TRUE
all(sapply(ncgY, is.gene)) # TRUE
all(sapply(ncgX, is.gene)) # TRUE
syngy <- intersect(syngo, ncgY)
length(syngy) # 13  (result present in the manuscript)
length(ncgY) # 207
syngx <- intersect(syngo, ncgX)
length(syngx) # 9  (result present in the manuscript)
length(ncgX) # 178
(ngsyn <- length(intersect(genes, syngo))) # 338

prop.test(c(length(syngy), length(syngx)),
          rep(ngsyn, 2),  alternative = "greater")
          # P=0.26  (result present in the manuscript)
