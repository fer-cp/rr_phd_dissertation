# Filtering our data with GTEx, as explained in the supplement.
# See the README of the reproducibility files for more details.

#### "brain" variants ####
files <- dir("gtex_v7_signifpairs/brain")
VB <- character(0)
smash <- function(x) unlist(strsplit(x, "_"))
monolit <- function(x) nchar(x) == 1
for (i in 1:length(files)) {
  fi <- files[i]
  D <- read.delim(paste("gtex_v7_signifpairs/brain/", fi, sep = ""))
  vars <- unique(D[, 1])
  rm(D)
  nv <- length(vars)
  V <- t(sapply(vars, smash))
  if (!all.equal(dim(V), c(nv, 5))) cat("Problems...\n")
  x3 <- sapply(V[, 3], monolit)
  x4 <- sapply(V[, 4], monolit)
  V <- V[x3 & x4, 1:4]
  fifi <- unlist(strsplit(fi, ".signifpairs.txt"))
  V <- cbind(V, rep(fifi, nrow(V)))
  rownames(V) <- NULL
  colnames(V) <- c("chr", "pos", "wt", "minor", "tissue")
  VB <- rbind(VB, V)
  rm(list = c("V", "vars", "x3", "x4"))
  # View(VB)
  cat("Done for tissue #", i, "out of", length(files), "\n")
}
View(VB)
length(unique(VB[, "tissue"])) == length(files)
unique(VB[, "tissue"])
dim(VB) # 2 791 641 x 5
write.table(VB, "snps_brain_nonexclusive.dat", col.names = F, row.names = F)
rm(VB)


#### "nonbrain" variants ####
# Much slower than "brain"
rm(list = ls())
files <- dir("gtex_v7_signifpairs/nonbrain")
VNB <- character(0)
smash <- function(x) unlist(strsplit(x, "_"))
monolit <- function(x) nchar(x) == 1
for (i in 1:length(files)) {
  fi <- files[i]
  D <- read.delim(paste("gtex_v7_signifpairs/nonbrain/", fi, sep = ""))
  vars <- unique(D[, 1])
  rm(D)
  nv <- length(vars)
  V <- t(sapply(vars, smash))
  if (!all.equal(dim(V), c(nv, 5))) cat("Problems...\n")
  x3 <- sapply(V[, 3], monolit)
  x4 <- sapply(V[, 4], monolit)
  # WW<-cbind(V[,3:4],x3&x4);View(WW)
  V <- V[x3 & x4, 1:4]
  fifi <- unlist(strsplit(fi, ".signifpairs.txt"))
  V <- cbind(V, rep(fifi, nrow(V)))
  rownames(V) <- NULL
  colnames(V) <- c("chr", "pos", "wt", "minor", "tissue")
  VNB <- rbind(VNB, V)
  rm(list = c("V", "vars", "x3", "x4"))
  # View(VNB)
  cat("Done for tissue #", i, "out of", length(files), "\n")
}
View(VNB)
length(unique(VNB[, "tissue"])) == length(files)
unique(VNB[, "tissue"])
dim(VNB) # 18 153 791 x 5
write.table(VNB, "snps_nonbrain.dat", col.names = F, row.names = F)
rm(VNB)

#### "brain" minus "nonbrain" ####
B <- read.table("snps_brain_nonexclusive.dat")
View(B[1:5, ])
cp <- apply(B[, 1:2], 1, paste, collapse = ":")
wospace <- function(x) paste(unlist(strsplit(x, " ")), collapse = "")
cp <- sapply(cp, wospace)
B <- cbind(B[, 1:2], cp, B[, 3:5])
rm(cp)
View(B[1:1000, ])
#
NB <- read.table("snps_nonbrain.dat")
View(NB[1:5, ])
cp <- apply(NB[, 1:2], 1, paste, collapse = ":")
wospace <- function(x) paste(unlist(strsplit(x, " ")), collapse = "")
cp <- sapply(cp, wospace)
NB <- cbind(NB[, 1:2], cp, NB[, 3:5])
rm(cp)
View(NB[1:1000, ])
#
cpB <- B[, 3]
cpNB <- NB[, 3]
excl <- is.element(cpB, cpNB)
length(cpB) # 2791641
length(excl) # 2791641
sum(excl) # 2658740
keep <- which(!excl)
length(keep) # 132901
length(B[keep, 3]) # 132901
length(unique(B[keep, 3])) # 97913

rm(NB)
B <- B[keep, ]
dim(B) # 132901 x 6
write.table(B, "snps_brain_only.dat", col.names = F, row.names = F)

#### High LD ####
x <- read.table("highld.txt")
x <- x[, c(2:4, 1)]
x <- rbind(x, c(23, 43200000, "45000000", "Region25"))
write.table(x, "highld_redux.txt", col.names = F, row.names = F, quote = F)

#### Excluding/extracting ####
snps_brain_only <- read.table("snps_brain_only.dat", header = F)
snpstoextract <- unique(snps_brain_only[, 3])
length(snpstoextract) # 97913
write.table(snpstoextract, "snpstoextract.txt",
  quote = F,
  col.names = F, row.names = F
)
toex <- read.table("snpstoextract.txt", quote = "", header = F)
# class(toex);head(toex)
length(drop(as.matrix(toex))) # 97913
system2("plink", c(
  "--bfile gwas_scz",
  "--extract snpstoextract.txt",
  "--exclude range highld_redux.txt",
  "--autosome",
  "--make-bed",
  "--out gwas_scz_brainonly"
)) # 56395


#### Pruning R^2<0.1 ####
system2("plink", c(
  "--bfile gwas_scz_brainonly",
  "--indep-pairwise 500 1 0.1",
  "--out survivors_pruning"
)) # 6456
system2("plink", c(
  "--bfile gwas_scz_brainonly",
  "--extract survivors_pruning.prune.in",
  "--make-bed",
  "--out data_final"
))

#### Translating into 0,1,2 ####
system2("plink", c(
  "--bfile data_final",
  "--recode A",
  "--out data_final_012"
))
# No alosomes => No warnings about heterozygosity
# Now we continue in the main script for Experiment II


#### Reading the RAW file ####
D <- read.delim("data_final_012.raw",  header = T, sep = " ", skip = 0 )
dim(D) # n=585+573=1158 x L+6=6456+6=6462
D[1:6, 1:6]
pheno <- D[, 6]
length(pheno)
table(pheno) # 1:=CO ; 2:=CA
D <- as.matrix(D[, -(1:6)])
dim(D) # n=585+573=1158 x L=6456
D[1:4, 1:6]
names <- colnames(D)
head(names)
names <- gsub("[.]", ":", gsub("X", "", names))
head(names)
for (B in c("A", "T", "G", "C")) {
  names <- gsub(paste("_", B, sep = ""), "", names)
}
head(names)


## No missing data:
all.equal(sort(unique(as.numeric(D)), na.last = T), 0:2) # TRUE


#### SNP name blinding ####
colnames(D) <- NULL # Confidentiality
class(D)
D[1:5, 1:5]
length(names) # L=6456
write.table(t(names), "snp_names.dat",
  row.names = F, col.names = F, quote = T, eol = "\r\n" )
nread <- unname(drop(as.matrix(read.table("snp_names.dat"))))
identical(nread, names) # TRUE
rm(names)
rm(nread)

#### Phenotype ####
identical(sort(unique(pheno)), 1:2) # 1:=CO ; 2:=CA
head(pheno)
case <- pheno == 2
sum(case)
sum(!case) # nCA=585; nCO=573
rm(pheno)

#### Blind matrices ####
# X=case, Y=control
X <- D[case, ]
Y <- D[!case, ]
dim(X) # nCA=585 x L=6456
dim(Y) # nCO=573 x L=6456
rm(D)
rm(case)

## Exporting the blind matrices:
write.table(X, "Matrix_X.dat", row.names = F, col.names = F)
Xread <- as.matrix(read.table("Matrix_X.dat"))
all(dim(X) == dim(Xread)) & all(X == Xread) # TRUE
unique(c(dim(X), dim(Xread))) # nCA=585 x L=6456
rm(Xread)

write.table(Y, "Matrix_Y.dat", row.names = F, col.names = F)
Yread <- as.matrix(read.table("Matrix_Y.dat"))
all(dim(Y) == dim(Yread)) & all(Y == Yread) # TRUE
unique(c(dim(Y), dim(Yread))) # nCO=573 x L=6456
rm(Yread)

### Chromosome and position ####
chr_pos <- read.table("data_final.bim")
chr_pos <- as.matrix(chr_pos[, c(1, 4)])
colnames(chr_pos) <- NULL
dim(chr_pos) # 6454 x 2
head(chr_pos)
tail(chr_pos)
write.table(chr_pos, "chr_pos.dat", col.names = F, row.names = F)
