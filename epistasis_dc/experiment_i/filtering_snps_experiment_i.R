# One would run the present script with the "missSNP_frq005" PLINK
# files (*.bed, *.bim, and *.fam) in the same folder.
# "plink.exe" should also be at the same location.

# Note that we are not allowed to share the original data from the
# whole GWA study. Therefore, the present script is only an illustration
# of the commands that one would need to use in PLINK to perform the
# filters we describe in the supplement.

## Removing high LD and alosomes:
system2("plink", c( "--bfile missSNP_frq005",
  "--exclude range highld_redux.txt",
  # "--autosome", # All are autosomal
  "--make-bed",
  "--out miss1"
)) # 18268

## Pruning R^2<0.1
system2("plink", c( "--bfile miss1",
  "--indep-pairwise 500 1 0.1",
  "--out survivors_pruning"
)) # 8030
system2("plink", c( "--bfile miss1",
  "--extract survivors_pruning.prune.in",
  "--make-bed",
  "--out miss_final"
))

## Translating into 0,1,2
system2("plink", c( "--bfile miss_final",
  "--recode A",
  "--out miss_final_012"
)) # 8030
# No alosomes => No warnings about heterozygosity


## Reading the RAW file:
D <- read.delim("miss_final_012.raw",
  header = T, sep = " ", skip = 0
)
dim(D) # n=585+573=1158 x L+6=8030+6=8036
D[1:9, 1:9]
pheno <- D[, 6]
length(pheno)
table(pheno) # 1:=CO ; 2:=CA
D <- as.matrix(D[, -(1:6)])
dim(D) # n=585+573=1158 x L=8030
D[1:9, 1:5]
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


## SNP name blinding:
colnames(D) <- NULL # Confidentiality
class(D)
D[1:5, 1:5]
length(names) # L=8030
write.table(t(names), "miss_snp_names.dat",
  row.names = F, col.names = F, quote = T, eol = "\r\n"
)
nread <- unname(drop(as.matrix(read.table("miss_snp_names.dat"))))
identical(nread, names) # TRUE
length(nread) # L=8030
rm(nread)

## Phenotype:
identical(sort(unique(pheno)), 1:2) # 1:=CO ; 2:=CA
head(pheno)
case <- pheno == 2
sum(case)
sum(!case) # nCA=585; nCO=573

## Blind matrices:
# X=case, Y=control (as I always do)
X <- D[case, ]
Y <- D[!case, ]
dim(X) # nCA=585 x L=8030
dim(Y) # nCO=573 x L=8030
rm(D)
rm(case)

## Exporting the blind matrices:
write.table(X, "Matrix_X.dat", row.names = F, col.names = F)
Xread <- as.matrix(read.table("Matrix_X.dat"))
all(dim(X) == dim(Xread)) & all(X == Xread) # TRUE
unique(c(dim(X), dim(Xread))) # nCA=585 x L=8030
rm(Xread)

write.table(Y, "Matrix_Y.dat", row.names = F, col.names = F)
Yread <- as.matrix(read.table("Matrix_Y.dat"))
all(dim(Y) == dim(Yread)) & all(Y == Yread) # TRUE
unique(c(dim(Y), dim(Yread))) # nCO=573 x L=8030
rm(Yread)

## Chromosome and position:
chr_pos <- read.table("miss_final.bim")
chr_pos <- as.matrix(chr_pos[, c(1, 4)])
colnames(chr_pos) <- NULL
write.table(chr_pos, "chr_pos.dat", row.names = F, col.names = F)

# We then did some minor changes in "chr_pos.dat" for the purpose
# of sharing our reproducibility material, in a way that preserves
# anonymity of sampled individuals and maintains the results we present.
