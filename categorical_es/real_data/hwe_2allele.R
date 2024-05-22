rm(list = ls())

require(data.table)
source("../../simulations/test_functions_ct_dcov_v4.R")

bim <- fread("hwe_example.bim")
interesting <- drop(as.matrix(read.table("pgc3_snps.txt")))

system2("plink", c(
  "--bfile hwe_example",
  "--recode A",
  "--out hwe_example"
))
system2("plink", c(
  "--bfile hwe_example",
  "--recode",
  "--tab",
  "--out hwe_example"
))
View(bim)
ped <- fread("hwe_example.ped")
View(ped)
(N <- nrow(ped)) # 427
ped <- as.matrix(ped)
ped <- ped[, -c(1:6)]
ncol(ped) # 11



## SNP: rs9545047 / 13:79859456:A:C
snp <- "rs9545047"
maf <- .41
snp %in% interesting # TRUE
t <- 1 - maf
gg <- ped[, 7]
(o <- table(gg))
p0 <- c(t^2, 2 * t * (1 - t), (1 - t)^2)
(e <- N * p0)
test.e(o, p0) # 0.027
test.chisq.gof(o, p0) # 0.027
