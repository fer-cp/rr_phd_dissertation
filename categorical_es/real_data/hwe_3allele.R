rm(list = ls())

require(data.table)
source("../../simulations/test_functions_ct_dcov.R")

N <- 427
interesting <- drop(as.matrix(read.table(
  "pgc3_snps.txt"
)))

## SNP: rs2594292 / 1:11618451 (GRCh37)
# The allele frequencies are A>G>T
snp <- "rs2594292"
snp %in% interesting # FALSE

# gnomAD European non-Finnish
t1 <- .69 # f(A)=0.69
t2 <- .26 # f(G)=0.26
t3 <- .05 # f(T)=0.05
t1 + t2 + t3 == 1 # TRUE
(o <- c(AA = 214, AG = 148, GG = 34, AT = 16, TG = 15, TT = 0))
sum(o) == N # TRUE (427)
p0 <- c(t1^2, 2 * t1 * t2, t2^2, 2 * t1 * t3, 2 * t2 * t3, t3^2)
(e <- N * p0)
test.e(o, p0) # 0.2365832
test.chisq.gof(o, p0) # 0.06865777
