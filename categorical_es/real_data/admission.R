rm(list = ls())
source("../../simulations/test_functions_ct_dcov.R")
source("../../simulations/bs_test_functions.R")
# load("admission_data.Rdata")
# datos.continx
# colnames(datos.continx)<-c("Chronicity","PRS")
# write.table(datos.continx,"admission_data.txt",row.names=F)
admission <- read.table("admission_data.txt", header = T)
# View(tail(admission))
attach(admission)
qq <- findInterval(PRS, quantile(x = PRS, probs = seq(0, 1, .1)),
  rightmost.closed = T
)
(ctp <- table(Chronicity, qq, dnn = c("Chronicity", "DecilePRS")))
cor.test(Chronicity, qq, method = "ke") # 0.021
as.matrix(ctp)
dimnames(ctp) <- NULL
ctp
dim(ctp) # 4 x 10
test.chisq(ctp) # 0.1881408
chisq.test(ctp) # Same as above, but warning that the CS approx may be wrong
set.seed(1)
ChisqPerm(ctp, 999) # 0.178
set.seed(1)
chiperm(ctp) # 0.178
set.seed(1)
chisq.test(ctp, simulate.p.value = T, B = 999)$p.value # 0.178
test.g(ctp) # 0.2083935
pchisq(GStat(ctp), df = (4 - 1) * (10 - 1), lower.tail = F) # Same as above
set.seed(1)
fisher.test(ctp, simulate.p.value = T, B = 999)$p.value # 0.223
# fisher.test(ctp,B=999)$p.value # Explodes
USPDiscrete(ctp, 999)$p.value # 0.211
test.dcov(ctp) # 0.186906


qq <- findInterval(PRS, quantile(x = PRS, probs = seq(0, 1, .2)),
  rightmost.closed = T
)
(ctp <- table(Chronicity, qq, dnn = c("Chronicity", "DecilePRS")))
cor.test(Chronicity, qq, method = "ke") # 0.04
as.matrix(ctp)
dimnames(ctp) <- NULL
ctp
dim(ctp) # 4 x 5
test.chisq(ctp) # 0.1897693
chisq.test(ctp) # Same as above, but warning that the CS approx may be wrong
set.seed(1)
ChisqPerm(ctp, 999) # 0.192
set.seed(1)
chiperm(ctp) # Same
set.seed(1)
chisq.test(ctp, simulate.p.value = T, B = 999)$p.value # Same
test.g(ctp) # 0.1648958
pchisq(GStat(ctp), df = (nrow(ctp) - 1) * (ncol(ctp) - 1), lower.tail = F) # Same as above
set.seed(1)
fisher.test(ctp, simulate.p.value = T, B = 999)$p.value # 0.183
# fisher.test(ctp,B=999)$p.value # Explodes
USPDiscrete(ctp, 999)$p.value # 0.257
test.dcov(ctp) # 0.2463434

qq <- findInterval(PRS, quantile(x = PRS, probs = seq(0, 1, 1 / 3)),
  rightmost.closed = T
)
(ctp <- table(Chronicity, qq, dnn = c("Chronicity", "DecilePRS")))
cor.test(Chronicity, qq, method = "ke") # 0.02784
as.matrix(ctp)
dimnames(ctp) <- NULL
ctp
dim(ctp) # 4 x 5
test.chisq(ctp) # 0.1897693
chisq.test(ctp) # Same as above, but warning that the CS approx may be wrong
set.seed(1)
ChisqPerm(ctp, 999) # 0.192
set.seed(1)
chiperm(ctp) # Same
set.seed(1)
chisq.test(ctp, simulate.p.value = T, B = 999)$p.value # Same
test.g(ctp) # 0.1648958
pchisq(GStat(ctp), df = (nrow(ctp) - 1) * (ncol(ctp) - 1), lower.tail = F) # Same as above
set.seed(1)
fisher.test(ctp, simulate.p.value = T, B = 999)$p.value # 0.183
# fisher.test(ctp,B=999)$p.value # Explodes
USPDiscrete(ctp, 999)$p.value # 0.257
test.dcov(ctp) # 0.2463434

# Example Agresti (Table 2.6):
# ctp<-matrix(c(17066, 14464, 788,126,37,48,38,5,1,1),5,2) # Best with CHS
# ctp
# test.chisq(ctp) # 0.0167514
# chisq.test(ctp) # 0.01675
# set.seed(1);ChisqPerm(ctp,999) # 0.046
# set.seed(1);chiperm(ctp) # 0.046
# set.seed(1);chisq.test(ctp,simulate.p.value=T,B=999)$p.value # 0.046
# test.g(ctp) # 0.1845623
# pchisq(GStat(ctp), df=(nrow(ctp)-1)*(ncol(ctp)-1), lower.tail=F) # Same as above
# set.seed(1);fisher.test(ctp,simulate.p.value=T,B=999)$p.value # 0.048
# # fisher.test(ctp,B=999)$p.value # Explodes
# USPDiscrete(ctp,999)$p.value # 0.614
# test.dcov(ctp) # 0.5450086
