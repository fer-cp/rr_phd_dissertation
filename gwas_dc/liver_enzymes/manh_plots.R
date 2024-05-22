# rm(list=ls())
require(data.table) # fread
require(qqman) # manhattan

chrs <- matrix(c(rep(0, 9), 1:9), ncol = 2)
(chrs <- apply(chrs, 1, paste, collapse = ""))
(chrs <- c(chrs, 10:22))
bim <- fread("trinity.bim")

Agesx <- T
anova <- F # Now ANOVA is compatible with AgeSx
common <- F

#### Producing Manhattan plots ####
for (enz in c("ggt")) { # "alt","ast",
  brange <- 1:4
  if (anova) brange <- 0
  for (bb in brange) {
    l <- character(0)
    tit <- paste("DC-snptest for enzyme ", toupper(enz), " with b=", bb, sep = "")
    if (anova) tit <- paste("ANOVA test for enzyme ", toupper(enz), sep = "")
    if (Agesx) tit <- paste(tit, ",\n correcting for Sex and Age", sep = "")
    if (common) tit <- paste(tit, ",\n (only for SNPs with MAF>0.1)", sep = "")
    # cat("Creating the matrix for plotting",enz,"with b = ",bb,"\n")
    for (i in chrs) {
      ff <- paste("pv/pv_", enz, "_chr", i, "_b", bb, # ".txt",
        sep = ""
      )
      if (Agesx) ff <- paste("pv_agesx/pv_", enz, "_chr", i, "_b", bb, "_agesx.txt", sep = "")
      if (anova) ff <- paste("pv_competing/pv_", enz, "_chr", i, "_anova_agesx.txt", sep = "")
      pp <- as.matrix(fread(ff)) # Instant
      dim(pp) # p x 2
      p <- nrow(pp)
      snps <- pp[, 1]
      pv <- as.numeric(pp[, 2])
      pv[is.na(pv)] <- 1 # I do not deal with NA's for plotting, etc.
      l <- rbind(l, cbind(snps, pv, rep(i, p)))
    }
    MM <- cbind(l, bim[, -3])
    MM <- as.data.frame(MM)
    colnames(MM) <- c("SNP_MA", "P", "CHR_bis", "CHR", "SNP", "BP", "MA", "mA")
    MM$P <- as.numeric(MM$P)
    nrow(MM) # 757577
    if (common) {
      cs <- drop(as.matrix(read.table("common/cs.txt")))
      ics <- is.element(MM$SNP, cs)
      sum(ics) # 592809
      MM.old <- MM
      MM <- MM[ics, ]
    }
    if (!anova) cat("GWS SNPs (P<1e-5) with b = ", bb, " for enzyme", enz, ":\n") # 5e-8
    if (anova) cat("GWS SNPs (P<1e-5) for ANOVA for enzyme", enz, ":\n") # 5e-8
    if (common) cat("(only common SNPs are displayed)\n")
    outp <- MM[MM$P < 1e-5, c(4:8, 2)]
    print(outp[order(outp[, 6]), ]) # 5e-8
    # fm<-paste("manh/mp_",enz,"_b",bb,".pdf",sep="")
    fm <- paste("manh/mp_", enz, "_b", bb, sep = "")
    if (anova) fm <- paste(fm, "_anova", sep = "")
    if (Agesx) fm <- paste(fm, "_agesx", sep = "")
    if (common) fm <- paste(fm, "_common", sep = "")
    fm <- paste(fm, ".png", sep = "")
    # pdf(fm,width=7,height=7/1.54)
    # tit=paste("b = ",bb,sep="")
    # png(fm,
    #     bg="transparent",type="windows",width=1587,height=934)
    # 1477 x 1207
    # if(!anova)cat("Plotting for enzyme",enz,"with b = ",bb,"\n")
    # if(anova)cat("Plotting for enzyme",enz,"with method anova\n")
    # par(mfrow=c(2,2))
    # manhattan(MM[,c(4,6,2,5)],main=tit)
    # Blue line: 1E-5; Red line: 5E-8
    # dev.off()
  }
}

# Maybe for mdofying the height of the red line
# 0.05/nrow(MM) # 6.6E-8
