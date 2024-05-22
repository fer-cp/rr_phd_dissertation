rm(list = ls())

#### Auxiliary functions ####
require(data.table) # fread
source("sim_functions_snp_pheno.R") # Testing
Mode <- function(x) { # Statistical mode
  x <- x[!is.na(x)]
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}
fill <- function(v) { # Filler function
  m <- Mode(v)
  v[is.na(v)] <- m
  if (!is.element(m, 0:2)) cat("ERROR!\n")
  return(v)
}

#### Phenotype matrix ####
ph <- as.matrix(read.table("../pheno/pheno_matrix.txt",
  header = T, sep = ""
))
dim(ph) # 2407    5

#### Iterating for each chromosome ####
chrs <- matrix(c(rep(0, 9), 1:9), ncol = 2)
(chrs <- apply(chrs, 1, paste, collapse = ""))
(chrs <- c(chrs, 10:22))

for (i in chrs) {
  cat("Analysing chromosome", i, "\n")
  l_chr <- paste("--chr ", i, sep = "") # Extra zero is totally fine
  l_out <- paste("--out chr_files/trinity_chr", i, "_012", sep = "")
  system2("plink", c(
    "--bfile trinity",
    # l_chr,
    "--recode A",
    l_out
  ))
  fi <- paste("chr_files/trinity_chr", i, "_012.raw", sep = "")
  cat("Reading the RAW file for chromosome", i, "(15-60 s).\n")
  D <- as.matrix(fread(fi)) # 15-60 seconds
  subj <- D[, 1]
  ls <- length(subj)
  D <- D[, -(1:6)]
  snps <- colnames(D)
  cat("Converting to numeric (15 seconds)\n")
  D <- matrix(as.numeric(D), nrow = nrow(D), ncol = ncol(D))
  colnames(D) <- snps
  #
  ## Phenotype vector:
  for (enz in c("alt", "ast", "ggt")) {
    cat("Reading phenotype", enz, "for chromosome", i, "\n")
    if (enz == "ast") {
      Y.orig <- ph[, 1] # 1:=AST; 2:=ALT; 3:=GGT
    } else if (enz == "alt") {
      Y.orig <- ph[, 2] # 1:=AST; 2:=ALT; 3:=GGT
    } else {
      Y.orig <- ph[, 3] # 1:=AST; 2:=ALT; 3:=GGT
    }
    sx <- ph[, 4] - 1 # 0=male; 1=female
    sx <- sx[!is.na(Y.orig)]
    class(sx)
    age <- ph[, 5]
    age <- age[!is.na(Y.orig)]
    class(age)
    summary(age)
    Y.orig <- Y.orig[!is.na(Y.orig)]
    samy <- names(Y.orig)
    if (!all.equal(samy, names(sx))) cat("Oh, oh...\n")
    ly <- length(samy)
    #
    ## Intersecting with the individuals matching the phenotype samples:
    indiv <- intersect(samy, subj)
    li <- length(indiv)
    ii <- subj %in% indiv
    if (length(ii) != ls | sum(ii) != li) cat("ii wrong\n")
    X <- D[ii, ]
    if (nrow(X) != li) cat("X wrong\n")
    iy <- samy %in% indiv
    if (length(iy) != ly | sum(iy) != li) cat("iy wrong\n")
    if (any(Y.orig <= 0)) cat("Y not positive\n")
    Y <- log10(Y.orig[iy]) # Do not forget the log10!
    sx <- sx[iy]
    age <- age[iy]
    if (length(Y) != li) cat("Y wrong\n")
    if (length(sx) != li) cat("sx wrong\n")
    if (length(age) != li) cat("age wrong\n")
    CM <- cbind(sx, age)
    dim(CM) # 2k~ x 2
    #
    ## Filling:
    if (!any(is.na(X))) cat("Nothing to fill.\n")
    cat("Filling (15 seconds)\n")
    X <- apply(X, 2, fill)
    if (nrow(X) != li) cat("X wrong...\n")
    if (any(is.na(X))) cat("Filling went wrong.\n")
    ## Testing:
    for (bb in 0:4) {
      cat("Testing chr", i, "for", enz, "with b =", bb, "(15 seconds)\n")
      test <- snptest(X, Y, b = bb, covariates = NULL) # CM
      # Runs in 1 second per 5k SNPs
      # 25-50x slower with Sex as a covariate
      pv <- test$p.recalc
      if (length(pv) != ncol(D)) cat("Incorrect number of p-values.\n")
      if (any(is.na(pv))) {
        cat("Missing p-values.\n")
        isn <- which(is.na(pv))
        ww <- rep(F, length(isn))
        for (u in 1:length(isn)) {
          if (length(unique(X[, isn[u]])) == 1) ww[u] <- T
        }
        if (all(ww)) {
          cat("All because of a constant column\n")
          cat("It happened", length(isn), "times\n")
        } else {
          cat("Check what happened with SNPs number\n")
          print(isn[!ww])
        }
      }
      if (!all(pv > 0 & pv < 1, na.rm = T)) cat("Incorrect range of p-values.\n")
      ff <- paste("pv/pv_", enz, "_chr", i, "_b", bb # ,"_agesx"
        , ".txt",
        sep = ""
      )
      cat("Writing the file (quick)\n")
      names(pv) <- snps
      write.table(pv, ff, quote = F, col.names = F)
      # I write SNP names; also NA p-values
    }
  }
}
