lettrify <- function(v) { # Generates an example genotype from 0-1-2 data
  v <- replace(v, v == 0, "A A")
  v <- replace(v, v == 1, "A T")
  v <- replace(v, v == 2, "T T")
  return(v)
}
make_ped <- function(X, Y, finame = "fi/current_boost_iteration",
                     silent = T) {
  # Creates a PED file that PLINK's implementation of BOOST can read
  n1 <- nrow(X)
  n2 <- nrow(Y)
  n <- n1 + n2
  if (ncol(X) != 2) cat("Trouble...\n")
  if (ncol(Y) != 2) cat("Trouble...\n")
  ped <- matrix(NA, n, 8)
  ped[, 1] <- 1:n # family
  ped[, 2] <- 1:n # individual
  ped[, 3:4] <- 0 # maternal, paternal
  ped[, 5] <- sample(1:2, n, replace = T) # sex
  ped[, 6] <- c(rep(2, n1), rep(1, n2)) # affection status
  ped[, 7] <- lettrify(c(X[, 1], Y[, 1]))
  ped[, 8] <- lettrify(c(X[, 2], Y[, 2]))
  fiped <- paste(finame, ".ped", sep = "")
  write.table(ped, fiped, quote = F, col.names = F, row.names = F)
  if (!silent) cat("I wrote file", fiped, "\n")
  #
  map <- matrix(NA, 2, 3)
  map[, 1] <- 22 # chr
  f <- function(x) paste("rs000", x, sep = "")
  map[, 2] <- sapply(1:2, f) # SNP identifiers
  g <- function(x) paste("12345", x, sep = "")
  bp <- sapply(1:2, g)
  h <- function(x) paste("0", x, collapse = " ")
  map[, 3] <- sapply(bp, h) # positions
  fimap <- paste(finame, ".map", sep = "")
  write.table(map, fimap, quote = F, col.names = F, row.names = F)
  if (!silent) cat("I wrote file", fimap, "\n")
}

test_boost <- function(X, Y, finame = "fi/current_boost_iteration",
                       silent = T, alpha = .05) {
  # Parser for testing one single SNP pair with PLINK's implementation of BOOST
  make_ped(X, Y, finame = finame, silent = silent)
  system2("plink", c(
    paste("--file", finame, collapse = " "),
    "--make-bed",
    paste("--out", finame, collapse = " ")
  ), stdout = !silent)
  system2("plink", c(
    paste("--bfile", finame, collapse = " "),
    "--fast-epistasis boost",
    paste("--epi1", alpha, sep = " "),
    paste("--out", finame, collapse = " ")
  ), stdout = !silent)
  D <- read.table(paste(finame, ".epi.cc", sep = ""), header = T)
  pval <- NA
  if (nrow(D) >= 1) {
    pval <- D$P
    if (pval < 0 | pval > 1) cat("Trouble...\n")
  }
  if (nrow(D) == 0) pval <- 1
  return(pval)
}
