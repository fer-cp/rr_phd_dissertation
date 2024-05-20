## Chromosome and position:
chr_pos <- as.matrix(read.table("chr_pos.dat"))

## Functions for testing "nearness":
samechr <- function(h, cp = chr_pos) { # Do the 2 SNPs lay on the same chromosome?
  if (length(h) != 2) {
    return(NA)
    cat("Trouble...\n")
  }
  return(length(unique(cp[h, 1])) == 1)
}
nearby <- function(h, cp = chr_pos) { # Are the two SNPs < 1 Mb apart?
  if (length(h) != 2) {
    return(NA)
    cat("Trouble...\n")
  }
  if (samechr(h) & abs(cp[h[1], 2] - cp[h[2], 2]) < 1E6) {
    return(T)
  }
  return(F)
}
