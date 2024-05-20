p <- 6456
v <- 1:p
psums <- v * (v - 1) / 2
rm(v)

itoij <- function(k, p = 6456, ps = psums) {
  # Coordinates of vector to upper tri matrix (empty diagonal!)
  if (k == 1) {
    return(c(1, 2))
  }
  j <- max(which(ps < k)) + 1
  i <- k - ps[j - 1]
  return(c(i, j))
}

ijtoi <- function(v) {
  # The inverse function of itoij
  i <- v[1]
  j <- v[2]
  if (length(v) != 2) {
    cat("No. of indices != 2.\n")
    return(NA)
  }
  if (i > j) {
    cat("Out of triu.\n")
    return(NA)
  }
  return(i + (j - 1) * (j - 2) / 2)
}
