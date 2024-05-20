# Here we define function "pow," which we will use in every simulation

#### Auxiliary function for time management ####
hms <- function(s) {
  s <- round(s)
  if (s <= 60) {
    return(paste(s, "s"))
  } else if (s <= 3600) {
    return(paste(s %/% 60, "min", s %% 60, "s"))
  } else {
    h <- s %/% 3600
    s <- s %% 3600
    m <- s %/% 60
    s <- s %% 60
    if (m != 0) {
      return(paste(h, "h", m, "min", s %% 60, "s"))
    } else {
      return(paste(h, "h", s %% 60, "s"))
    }
  }
}

source("../sim_functions_snp_snp.R") # Auxiliary functions for testing with our method
source("aux_functions_boost.R") # Auxiliary functions for testing with BOOST

#### FUNCTION FOR SIMULATION ####
pow <- function(repl = 1E4, nbasal = 500, modelX = "qexp", expoX = 1, gX = 1,
                modelY = "qexp", expoY = 1, gY = 1, alpha.nom = .05,
                N = 1E3, ranseed = 1, determ = T, silent = T, append = T,
                file = "pow.dat", FWER = "No correct",
                method = "dcov") {
  if (!determ) {
    ranseed <- NA
  } else {
    set.seed(ranseed)
  }

  if (expoX < 1 | expoY < 1) {
    cat("expo must be >=1\n")
    return(NA)
  }
  if (gX < 0 | gY < 0) {
    cat("g must be >=0\n")
    return(NA)
  }
  if (gX > 1 | gY > 1) {
    cat("g must be <=1\n")
    return(NA)
  }

  #### Initial parameters of the hypothesis test ####
  # FWER correction:
  if (FWER == "Sidak") {
    alpha <- 1 - sqrt(1 - alpha.nom) # Sidak
  } else if (FWER == "Bonferroni") {
    alpha <- alpha.nom / 2 # Bonferroni
  } else if (FWER == "No correct") {
    alpha <- alpha.nom
  }


  #### Model generation ####
  if (!silent) cat("Generating model...\n")
  n1 <- nbasal + sample(-10:10, 1)
  n2 <- nbasal + sample(-10:10, 1)
  # Sample sizes with a small perturbation,
  # so that they are not identical

  nrej <- 0

  if (!silent) {
    cat("Starting the replicate loop.\n")
    t0 <- Sys.time()
  }

  for (w in 1:repl) { # Replicates loop

    Xmat <- matrix(NA, n1, 2) # Column 1: X_i. Column 2: X_j
    Ymat <- matrix(NA, n2, 2) # Column 1: Y_i. Column 2: Y_j

    #### Contingency table for cases (X) ####
    tableX <- matrix(NA, 3, 3)
    theta1 <- sample((5:20) / 100, 1)
    pr <- (1 - theta1)^2
    q <- 2 * theta1 * (1 - theta1) # HWE
    theta2 <- sample((5:20) / 100, 1)
    r <- (1 - theta2)^2
    s <- 2 * theta2 * (1 - theta2) # HWE
    tableX <- c(pr, q, 1 - pr - q) %o% c(r, s, 1 - r - s)
    if (modelX == "indep") {
      # Nothing else to do, in this case
    } else if (modelX == "qexp") {
      tableX[1, 1] <- pr * r - q * s + q^expoX * s
      tableX[1, 2] <- pr * s - q^expoX * s + q * s
      tableX[2, 1] <- q * r - q^expoX * s + q * s
      tableX[2, 2] <- q^expoX * s
    } else if (modelX == "rexp") {
      tableX[2, 2] <- q * s + (1 - pr - q)^expoX * (1 - r - s) - (1 - pr - q) * (1 - r - s)
      tableX[2, 3] <- q * (1 - r - s) + (1 - pr - q) * (1 - r - s) - (1 - pr - q)^expoX * (1 - r - s)
      tableX[3, 2] <- (1 - pr - q) * s + (1 - pr - q) * (1 - r - s) - (1 - pr - q)^expoX * (1 - r - s)
      tableX[3, 3] <- (1 - pr - q)^expoX * (1 - r - s)
    } else if (modelX == "qmult") {
      tableX[1, 1] <- pr * r - q * s + gX * q * s
      tableX[1, 2] <- pr * s - gX * q * s + q * s
      tableX[2, 1] <- q * r - gX * q * s + q * s
      tableX[2, 2] <- gX * q * s
    } else if (modelX == "rmult") {
      tableX[2, 2] <- q * s + gX * (1 - pr - q) * (1 - r - s) - (1 - pr - q) * (1 - r - s)
      tableX[2, 3] <- q * (1 - r - s) + (1 - pr - q) * (1 - r - s) - gX * (1 - pr - q) * (1 - r - s)
      tableX[3, 2] <- (1 - pr - q) * s + (1 - pr - q) * (1 - r - s) - gX * (1 - pr - q) * (1 - r - s)
      tableX[3, 3] <- gX * (1 - pr - q) * (1 - r - s)
    }

    #### Filling up the matrix for cases (X) ####
    lett <- sample(LETTERS[1:9], n1, T, as.numeric(tableX))
    Xmat[lett == "A", ] <- 0
    Xmat[lett == "B", 1] <- 1
    Xmat[lett == "B", 2] <- 0
    Xmat[lett == "C", 1] <- 2
    Xmat[lett == "C", 2] <- 0
    Xmat[lett == "D", 1] <- 0
    Xmat[lett == "D", 2] <- 1
    Xmat[lett == "E", ] <- 1
    Xmat[lett == "F", 1] <- 2
    Xmat[lett == "F", 2] <- 1
    Xmat[lett == "G", 1] <- 0
    Xmat[lett == "G", 2] <- 2
    Xmat[lett == "H", 1] <- 1
    Xmat[lett == "H", 2] <- 2
    Xmat[lett == "I", ] <- 2


    #### Contingency table for controls (Y) ####
    # I assume marginal frequencies don't differ with respect to X.
    tableY <- c(pr, q, 1 - pr - q) %o% c(r, s, 1 - r - s)
    if (modelY == "indep") {
      # Nothing else to do, in this case
    } else if (modelY == "qexp") {
      tableY[1, 1] <- pr * r - q * s + q^expoY * s
      tableY[1, 2] <- pr * s - q^expoY * s + q * s
      tableY[2, 1] <- q * r - q^expoY * s + q * s
      tableY[2, 2] <- q^expoY * s
    } else if (modelY == "rexp") {
      tableY[2, 2] <- q * s + (1 - pr - q)^expoY * (1 - r - s) - (1 - pr - q) * (1 - r - s)
      tableY[2, 3] <- q * (1 - r - s) + (1 - pr - q) * (1 - r - s) - (1 - pr - q)^expoY * (1 - r - s)
      tableY[3, 2] <- (1 - pr - q) * s + (1 - pr - q) * (1 - r - s) - (1 - pr - q)^expoY * (1 - r - s)
      tableY[3, 3] <- (1 - pr - q)^expoY * (1 - r - s)
    } else if (modelY == "qmult") {
      tableY[1, 1] <- pr * r - q * s + gY * q * s
      tableY[1, 2] <- pr * s - gY * q * s + q * s
      tableY[2, 1] <- q * r - gY * q * s + q * s
      tableY[2, 2] <- gY * q * s
    } else if (modelY == "rmult") {
      tableY[2, 2] <- q * s + gY * (1 - pr - q) * (1 - r - s) - (1 - pr - q) * (1 - r - s)
      tableY[2, 3] <- q * (1 - r - s) + (1 - pr - q) * (1 - r - s) - gY * (1 - pr - q) * (1 - r - s)
      tableY[3, 2] <- (1 - pr - q) * s + (1 - pr - q) * (1 - r - s) - gY * (1 - pr - q) * (1 - r - s)
      tableY[3, 3] <- gY * (1 - pr - q) * (1 - r - s)
    }

    #### Filling up the matrix for controls (Y) ####
    lett <- sample(LETTERS[1:9], n1, T, as.numeric(tableY))
    Ymat[lett == "A", ] <- 0
    Ymat[lett == "B", 1] <- 1
    Ymat[lett == "B", 2] <- 0
    Ymat[lett == "C", 1] <- 2
    Ymat[lett == "C", 2] <- 0
    Ymat[lett == "D", 1] <- 0
    Ymat[lett == "D", 2] <- 1
    Ymat[lett == "E", ] <- 1
    Ymat[lett == "F", 1] <- 2
    Ymat[lett == "F", 2] <- 1
    Ymat[lett == "G", 1] <- 0
    Ymat[lett == "G", 2] <- 2
    Ymat[lett == "H", 1] <- 1
    Ymat[lett == "H", 2] <- 2
    Ymat[lett == "I", ] <- 2

    if (method == "dcov") {
      #### Test for Xmat (cases) ####
      pvX <- test1(Xmat)
      if (unique(dim(pvX)) != 2) cat("Trouble...\n")
      RejectX <- pvX[1, 2] < alpha

      #### Test for Ymat (controls) ####
      pvY <- test1(Ymat)
      if (unique(dim(pvY)) != 2) cat("Trouble...\n")
      RejectY <- pvY[1, 2] < alpha

      #### Second stage of the testing procedure ####
      Reject <- RejectX & !RejectY
    } else if (method == "boost") {
      pv <- test_boost(Xmat, Ymat, alpha = alpha.nom)
      Reject <- pv < alpha.nom
    }
    nrej <- nrej + Reject

    #### Estimation of the remaining time ####
    if (!silent) {
      t1 <- Sys.time()
      Deltat <- as.numeric(difftime(t1, t0, units = "secs"))
      prop <- w / repl
      cat(
        100 * prop, "% complete.", "Time spent:", hms(Deltat),
        ". Remaining time:", hms(Deltat / prop * (1 - prop)), "\n"
      )
    }
  } # End of loop w=1:repl


  #### RESULTS ANALYSIS ####
  beta <- nrej / repl

  #### OUTPUT TO file ####
  paramX <- NA -> paramY # Prevents an ancestral bug.
  # Columns are: modelX, paramX, modelY, paramY, alpha.nom, beta
  if (modelX %in% c("qexp", "rexp")) paramX <- expoX
  if (modelX %in% c("qmult", "rmult")) paramX <- gX
  if (modelY %in% c("qexp", "rexp")) paramY <- expoY
  if (modelY %in% c("qmult", "rmult")) paramY <- gY
  write.table(t(c(modelX, paramX, modelY, paramY, alpha.nom, beta)), file,
    row.names = F, col.names = F, quote = F, eol = "\r\n", append = append
  )
  cat("The results file is called", file, "and can be found in:\n")
  print(getwd())
} # End of function "pow"
