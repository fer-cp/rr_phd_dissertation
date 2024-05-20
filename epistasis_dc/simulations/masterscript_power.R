# This file generates the data tables for the calibration plots,
# as well and the tables and plots for the power comparison with BOOST.

source("power_simulation_functions.R")

#### Running "pow" for alpha-calibration ####
for (a in c(1, 2, 5, 10) / 100) {
  pow( repl = 1E3, modelX = "indep", modelY = "indep",
       file = "Calibr_indep_alpha.dat", alpha.nom = a  )
}

for (a in c(1, 2, 5, 10) / 100) {
  pow( repl = 1E3, modelX = "rexp", expoX = 10, modelY = "rexp", expoY = 10,
       file = "Calibr_rexp_e10_alpha.dat", alpha.nom = a )
}

#### Power estimation and plotting for indep/qexp ####
finame <- "Power_indep_qexp"
for (e in c(1 + (0:10) / 10, 2:10)) {
  pow( modelX = "qexp", modelY = "indep", expoX = e,
       file = paste(finame, ".dat", sep = "")  )
}
## Plotting:
resu <- read.table(paste(finame, ".dat", sep = ""))
pdf(paste(finame, ".pdf", sep = ""))
xx <- resu[, 2]
yy <- resu[, 6]
plot(xx, yy,
  type = "b", col = "blue", axes = F,
  xlab = "e",
  ylab = "Empirical power", ylim = c(0, 1)
)
axis(2, seq(0, 1, .1))
axis(1, c(1:10))
dev.off()


#### Power estimation and plotting for indep/qmult ####
finame <- "Power_indep_qmult"
for (g in (1:10) / 10) {
  pow( modelX = "qmult", modelY = "indep", gX = g,
       file = paste(finame, ".dat", sep = "")  )
}

#### SAME, BUT FOR BOOST ####
#### Running "pow" for alpha-calibration ####
for (a in c(1, 2, 5, 10) / 100) {
  pow( repl = 1E3, modelX = "indep", modelY = "indep", silent = F,
       file = "Calibr_indep_alpha_BOOST.dat", alpha.nom = a,
       method = "boost"  )
}

for (a in c(1, 2, 5, 10) / 100) {
  pow( repl = 1E3, modelX = "rexp", expoX = 10, modelY = "rexp", expoY = 10,
       file = "Calibr_rexp_e10_alpha_BOOST.dat", alpha.nom = a,
       silent = F, method = "boost"  )
}

## To plot calibration, go to: plotting_calibration.R

#### Power estimation and plotting for indep/qexp ####
finame <- "Power_indep_qexp_BOOST"
for (e in c(1:10, 1 + (1:9) / 10)) {
  pow( repl = 1E3, modelX = "qexp", modelY = "indep", expoX = e,
       file = paste(finame, ".dat", sep = ""), silent = F, method = "boost"  )
}
## Plotting:
resu <- read.table(paste(finame, ".dat", sep = ""))
pdf(paste(finame, ".pdf", sep = ""))
xx <- resu[, 2]
yy <- resu[, 6]
yy <- yy[order(xx)]
xx <- sort(xx)
plot(xx, yy,
  type = "b", col = "red", axes = F,
  xlab = "e",
  ylab = "Empirical power", ylim = c(0, 1)
)
resu <- read.table(paste(woboost(finame), ".dat", sep = ""))
xx <- resu[, 2]
yy <- resu[, 6]
yy <- yy[order(xx)]
xx <- sort(xx)
lines(xx, yy,
  type = "b", col = "blue",
  xlab = "e",
  ylab = "Empirical power", ylim = c(0, 1)
)
axis(2, seq(0, 1, .1))
axis(1, c(1:10))
dev.off()


#### Power estimation and plotting for qmult/qexp ####
finame <- "Power_indep_qmult_BOOST"
for (g in (1:10) / 10) {
  pow( repl = 1E3, modelX = "qmult", modelY = "indep", gX = g,
       file = paste(finame, ".dat", sep = ""), silent = F, method = "boost" )
}
## Plotting:
resu <- read.table(paste(finame, ".dat", sep = ""))
pdf(paste(finame, ".pdf", sep = ""))
xx <- resu[, 2]
yy <- resu[, 6]
plot(xx, yy,
  type = "b", col = "red", axes = F,
  xlab = "g",
  ylab = "Empirical power", ylim = c(0, 1)
)
resu <- read.table(paste(woboost(finame), ".dat", sep = ""))
xx <- resu[, 2]
yy <- resu[, 6]
yy <- yy[order(xx)]
xx <- sort(xx)
lines(xx, yy,
  type = "b", col = "blue",
  xlab = "e",
  ylab = "Empirical power", ylim = c(0, 1)
)
axis(2, seq(0, 1, .1)) # Y-axis
axis(1, (1:10) / 10) # X-axis
dev.off()
