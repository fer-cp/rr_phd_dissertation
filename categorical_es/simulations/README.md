# Reproducing the simulation study

This chapter of the dissertation proposes testing procedures for two separate problemas with categorical data: independence of two variables, and goodness of fit of one variable to a given distribution.

There is an R script called _test_functions_ct_dcov.R_ which provides the testing functions necessary for both problems, and then the simulations for each of the two is organised in different folders.

For independence, the numerical results are generated in _simu_indep_with_plots.R_ , which also provides plots for the power curve comparison of our methodology with competitors. The figures related to the type I error control can be generated with _plotting_calibration_methods_indep.R_ .

For goodness of fit, the numerical results are generated in _simu_gof.R_ . Power plots are created by sourcing _plotting_power_gof.R_. The figures related to the type I error control can be generated with _plotting_calibration_gof.R_ .

R packages used:
* CompQuadForm (evaluation of the distribution function of quadratic forms of Gaussian variables);
* ggplot2 (graphic functions that expand those in R by default).

