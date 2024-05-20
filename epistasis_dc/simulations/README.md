# Reproducing our simulation study

By simply running the R script [masterscript_power.R],
one generates the data tables (as _*.dat_ files) necessary to produce the plots
that we display as a result of our simulation study (Section 4 of the preprint; Section 3.5 in the dissertation).

The power plots (which include the comparison with competing method BOOST, in a different colour)
are directly generated when running the masterscript.

The code for the calibration plots is a bit more cumbersome, due to the confidence band,
so we split it to a separate script. Please run [plotting_calibration.R] to obtain the corresponding plots
(Figures 2 of the paper; Figure 3.2 in the dissertation).

In order to obtain plots or results for other models, one should either perform small manual edits
in the scripts, or run the simulation functions with different values of the parameters.
