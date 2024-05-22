# Reproducing the simulation study

The numerical results for simulations of type I error and power can be re-run by sourcing the R files with self-explanatory names in this folder (of the form _typeI*.R_ and _powersimu*.R_). The computation times are estimated in the script _comptime.R_.

All the graphics in the main body of the dissertation can be reproduced by first running the numerical results and then using the plotting configuration in _plots.R_ .

There is a script with testing functions, which is used every time that a numerical result for our methodology is generated. It is called _sim_functions_snp_pheno.R_ and it calls a Python file named _pvalue_python.py_ for the evaluation of _p_-values with the library _mpmath_.

The following R packages are used:
* AssocTests (for comparing with pre-existing competing tests);
* parallel (allows for multi-thread or multi-core computations, whenever the hardware meets these needs);
* microbenchmark (measuring times).
