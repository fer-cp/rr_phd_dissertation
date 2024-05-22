# Reproducing the real data analysis

The scripts in this folder correspond to the example of hepatic enzymes we study in the main body of the dissertation.
This data is available through dbGaP to anyone who fulfills their strict requirements on information security, and the agreement
we have signed does not allow for sharing the data with third parties. Therefore, the software we here share can potentially be used
with that or other similar data, but we do not provide any specific files for it. Once more, it is an option to use the toy example that
comes with every release of PLINK (the genetic software that we again use in this application).

The individual numerical results for each SNP are obtained by running the script _enzymes.R_, which again depends on R functions that call Python
for the computation of _p_-values. Once this has been run, Manhattan plots can be created by means of _manh_plots.R_.

R packages used:
* coga (for the generalised _F_ distribution);
* qqman (Manhattan plots);
* data.table (fast and efficient reading and writing of external files).


