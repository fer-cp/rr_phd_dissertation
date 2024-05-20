# Reproducing Experiment I
We use our full schizophrenia GWAS database,
as described in Appendix F of the supplement of the preprint,
and Section 3.6.1 of the dissertation.

We assume that we have a triplet of PLINK files (\*.bed, \*.bim, \*.fam) within the `experiment_i` folder,
which must be set as our current working directory.

We are not allowed to share our original PLINK files, but one can run the script [_filtering_snps_experiment_i.R_](https://github.com/fer-cp/rr_phd_dissertation/blob/main/epistasis_dc/experiment_i/filtering_snps_experiment_i.R) to obtain
the matrices with the observations for cases and controls with the same filters that we describe in the
supplement. One can do so, for example, with the `toy.ped` example supplied at the [PLINK demo](https://www.cog-genomics.org/plink/1.9/).

In our case, the filters yield L=8030 SNPs for Experiment I. We indicate as comments of the
*.R file the results of each relevant step.

One should run [masterscript_experiment_i.R] to reproduce Experiment I. It uses as input the data
from the 8030 SNPs for cases and controls (Matrix_X.dat and Matrix_Y.dat), as well as the SNP IDs
in the chromosome-position format (with some alterations for the sake of anonymity of sampled
individuals, in order to make this data shareable; such modifications do not influence the results we
present). The latter can be found in chr_pos.dat. Every relevant result has been written down as a
comment in the *.R file. We recommend using the search function with the query “result present in
the manuscript” to find the exact lines of code that replicate every numerical result for Experiment I
that is cited in the main manuscript.

At many points of the script, we generate intermediate result files, in order to ease running only parts
of it. We do this in light of the moderately long running times of some segments, but it is also
feasible to run the entire script within reasonable time in any modern desktop computer.
Please note that it is necessary to set as the working directory the location of the masterscript R file
before running it.
