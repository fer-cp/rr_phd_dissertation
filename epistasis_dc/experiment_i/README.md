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
