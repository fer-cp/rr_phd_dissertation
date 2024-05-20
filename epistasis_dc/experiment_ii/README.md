# Reproducing Experiment II
As in Experiment I, we begin with the full GWAS database. We assume that we have a triplet of
PLINK files (_*.bed, *.bim, *.fam_) within the experiment_ii folder. We are not allowed to share ours,
but one can run the script filtering_snps_experiment_ii.R to obtain the matrices with the observations
for cases and controls with the same filters that we describe in the supplement. In our case, this
yields L=6456 SNPs for Experiment II. We indicate in the comments of the _*.R_ file the results of the
relevant steps.

An important caveat is that we do not attach the GTEx files necessary for running this script. For
obtaining them, one should visit the GTEx Portal on https://www.gtexportal.org/home/datasets ,
download the file GTEx_Analysis_v7_eQTL.tar.gz (single tissue cis-eQTL data for GTEx Analysis
V7, dbGAP accession phs000424.v7.p2 ), unzip it and place all the _Brain\_*.signifpairs.txt_ files
within the folder [experiment_ii/gtex_v7_signifpairs/brain] , and all the remaining
_\*.signifpairs.txt_ files (the ones not beginning with _Brain\_\*_) in the analogous nonbrain folder.

To replicate Experiment II, the reader is kindly asked to run masterscript_experiment_ii.R. All the
observations we made for the masterscript of Experiment I also apply to this one.
Source of data.


Due to confidentiality, we are not allowed to share the whole database, so we will directly start
the reproducibility scripts with the two subsets that are used on Experiment I (8030 SNPs) and
Experiment II (6456) SNPs. We do provide scripts with the necessary code for filtering GWAS
datasets in the same way we did in each experiment (file names
filtering_snps_experiment_i[i].R)
Please note that some SNP IDs have been altered for the purpose of sharing our code for
reproducibility, in a way
