# Reproducing the real data application

As with the simulations, here we do everything twice, once for the independence test and another time for that of goodness of fit.

Unlike in other chapters of the dissertation, the real data applications that we present here do not involve individual-level genotype data, so there is no privacy concern. Therefore, the real data examples here can be run fully.

For independence, we provide the data for the example on admission history of schizophrenia patients in _admission_data.txt_ . It is then analysed in _admission.R_ . Relevant results and intermediate steps are marked as comments in that script.

For goodness of fit, the data of the allelic frequencies is typed out inside the corresponding R scripts, and the external data file _pgc3_snps.txt_ provides a list of SNPs known to be associated with schizophrenia, against which we check the variants we consider in each example. The testing for HWE in a biallelic locus is carried out in _hwe_2allele.R_ , whereas the triallelic setting is dealt with in _hwe_3allele.R_ .
