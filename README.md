# Repro repo for: _Nonparametric Independence Tests in High-Dimensional Settings, with Applications to the Genetics of Complex Disease_

Resproducible research material for the numerical examples of the PhD dissertation with title:

_Nonparametric Independence Tests in High-Dimensional Settings, with Applications to the Genetics of Complex Disease_

# System requirements

Software in this folder mostly relies on [R](https://cran.r-project.org) and [R packages](https://cran.r-project.org/web/packages/available_packages_by_name.html) developed by various authors. We recommend R version 4.3.1+, in Windows 10+, for running our scripts. Details on the OS used, as well as the version of each package, are provided below.

The applications to genetics depend on [PLINK v1.9](https://www.cog-genomics.org/plink/1.9/), by summoning `plink.exe` from the R scripts. The \*.exe file is expected to have been downloaded from the PLINK website into the current working directory of R. Users of OSs other than Windows should adapt the command line for calling PLINK to the requirements of their system.


# General information on the R version used
After running the `sessionInfo()` command, we get the information on our R session that we show below. We note that we do not expect any important deviation in the functioning of our scripts with incremental updates to R and its packages in the near future, but it is a good practice to document the session information regardless.
```
R version 4.3.1 (2023-06-16 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)
Matrix products: default
locale:
[1] LC_COLLATE=Spanish_Spain.utf8 LC_CTYPE=Spanish_Spain.utf8
[3] LC_MONETARY=Spanish_Spain.utf8 LC_NUMERIC=C
[5] LC_TIME=Spanish_Spain.utf8
time zone: Europe/Madrid
tzcode source: internal
attached base packages:
[1] stats graphics grDevices utils datasets methods base
other attached packages:
[1] ggplot2_3.4.3 momentchi2_0.1.5 coga_1.2.0
loaded via a namespace (and not attached):
[1] utf8_1.2.3 R6_2.5.1 magrittr_2.0.3 gtable_0.3.4 glue_1.6.2
[6] tibble_3.2.1 pkgconfig_2.0.3 lifecycle_1.0.3 cli_3.6.1 fansi_1.0.5
[11] scales_1.2.1 grid_4.3.1 vctrs_0.6.3 withr_2.5.1 compiler_4.3.1
[16] tools_4.3.1 pillar_1.9.0 munsell_0.5.0 Rcpp_1.0.11 colorspace_2.1-0
[21] rlang_1.1.1
```
