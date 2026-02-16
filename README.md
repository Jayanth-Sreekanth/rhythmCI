
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rhythmCI

<!-- badges: start -->

<!-- badges: end -->

A toolkit covering complete circadian data analysis workflow primarily
for Nanostring nCounter data, encompassing everything from importing and
processing to normalization, analysis and visualization. The regression
based analysis for estimating the confidence interval of amplitude is
broadly applicable, not only for gene-expression data, but also for
other biological data. An informative visualization of results accompany
the analysis.

## Installation

You can install the development version of rhythmCI from
[GitHub](https://github.com/Jayanth-Sreekanth/rhythmCI) with:

``` r
# install.packages("pak")
pak::pak("Jayanth-Sreekanth/rhythmCI", build = TRUE)
```

Alternatively, installation can be done with:

``` r
# install.packages("remotes")
remotes::install_github("Jayanth-Sreekanth/rhythmCI", build_vignettes = TRUE)
```

## Vignettes and Documentation

For detailed examples and workflows, refer to the vignettes:

``` r
# Open the main vignette
vignette("rhythmCI", package="rhythmCI")
```
