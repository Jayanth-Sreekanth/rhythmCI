# Matrix-wise operation - main function

`rhy_func_mat` performs harmonic regression and returns the CI of
amplitude, rhythmic status and mesor. Note that the mesor is estimated
before bootstrapping

## Usage

``` r
rhy_func_mat(df, ttl_sub, tau, alpha, fc_cutoff, num_boot)
```

## Arguments

- df:

  Matrix of value, columns corresponding to variables and rows
  corresponding to samples w.r.t time

- ttl_sub:

  Total number of subjects

- tau:

  The period of one oscillation

- alpha:

  The significance level to estimate the confidence level. Default is
  0.05

- fc_cutoff:

  The amplitude cutoff to define the rhythmic status of variables (genes
  and other)

- num_boot:

  The number of bootstrap samples to be generated

## Value

A list with the CI of amplitude, rhythmic status and mesor for each
variable separated according to subjects

## Author

Jayanth Sreekanth
