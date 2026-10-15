# Column-wise operation - sub function

`fit_fun` performs harmonic regression, calls bootstrapping, calculates
the amplitude CI

## Usage

``` r
fit_fun(nm, df_dyn, n_col, ttl_sub, tau, alpha, fc_cutoff, num_boot)
```

## Arguments

- nm:

  Column number

- df_dyn:

  The matrix with columns having more than three time-points

- n_col:

  Number of columns in the matrix

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
