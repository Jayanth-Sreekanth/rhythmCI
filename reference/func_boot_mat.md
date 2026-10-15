# Matrix-wise operation - sub function to estimate the amplitude during every bootstrap instance

`func_boot_mat` adds the sampled residuals with the obtained best fit.
Then performs harmonic regression and returns the amplitude for each
bootstrap instance.

## Usage

``` r
func_boot_mat(df_resi, n, lm_fitted, tau)
```

## Arguments

- df_resi:

  The sampled matrix of residuals in an instance

- n:

  Row number in the matrix

- lm_fitted:

  The fitted value of initial harmonic regression

- tau:

  The period of one oscillation

## Value

An amplitude estimate for the instance

## Author

Jayanth Sreekanth
