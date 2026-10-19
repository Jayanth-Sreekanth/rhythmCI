# Rhythmicity analysis using Harmonic Regression and bootstrapping

`rhythmicity` applies harmonic regression to the variables and
calculates confidence intervals using bootstrapping. The intervals are
corrected using bonferroni method in the end.

## Usage

``` r
rhythmicity(
  input,
  operation,
  seed = NA,
  tau = 24,
  num_boot,
  alpha = 0.05,
  fc_cutoff
)
```

## Arguments

- input:

  Matrix, where columns correspond to variables and rows correspond to
  samples w.r.t time.

- operation:

  Anyone of the following should be provided: "matrix" or "column". Do
  not use "matrix" mode when there are missing values in the matrix,
  instead use "column" mode.

- seed:

  Is used to reproduce the same results across multiple runs.

- tau:

  The period of oscillation.

- num_boot:

  The number of bootstrap samples to be generated.

- alpha:

  The significance level to estimate the confidence level. Default is
  0.05.

- fc_cutoff:

  The amplitude cutoff to define the rhythmic status of variables (genes
  and other).

## Value

A list for each subject containing the following: a lower (ci.lower) and
upper (ci.upper) value of the confidence interval of the amplitude, the
rhythmic status of each variable (rhythmicity), mesor from harmonic
regression (hr.mesor), amplitude from harmonic regression
(hr.amplitude), and phase from harmonic regression (hr.phase).

## Details

A matrix with samples across time in rows and with columns of variables
are to be provided to perform harmonic regression. After regressing the
data with cosine functions, the best fit and the residuals are obtained.
The residuals are then sampled with replacement using bootstrapping. In
each instance of bootstrapping, the sampled residuals are summed to the
best fit yielding a bootstrap estimate of amplitude. The Confidence
Interval (CI) for amplitude estimates will be computed using the
provided alpha value (also includes Bonferroni correction). Using the CI
and fc_cutoff the rhythmic status of each variables are defined as
follows: rhythmic (if the lower and upper CI lies above the cutoff),
arrhythmic (if the lower and upper CI lies below or is equal to the
cutoff) and inconclusive (if the lower and upper CI straddles the
cutoff). The obtained mesor (the average value), amplitude, and phase is
estimated before bootstrapping.

Please note that at least 4 samples are needed per subject. Only the
variables with at least 4 samples will be analyzed. The column-wise
operation will need more time than matrix-wise operation.

## Author

Jayanth Sreekanth, Bharath Ananthasubramaniam

## Examples

``` r

# Load the data that comes along with the package and run rhythmicity.
# This data was generated using simphony.

data(simdf)

suppressWarnings(

output_hr <- rhythmicity(input = simdf,
                         operation = "matrix",
                         seed = 7,
                         tau = 24, #value in general for humans
                         num_boot = 999,
                         alpha = 0.05,
                         fc_cutoff = 0.5)

)
```
