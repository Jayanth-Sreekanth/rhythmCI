# 

     Quality control - limit of detection 

nSolver limit of detection flag identified in the nSolver manual Scripts
copied from
https://github.com/bhattacharya-a-bt/NanoNormIter/tree/master/R Citing
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8138885/

## Usage

``` r
limitOfDetectionQC(rcc, numSD)
```

## Arguments

- rcc:

  loaded and valid rcc file

- numSD:

  number of standard deviations to calibrate the LOD

## Value

A character of 'Flag' or 'No flag' for limit of detection
