# Create a RLE plot customized to time-series data

`create_RLEplot` generates a RLE plot with samples colored according to
subjects.

## Usage

``` r
create_RLEplot(exp_matrix, ylim)
```

## Arguments

- exp_matrix:

  Matrix of the raw expression of genes or log2 transformed normalized
  genes.

- ylim:

  Numeric vector of length 2 for the y-coordinate range.

## Value

A RLE plot

## Details

Only a raw expression matrix of genes or a log2 transformed matrix of
normalized genes with properly formatted column names is needed to
create a plot.

## Author

Jayanth Sreekanth

## Examples

``` r

# Path of the sample RCC files
path_to_rcc <- system.file("extdata", "GSE146530", package = "rhythmCI")
meta <- system.file("extdata", "meta_GSE146530.csv", package = "rhythmCI")

# Load and No normalization -----------------------------------------------------------------
suppressWarnings(

exp_df1 <- loadnormalizercc(dir = path_to_rcc, meta_dir = meta,
                            hkgenes = c("LOC106606270", "LOC106600007"),
                            file_number = FALSE
                            )
)

df1_raw <- exp_df1$rawdata
df1_pdata <- exp_df1$pData
df1_fdata <- exp_df1$fData
df1_endblod <- exp_df1$num_endogenous_blod

# Perform nSolver normalization ------------------------------------------------------------
suppressMessages(

exp_df3 <- normalizeruvg(rawdata = df1_raw,
  fData = df1_fdata,
  pData = df1_pdata,
  hkgenes = c("LOC106606270", "LOC106600007"),
  exclude = NULL,
  k=3)
)

df3_set <- exp_df3$set
df3_vsd <- exp_df3$vsd
df3_ruvgnormlog2 <- exp_df3$ruvgnormlog2

create_RLEplot(df3_ruvgnormlog2, ylim = c(-2,2))

```
