# 

     Read an .RCC file 

Parse an .RCC file into a list with each part of the file (Header,
Sample_Attributes, Lane_Attributes, Code_Summary, etc) stored as a
vector or data frame.

## Usage

``` r
readRcc(rcc, removeSpikeInLabels = TRUE)
```

## Arguments

- rcc:

  Path to the .RCC file.

- removeSpikeInLabels:

  Logical. If TRUE (the default), RNA “spike-in” input labels (if any)
  in the GeneName for positive and negative control probes will be
  removed.

## Value

A list where each element holds the contents of one part of the .RCC
file (Header, Sample_Attributes, Lane_Attributes, Code_Summary, etc) as
a vector or data frame.

## Author

NanoStringQCPro
