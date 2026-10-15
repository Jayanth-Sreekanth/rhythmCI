# getSpikeInInput

Gets the RNA “spike-in” input levels for positive and negative control
probes from the label in their GeneName. Note that this is a helper
function for readRlf() and elsewhere and is not intended for external
use.

## Usage

``` r
getSpikeInInput(CodeClass, GeneName, nonCtrlProbeVal = NA)
```

## Arguments

- CodeClass:

  Character vector with code classes for each probe.

- GeneName:

  Character vector with gene names/symbols for each probe.

- nonCtrlProbeVal:

  Value to assign as the spike-in input for the non-control probes.

## Value

A data frame with the input CodeClass and GeneName but where the latter
has been split into two columns: one showing the GeneName for each probe
with spike-in input labels removed – and another with the spike-in input
levels.

## Author

NanoStringQCPro
