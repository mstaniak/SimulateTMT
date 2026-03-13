# SimulateTMT: R package for simulating feature-level mass spectrometry data 

## Overview of the package

`SimulateTMT` package provides functions that enable simulating feature-level data from mass spectrometry-based proteomics experiment. 
Output follows the `MSstats` format and can be used to evaluate protein quantification methods. Observations corresponding to both unique and shared peptides can be created,
and all observations are perturbation of user-defined protein-level summaries. Hence, both estimation of protein-level abundances and of log-fold changes can be evaluated. 
Various experimental designs and structures of peptide-protein networks can be simulated. 

## Installation:

```
if (!require(devtools)) {
    install.packages("devtools")
}
devtools::install_github("Vitek-Lab/MSstatsWeightedSummary")
```

