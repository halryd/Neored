000\_make\_directories.R
================
xrydbh
2022-02-04

``` r
library(tidyverse)
library(readxl)

# Three version of sample names are used: 
# Unit_no like: 19900088
# Plf_unique_no like: MS_01 or PEA_01
# SC_Donor High_1
# A conversion key is made with this script

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).



# Create output directory structure, if not allready existing
dir.create("RData")
```

    ## Warning in dir.create("RData"): 'RData' already exists

``` r
dir.create("RData/reactome_results")
```

    ## Warning in dir.create("RData/reactome_results"): 'RData/reactome_results'
    ## already exists

``` r
dir.create("out_r/Reactome", recursive = TRUE)
```

    ## Warning in dir.create("out_r/Reactome", recursive = TRUE): 'out_r/Reactome'
    ## already exists

``` r
dir.create("out_r/GO_BP", recursive = TRUE)
```

    ## Warning in dir.create("out_r/GO_BP", recursive = TRUE): 'out_r/GO_BP' already
    ## exists

``` r
dir.create("out_r/WGCNA/ms", recursive = TRUE)
```

    ## Warning in dir.create("out_r/WGCNA/ms", recursive = TRUE): 'out_r/WGCNA/ms'
    ## already exists

``` r
dir.create("out_r/WGCNA/pea", recursive = TRUE)
```

    ## Warning in dir.create("out_r/WGCNA/pea", recursive = TRUE): 'out_r/WGCNA/pea'
    ## already exists
