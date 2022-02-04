0\_install\_packages.R
================
xrydbh
2022-02-04

``` r
###########################
# Install form The Comprehensive R Archive Network (CRAN)
###########################

# First specify the packages of interest
baseR_packages <- c("tidyverse","writexl","broom","heatmap3","openxlsx","qgraph","rlist","sjmisc")

## Then load or install&load all
package.check <- lapply(
  baseR_packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE, repos="https://ftpmirror1.infania.net/mirror/CRAN/",source="mac.binary")
      library(x, character.only = TRUE)
    }
  }
)
```

    ## Loading required package: tidyverse

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
    ## ✓ tibble  3.1.6     ✓ dplyr   1.0.7
    ## ✓ tidyr   1.1.4     ✓ stringr 1.4.0
    ## ✓ readr   2.1.1     ✓ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    ## Loading required package: writexl

    ## Loading required package: broom

    ## Loading required package: heatmap3

    ## Loading required package: openxlsx

    ## Loading required package: qgraph

    ## Loading required package: rlist

    ## Loading required package: sjmisc

    ## 
    ## Attaching package: 'sjmisc'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     is_empty

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     replace_na

    ## The following object is masked from 'package:tibble':
    ## 
    ##     add_case

``` r
###########################
# install from Bioconductor
###########################

# First install BiocManager

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12",update = F)
```

    ## Bioconductor version 3.12 (BiocManager 1.30.16), R 4.0.5 (2021-03-31)

``` r
# Then specify the packages of interest
bioconductor_packages <- c("ReactomeGSA", "WGCNA", "EnhancedVolcano", "clusterProfiler", "pathview", "enrichplot","org.Hs.eg.db")

# Then install/update bioconductor_packages. This one has, as far as I know, to be commented out 
# when packages has been installed once to avoid spending time on re-installing
BiocManager::install(bioconductor_packages,update=F)
```

    ## Bioconductor version 3.12 (BiocManager 1.30.16), R 4.0.5 (2021-03-31)

    ## Warning: package(s) not installed when version(s) same as current; use `force = TRUE` to
    ##   re-install: 'ReactomeGSA' 'WGCNA' 'EnhancedVolcano' 'clusterProfiler'
    ##   'pathview' 'enrichplot' 'org.Hs.eg.db'

``` r
###########################
# install with Devtools
###########################

install.packages("devtools", repos="https://ftpmirror1.infania.net/mirror/CRAN/")
```

    ## 
    ## The downloaded binary packages are in
    ##  /var/folders/wc/0fht_gvd3js11d262_qdmzs00000gp/T//RtmpdnbGfS/downloaded_packages

``` r
devtools::install_github("stephenturner/annotables")
```

    ## Skipping install of 'annotables' from a github remote, the SHA1 (631423c3) has not changed since last install.
    ##   Use `force = TRUE` to force installation

``` r
#igraph??
#ggrepel??
baseR_packages
```

    ## [1] "tidyverse" "writexl"   "broom"     "heatmap3"  "openxlsx"  "qgraph"   
    ## [7] "rlist"     "sjmisc"

``` r
knitr::write_bib(c(baseR_packages,bioconductor_packages,"annotables"), "./data/neored_packages.bib")
```

    ## Warning in utils::citation(..., lib.loc = lib.loc): no date field in DESCRIPTION
    ## file of package 'annotables'

    ## Warning in utils::citation(..., lib.loc = lib.loc): no date field in DESCRIPTION
    ## file of package 'org.Hs.eg.db'
