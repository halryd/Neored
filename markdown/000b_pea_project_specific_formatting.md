000b\_pea\_project\_specific\_formatting.R
================
xrydbh
2022-02-04

``` r
library(OlinkAnalyze)
```

    ## Loading required package: car

    ## Loading required package: carData

    ## 
    ## Attaching package: 'car'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     recode

    ## The following object is masked from 'package:purrr':
    ## 
    ##     some

    ## Loading required package: emmeans

    ## Loading required package: ggfortify

    ## Loading required package: ggrepel

    ## Loading required package: lmerTest

    ## Loading required package: lme4

    ## Loading required package: Matrix

    ## 
    ## Attaching package: 'Matrix'

    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack

    ## 
    ## Attaching package: 'lmerTest'

    ## The following object is masked from 'package:lme4':
    ## 
    ##     lmer

    ## The following object is masked from 'package:stats':
    ## 
    ##     step

    ## Warning: replacing previous import 'car::recode' by 'dplyr::recode' when loading
    ## 'OlinkAnalyze'

``` r
library(tidyverse)
library(writexl)
library(broom)
library(rlist)



# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).

input_path <- "./data/pea.xlsx"
missing_values_allowed_for_prot <- 0 #3 for with_at_least_5_vals (_in_high_or_low) or 
filter_missing_val_both_groups_or_either_or <- "either_or" #"either_or" # or "both"

# Read in prepared column header dfs
PEA_column_head.0.1 <- readxl::read_excel("./data/PEA_column_head.xlsx") 


# Reading Olink NPX file exported from NPX Manager with read_NPX
NPX.0.0 <- read_NPX(filename = input_path) %>% 
  tibble::as_tibble() 
```

    ## Warning: `tbl_df()` was deprecated in dplyr 1.0.0.
    ## Please use `tibble::as_tibble()` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.

``` r
names(NPX.0.0)
```

    ##  [1] "SampleID"      "Index"         "OlinkID"       "UniProt"      
    ##  [5] "Assay"         "MissingFreq"   "Panel"         "Panel_Version"
    ##  [9] "PlateID"       "QC_Warning"    "LOD"           "NPX"          
    ## [13] "Normalization"

``` r
head(NPX.0.0)
```

    ## # A tibble: 6 × 13
    ##   SampleID   Index OlinkID UniProt Assay MissingFreq Panel Panel_Version PlateID
    ##   <chr>      <int> <chr>   <chr>   <chr> <chr>       <chr> <chr>         <chr>  
    ## 1 High_1         1 OID012… P42785  PRCP  0.27778000… Olin… v.3603        202010…
    ## 2 High_5         2 OID012… P42785  PRCP  0.27778000… Olin… v.3603        202010…
    ## 3 CONTROL_S…     3 OID012… P42785  PRCP  0.27778000… Olin… v.3603        202010…
    ## 4 Low_1          4 OID012… P42785  PRCP  0.27778000… Olin… v.3603        202010…
    ## 5 Low_5          5 OID012… P42785  PRCP  0.27778000… Olin… v.3603        202010…
    ## 6 CONTROL_S…     6 OID012… P42785  PRCP  0.27778000… Olin… v.3603        202010…
    ## # … with 4 more variables: QC_Warning <chr>, LOD <dbl>, NPX <dbl>,
    ## #   Normalization <chr>

``` r
# Missing UniprotId (for protein "NT-proBNP") is read in as the character "NA" and not as missing value NA 
NPX.0.0 %>% filter(UniProt=="NA")
```

    ## # A tibble: 36 × 13
    ##    SampleID  Index OlinkID UniProt Assay MissingFreq Panel Panel_Version PlateID
    ##    <chr>     <int> <chr>   <chr>   <chr> <chr>       <chr> <chr>         <chr>  
    ##  1 High_1        1 OID001… NA      NT-p… 0           Olin… v.6113        202010…
    ##  2 High_5        2 OID001… NA      NT-p… 0           Olin… v.6113        202010…
    ##  3 CONTROL_…     3 OID001… NA      NT-p… 0           Olin… v.6113        202010…
    ##  4 Low_1         4 OID001… NA      NT-p… 0           Olin… v.6113        202010…
    ##  5 Low_5         5 OID001… NA      NT-p… 0           Olin… v.6113        202010…
    ##  6 CONTROL_…     6 OID001… NA      NT-p… 0           Olin… v.6113        202010…
    ##  7 High_2        7 OID001… NA      NT-p… 0           Olin… v.6113        202010…
    ##  8 High_6        8 OID001… NA      NT-p… 0           Olin… v.6113        202010…
    ##  9 Low_2         9 OID001… NA      NT-p… 0           Olin… v.6113        202010…
    ## 10 Low_6        10 OID001… NA      NT-p… 0           Olin… v.6113        202010…
    ## # … with 26 more rows, and 4 more variables: QC_Warning <chr>, LOD <dbl>,
    ## #   NPX <dbl>, Normalization <chr>

``` r
# Replace any "NA" in df with NA to correct
NPX.0.0[NPX.0.0=="NA"] = NA

######################################
# Make wider and Check if any OlinkID is duplicated
######################################
# Make wider by spreading SC_donor (Unit/SampleID)
NPX.wide.0.0 <- NPX.0.0 %>% 
  dplyr::select(-c(Index,QC_Warning)) %>% # Filtering out QC warnings OBS!! Sample 19900199 Panel METABOLISM has QC warning
  # Olink says:"data points should be treated with caution"
  pivot_wider(names_from = SampleID,values_from=NPX)

# Check if any OlinkID is duplicated - There is none
duplicated.prots <- NPX.wide.0.0 $OlinkID[NPX.wide.0.0$OlinkID %>% duplicated()] 
duplicated.prots
```

    ## character(0)

``` r
###################################### end make wider
# Remove control samples from myNPX data
NPX.0.1 <- NPX.0.0 %>% 
  filter(!str_detect(SampleID, "^CONTROL")) # Filter out any row in SampleID startng with "CONTROL"

######################################
# Reorder, add headers and missing gene/protein names
######################################
# Rename and reorder columns
NPX.0.2 <- NPX.0.1 %>% 
  dplyr::select(OlinkID,
                UniprotID ="UniProt",
                Gene.symbol ="Assay",
                Unit_no = SampleID,
                NPX,
                LOD,
                Panel)

# Add Plf_unique_no (PEA_01) and SC_Donor(High_1) in addition to Unit_no(19900088) and reorder
# SC_Donor(High_1) is used to divide data set into High and Low. Otherwise this step would not be needed
NPX.0.3 <- full_join(NPX.0.2,
                     PEA_column_head.0.1,
                     by="Unit_no") %>%
  dplyr::select(OlinkID,
                UniprotID,
                Gene.symbol,
                SC_Donor,
                Unit_no,
                Plf_unique_no,
                NPX,
                LOD,
                Panel) 


# filter NPX table to only show rows with missing UniprotID
# NT-proBNP, dosn have any  UniprotID, see manuscript
NPX.missing.unip <- NPX.0.3[NPX.0.3$UniprotID %>% 
                              is.na(),]

# get vector of Olink IDs with missing Uniprot
# NT-proBNP, finns med i de två panelerna "CARDIOVASCULAR III" och "METABOLISM" med separata IDs
OlIDsWithNoUp.0.0 <- NPX.0.3$OlinkID[NPX.0.3$UniprotID %>% 
                                       is.na()] %>% 
  unique()
# Get vector of Uniprot Ids with no matching Gene Symbol- there is none
UpIDsWithNoGS.0.0 <- NPX.0.3$UniprotID[NPX.0.3$Gene.symbol %>% 
                                         is.na()] %>% 
  unique()

# Get distinct combinations of genenames and platforms
distintGeneNameComb.0.0 <- NPX.0.3 %>% 
  dplyr::select(OlinkID,UniprotID,Gene.symbol,Panel) %>% 
  distinct()
# Get duplicated Uniprot IDs
distintGeneNameComb.0.0$UniprotID[distintGeneNameComb.0.0$UniprotID %>% 
                                    duplicated()]
```

    ## [1] NA

``` r
# Get duplicated gene symbols
distintGeneNameComb.0.0$Gene.symbol[distintGeneNameComb.0.0$Gene.symbol %>% 
                                      duplicated()]
```

    ## [1] "NT-proBNP"

``` r
# Get platforms where Uniprot ID is appears for the second time
distintGeneNameComb.0.0$Panel[distintGeneNameComb.0.0$Gene.symbol %>% 
                                duplicated()]
```

    ## [1] "Olink METABOLISM"

``` r
# The only case where two different Olink ID are used for the same gene name between platforms seems to be for 
# NT-proBNP
# This also seems to be the only protein that occurs on two platforms
######################################

############################## Get raw data, not NA filt

pea_raw.0.0 <- NPX.0.3 %>% 
  #  vars_to_sel <-()
  ungroup() %>% 
  dplyr::select(OlinkID,
                UniprotID,
                Gene.symbol,
                SC_Donor,
                NPX, # Here the unique column is selected
                LOD,
                Panel
  ) %>% 
  distinct() %>% 
  pivot_wider(names_from = SC_Donor,values_from = "NPX") %>% # make wider
  mutate_at(vars(-group_cols()), as.character) %>%   #make all columns as characters
  ungroup() %>% 
  select(c("OlinkID","UniprotID","Gene.symbol",PEA_column_head.0.1$SC_Donor,"LOD","Panel"))

# Save raw
fn = paste("./data/001_pea_raw.xlsx",sep="")
pea_raw.0.0  %>%
  writexl::write_xlsx(path = fn)


# initiate list
pea_numb_unique_prot_IDs <- list()
##############################
# Collect number of IDS info into list
##############################
### collect prot ID count for NPX.0.0
pea_numb_unique_prot_IDs[["input"]][["Gene.symbol"]] <- NPX.0.0$Assay %>% 
  unique() %>% 
  length()
pea_numb_unique_prot_IDs[["input"]][["Uniprot"]] <- NPX.0.0$UniProt %>% 
  unique() %>% 
  length()
pea_numb_unique_prot_IDs[["input"]][["Olink"]] <- NPX.0.0$OlinkID %>% 
  unique() %>% 
  length()


### collect prot ID counts for NPX.0.3
after_proj_spec_format <- list()

after_proj_spec_format[["Gene.symbol"]] <- NPX.0.3$Gene.symbol %>% 
  unique() %>% 
  length()
after_proj_spec_format[["Uniprot"]] <- NPX.0.3$UniprotID %>% 
  unique() %>% 
  length()
after_proj_spec_format[["Olink"]] <- NPX.0.3$OlinkID  %>% 
  unique() %>% 
  length()

pea_numb_unique_prot_IDs <- pea_numb_unique_prot_IDs %>% list.append(after_proj_spec_format=after_proj_spec_format)


save(pea_numb_unique_prot_IDs,file = "./RData/000_pea_numb_unique_prot_IDs.RData")
```
