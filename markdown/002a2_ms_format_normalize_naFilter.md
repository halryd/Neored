002a2\_ms\_format\_normalize\_naFilter.R
================
xrydbh
2022-02-04

``` r
library(readxl)
library(writexl)
library(tidyverse)
library(rlist)


# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).



missing_values_allowed_for_prot <- 0 #3 for with_at_least_5_vals (_in_high_or_low) or 
filter_missing_val_both_groups_or_either_or <- "either_or" #"either_or" # or "both"


# read ms file
ms.raw.0.0 <- readxl::read_excel("./data/001_ms_raw.xlsx") 


# make longer

mass_spec.0.3.3 <- ms.raw.0.0 %>% 
  pivot_longer(cols=contains("_"),names_to="SC_Donor",values_to = "Norm_abundance") %>% 
  mutate(Norm_abundance=as.numeric(Norm_abundance))

############################################
# Count numb of NAs in High and Low
# group by proteins
############################################
# Separate column SC_Donor into "High low Group" and "Unit 1-8"
mass_spec.0.3.4 <- mass_spec.0.3.3 %>%
  separate(SC_Donor,
           into = c("Group","Unit"),
           sep = "_")

# Get missing values per gene and group
missing_values_per_gene.0.0 <- mass_spec.0.3.4 %>% 
  group_by(Gene.symbol,Group) %>%           # group by gene symbol and High/Low
  summarise(Missing_ab_val = sum(is.na(Norm_abundance))) %>%  # get number of NAs per gene and High/Low groups
  pivot_wider(names_from = Group,values_from =Missing_ab_val) %>% 
  dplyr::rename(Numb_NA_High=High,Numb_NA_Low=Low)
```

    ## `summarise()` has grouped output by 'Gene.symbol'. You can override using the
    ## `.groups` argument.

``` r
mass_spec.0.3.5 <- full_join(mass_spec.0.3.3,missing_values_per_gene.0.0,by="Gene.symbol")


########################################
# Normalize 
########################################
mass_spec.0.4 <- mass_spec.0.3.5 %>% 
  mutate(Log2_norm_ab = log2(Norm_abundance)) %>% 
  ungroup() %>%       # We dont want to standardize within each sample
  group_by(UniprotID) %>%  # We want to standardize within each gene
  mutate(Prot_median_l2_norm_ab=median(Log2_norm_ab,na.rm = T)) %>% 
  mutate(Prot_median_cent_l2_norm_ab=Log2_norm_ab-Prot_median_l2_norm_ab)

########################################
# Filter based on NA
########################################
## Filter based on NA according to missing_values_allowed_for_prot parameter and weather to filter in both groups or either or group
if (filter_missing_val_both_groups_or_either_or=="either_or"){
  
  mass_spec.0.5 <- mass_spec.0.4 %>% 
    filter(Numb_NA_High<=missing_values_allowed_for_prot | Numb_NA_Low<=missing_values_allowed_for_prot )
  print(paste("Abundance data filtered so that there are at least ",8-missing_values_allowed_for_prot," samples with values in High or in Low",sep="" ))
  
} else if (filter_missing_val_both_groups_or_either_or=="both") {
  mass_spec.0.5 <- mass_spec.0.4 %>% 
    filter(Numb_NA_High<=missing_values_allowed_for_prot & Numb_NA_Low<=missing_values_allowed_for_prot )
  print(paste("Abundance data filtered so that there are at least ",8-missing_values_allowed_for_prot," samples with values in High and in Low",sep="" ))
  
}else{
  print("No valid filtering alternative")
}
```

    ## [1] "Abundance data filtered so that there are at least 8 samples with values in High or in Low"

``` r
dim(mass_spec.0.5)
```

    ## [1] 11920    14

``` r
# Make wider for table output for Excel
###################################################
# Use mass_spec.0.8 (for mass_spec.0.9 and further gene.symbol has been removed)
# collect interesting data transformations/normalised variables (column/variables names) in vector


norm_ab_transformations.0.0 <- mass_spec.0.5 %>% # Use mass_spec.0.8 (for mass_spec.0.9 and further gene.symbol has been removed)
  dplyr::select(UniprotID,
                Gene.symbol,
                Gene.synonyms,
                SC_Donor,
                Prot_median_cent_l2_norm_ab, # Here the unique column is selected
                Numb_NA_High,
                Numb_NA_Low,
                No.peptides,
                PSMs,
                Unique.Peptides,
                MW.kDa,
  ) %>% 
  distinct() %>% 
  pivot_wider(names_from = SC_Donor,values_from = "Prot_median_cent_l2_norm_ab") %>% # make wider
  mutate_at(vars(-group_cols()), as.character) %>%  #make all columns as characters
  ungroup()




prot_count <- dim(norm_ab_transformations.0.0)[1]

cat(paste("The number of proteins after applying criteria of","allowed fraction of missing values is:",prot_count,sep="\n"))
```

    ## The number of proteins after applying criteria of
    ## allowed fraction of missing values is:
    ## 745

``` r
norm_ab_transformations_long <- mass_spec.0.5 %>% 
  ungroup()

# remove old file
file_to_remove <- paste("002_norm_ms_",missing_values_allowed_for_prot,"_missing_vals_in_",filter_missing_val_both_groups_or_either_or,sep="")
junk.0 <- dir(path="./out_r/",  pattern=file_to_remove) # ?dir
junk <- paste("./out_r/",junk.0,sep="")
file.remove(junk) # ?file.remove
```

    ## [1] TRUE

``` r
# Save normalized wide
fn = paste("./out_r/002_norm_ms_",missing_values_allowed_for_prot,"_missing_vals_in_",filter_missing_val_both_groups_or_either_or,"_",prot_count ,".xlsx",sep="")
norm_ab_transformations.0.0 %>%
  writexl::write_xlsx(path = fn)


###########write data files
save(norm_ab_transformations_long,norm_ab_transformations.0.0, file = paste("./RData/002a_norm_ms.RData", sep=""))


##############################
### collect prot ID counts for MS
##############################


after_norm <- list()

after_norm[["Gene.symbol"]] <- mass_spec.0.4$Gene.symbol %>% 
  unique() %>% 
  length()
after_norm[["Uniprot"]] <- mass_spec.0.4$UniprotID %>% 
  unique() %>% 
  length()

after_NA_filt <- list()

after_NA_filt[["Gene.symbol"]] <- mass_spec.0.5$Gene.symbol %>% 
  unique() %>% 
  length()
after_NA_filt[["Uniprot"]] <- mass_spec.0.5$UniprotID %>% 
  unique() %>% 
  length()


# load prot ID counts
prot_id_counts <- load("./RData/000_ms_numb_unique_prot_IDs.RData")


ms_numb_unique_prot_IDs <- ms_numb_unique_prot_IDs %>% 
  list.append(after_norm=after_norm) %>% 
  list.append(after_NA_filt=after_NA_filt)





##############################
#Turn ID count list into tibble
##############################
# Turn list into matrix
ms_numb_unique_prot_IDs_mat <- ms_numb_unique_prot_IDs %>% 
  unlist() %>% # recursive = FALSE, use.names = FALSE
  matrix(nrow=length(ms_numb_unique_prot_IDs), byrow=TRUE) 

# Assign col names
colnames(ms_numb_unique_prot_IDs_mat) <- names(ms_numb_unique_prot_IDs[[1]])


# Turn matrix into tibble add col with step info
ms_prot_count_table <- ms_numb_unique_prot_IDs_mat %>% 
  as_tibble() %>% 
  mutate(Step=names(ms_numb_unique_prot_IDs)) %>% 
  relocate(Step,.before = Gene.symbol)

ms_prot_count_table
```

    ## # A tibble: 4 × 3
    ##   Step                   Gene.symbol Uniprot
    ##   <chr>                        <int>   <int>
    ## 1 input                         1396    1419
    ## 2 after_proj_spec_format        1419    1419
    ## 3 after_norm                    1356    1356
    ## 4 after_NA_filt                  745     745

``` r
write_xlsx(ms_prot_count_table,"./out_r/002_ms_prot_count_table.xlsx")
```
