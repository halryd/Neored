002b\_pea\_format\_normalize\_naFilter.R
================
xrydbh
2022-02-04

``` r
library(tidyverse)
library(rlist)
library(writexl)

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).



# read pea file
pea.raw.0.0 <- readxl::read_excel("./data/001_pea_raw.xlsx") 

missing_values_allowed_for_prot <- 0 #3 for with_at_least_5_vals (_in_high_or_low) or 
filter_missing_val_both_groups_or_either_or <- "either_or" #"either_or" # or "both"


# make longer

NPX.0.3 <- pea.raw.0.0 %>% 
  pivot_longer(cols=contains("_"),names_to="SC_Donor",values_to = "NPX") %>% 
  mutate(NPX=as.numeric(NPX)) %>% 
  mutate(LOD=as.numeric(LOD)) 

######################################
# Count numb of NAs in High and Low
# group by proteins
#################################3####
# Separate column SC_Donor into "High low Group" and "Unit 1-8"
NPX.0.3.5 <- NPX.0.3 %>%
  separate(SC_Donor,
           into = c("Group","Unit"),
           sep = "_")

# Get number of NPX below LOD (missing values) per gene and group
numb_NPX_vals_below_LOD.0.0 <- NPX.0.3.5 %>% 
  group_by(Gene.symbol,Group) %>% 
  summarise(NPX_below_LOD = sum(NPX<LOD)) %>% 
  pivot_wider(names_from = Group,values_from =NPX_below_LOD) %>% 
  dplyr::rename(Numb_NA_High=High,Numb_NA_Low=Low)
```

    ## `summarise()` has grouped output by 'Gene.symbol'. You can override using the
    ## `.groups` argument.

``` r
NPX.0.3.1 <- full_join(NPX.0.3,numb_NPX_vals_below_LOD.0.0 ,by="Gene.symbol")



# Replace UniprotID NA for Gene.symbol NT-proBNP and OlinkID  OID00131 and OID01214 by  P16860
# Also calculate average between OID00131 and OID01214

#what OlinkID have Corresponding UniprotID as NA
#NPX.0.0 %>% filter(is.na(UniProt )) %>% pluck("OlinkID") %>% base::unique()


# Get average for the two NT-proBNP Olink variants
NT_proBNP <- NPX.0.3.1  %>% 
  filter(Gene.symbol=="NT-proBNP") %>% 
  group_by(SC_Donor) %>% 
  mutate(NPX=(mean(NPX))) %>% # Calculate mean
  ungroup() %>% 
  filter(OlinkID=="OID00131") # keep only on of the two laternative OlinkIDs

# Replace duplicates with average
NPX.0.4 <- NPX.0.3.1  %>% 
  filter(Gene.symbol!="NT-proBNP") %>%
  bind_rows(NT_proBNP)

# Replace UniprotID NA for NT-proBNP with "P16860"
NPX.0.4$UniprotID[is.na(NPX.0.4$UniprotID)] <- "P16860"

# Remove duplicate UniprotIDs "Q8NEV9,Q14213" (IL27a/IL27b) with only one ID
NPX.0.4$UniprotID[NPX.0.4$UniprotID=="Q8NEV9,Q14213"] <- "Q8NEV9"
# Remove duplicate UniprotIDs "Q11128,P21217" (FUT3/FUT5) with only one ID
NPX.0.4$UniprotID[NPX.0.4$UniprotID=="Q11128,P21217"] <- "Q11128"


#NPX.0.4 %>% pivot_wider(names_from=SC_Donor,values_from = NPX)
########################################
# Normalize
########################################

# Group by sample and calculate sum NPX for all proteins
NPX.0.5 <- NPX.0.4 %>% 
  group_by(SC_Donor) %>% # Maybe here we should also group by panel?
  mutate(Samp_sum_unlogged_npx=sum(2^NPX,na.rm = T)) %>% # The sum of unlogged NPX within each sample
  mutate(Samp_sum_npx=sum(NPX,na.rm = T)) # The sum of unlogged NPX within each sample


# Calculate unlogged NPX relative to unlogged sample sum
NPX.0.6 <- NPX.0.5 %>% 
  mutate(Rel_unlog_npx = 2^NPX/Samp_sum_unlogged_npx)  # Relative unlogged NPX by dividing unlogged NPX with sum of unlogged NPX within each sample



#dplyr::select(-Samp_sum_unlogged_npx)
# Calculating sum of NPX and sum of npx relative to sample sum
NPX.0.6.5 <- NPX.0.6 %>% 
  ungroup() %>% # maybe uncoment
  mutate(Sum_unlogged_npx=sum(2^NPX,na.rm = T)) %>% # The sum of unlogged NPX within each sample
  mutate(Sum_rel_unlogged_npx=sum(Rel_unlog_npx,na.rm = T)) %>% # The sum of relative unlogged NPX within each sample
  mutate(fraction_sum_pre_post_norm = Sum_unlogged_npx/Sum_rel_unlogged_npx) # Get the fraction (a constant for all proteins) of sample sum for nps/rel_npx

# Take log2

NPX.0.7 <- NPX.0.6.5 %>% 
  group_by(SC_Donor) %>%  # maybe uncoment
  mutate(Rel_npx = log2(Rel_unlog_npx*fraction_sum_pre_post_norm))  # Multiply all value with the fraction (a constant) Relative NPX by taking log2



# A note on why NPX/Samp_sum_npx is correct:
# we want to relate the NPX value to the total. We wnt to have it as a fraction of the total
# We dont want the fraction of the Unlogged NPX and the Unlogged sum of NPX
# If we want to go through these fraction have to subtract instead of divide
# But in the end we want NPX to be treated the same way as MS data and by doing fraction at 
# the unlogged scale we endup wit similar distribution for shared
  # Calculate standardized log2 relative NPX

NPX.0.7.5 <- NPX.0.7 %>% 
  ungroup() %>%       # We dont want to standardize within each sample
  group_by(OlinkID)  # We want to standardize within each gene

# Calculate median centered
NPX.0.8 <- NPX.0.7.5 %>% 
  # ungroup() %>% 
  # group_by(Unit_no) %>%  
  mutate(Prot_median_npx=median(Rel_npx,na.rm = T)) %>% 
  mutate(Prot_median_cent_rel_npx=Rel_npx-Prot_median_npx)


########################################
# Filter based on NA
########################################
## Filter based on NA according to missing_values_allowed_for_prot parameter and weather to filter in both groups or either or group
if (filter_missing_val_both_groups_or_either_or=="either_or"){
  
  NPX.0.9 <- NPX.0.8  %>% 
    filter(Numb_NA_High<=missing_values_allowed_for_prot | Numb_NA_Low<=missing_values_allowed_for_prot )
  print(paste("Abundance data filtered so that there are at least ",8-missing_values_allowed_for_prot," samples with values in High or in Low",sep="" ))
  
} else if (filter_missing_val_both_groups_or_either_or=="both") {
  NPX.0.9  <- NPX.0.8  %>% 
    filter(Numb_NA_High<=missing_values_allowed_for_prot & Numb_NA_Low<=missing_values_allowed_for_prot )
  print(paste("Abundance data filtered so that there are at least ",8-missing_values_allowed_for_prot," samples with values in High and in Low",sep="" ))
  
}else{
  print("No valid filtering alternative")
}
```

    ## [1] "Abundance data filtered so that there are at least 8 samples with values in High or in Low"

``` r
# Make wider for table output for Excel
###################################################
# # Read in prepared column header dfs
# headers = load(file = paste("./RData/001_alt_samp_names.RData", sep=""))
# headers

PEA_column_head.0.1 <- readxl::read_excel("./data/PEA_column_head.xlsx") 


npx_transformations.0.0 <- NPX.0.9 %>% 
  ungroup() %>% 
  dplyr::select(OlinkID,
                UniprotID,
                Gene.symbol,
                SC_Donor,
                Prot_median_cent_rel_npx, # Here the unique column is selected
                Numb_NA_High,
                Numb_NA_Low,
                Panel
  ) %>% 
  distinct() %>% 
  pivot_wider(names_from = SC_Donor,values_from = "Prot_median_cent_rel_npx") %>% # make wider
  mutate_at(vars(-group_cols()), as.character) %>%   #make all columns as characters
  ungroup() %>% 
  dplyr::select(c("OlinkID","UniprotID","Gene.symbol",PEA_column_head.0.1$SC_Donor,"Numb_NA_High","Numb_NA_Low",Panel))



###########

prot_count <- dim(npx_transformations.0.0)[1]

npx_transformations_long <- NPX.0.9 %>% 
  ungroup() %>% 
  arrange(OlinkID,UniprotID)


# Save normalized wide
fn = paste("./out_r/002_norm_pea_",missing_values_allowed_for_prot,"_missing_vals_in_",filter_missing_val_both_groups_or_either_or,"_",prot_count ,".xlsx",sep="")
npx_transformations.0.0 %>%
  writexl::write_xlsx(path = fn)

# save list of df as RData
save(npx_transformations_long,npx_transformations.0.0, file = paste("./RData/002b_norm_pea.RData", sep=""))



##############################
### collect prot ID counts for NPX.0.3
##############################

after_repl_duplicates <- list()

after_repl_duplicates[["Gene.symbol"]] <- NPX.0.4$Gene.symbol %>% 
  unique() %>% 
  length()
after_repl_duplicates[["Uniprot"]] <- NPX.0.4$UniprotID %>% 
  unique() %>% 
  length()
after_repl_duplicates[["Olink"]] <- NPX.0.4$OlinkID  %>% 
  unique() %>% 
  length()

after_norm <- list()

after_norm[["Gene.symbol"]] <- NPX.0.8$Gene.symbol %>% 
  unique() %>% 
  length()
after_norm[["Uniprot"]] <- NPX.0.8$UniprotID %>% 
  unique() %>% 
  length()
after_norm[["Olink"]] <- NPX.0.8$OlinkID  %>% 
  unique() %>% 
  length()

after_NA_filt <- list()

after_NA_filt[["Gene.symbol"]] <- NPX.0.9 $Gene.symbol %>% 
  unique() %>% 
  length()
after_NA_filt[["Uniprot"]] <- NPX.0.9 $UniprotID %>% 
  unique() %>% 
  length()
after_NA_filt[["Olink"]] <- NPX.0.9 $OlinkID  %>% 
  unique() %>% 
  length()



# load prot ID counts
prot_id_counts <- load("./RData/000_pea_numb_unique_prot_IDs.RData")


pea_numb_unique_prot_IDs <- pea_numb_unique_prot_IDs %>% 
  list.append(after_repl_duplicates=after_repl_duplicates) %>% 
  list.append(after_norm=after_norm) %>% 
  list.append(after_NA_filt=after_NA_filt)





##############################
#Turn ID count list into tibble
##############################
# Turn list into matrix
pea_numb_unique_prot_IDs_mat <- pea_numb_unique_prot_IDs %>% 
  unlist() %>% # recursive = FALSE, use.names = FALSE
  matrix(nrow=length(pea_numb_unique_prot_IDs), byrow=TRUE) 

# Assign col names
colnames(pea_numb_unique_prot_IDs_mat) <- names(pea_numb_unique_prot_IDs[[1]])


# Turn matrix into tibble add col with step info
pea_prot_count_table <- pea_numb_unique_prot_IDs_mat %>% 
  as_tibble() %>% 
  mutate(Step=names(pea_numb_unique_prot_IDs)) %>% 
  relocate(Step,.before = Gene.symbol)

pea_prot_count_table
```

    ## # A tibble: 5 × 4
    ##   Step                   Gene.symbol Uniprot Olink
    ##   <chr>                        <int>   <int> <int>
    ## 1 input                          459     459   460
    ## 2 after_proj_spec_format         459     459   460
    ## 3 after_repl_duplicates          459     458   459
    ## 4 after_norm                     459     458   459
    ## 5 after_NA_filt                  417     417   417

``` r
write_xlsx(pea_prot_count_table,"./out_r/002_pea_prot_count_table.xlsx")
```
