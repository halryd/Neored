library(readxl)
library(tidyverse)
library(rlist)
#library(ggpubr)

# testing

# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).


# specify file paths
input_path <- "./data/ms.xlsx"
missing_values_allowed_for_prot <- 0 #3 for with_at_least_5_vals (_in_high_or_low) or 
filter_missing_val_both_groups_or_either_or <- "either_or" #"either_or" # or "both"

# read ms file
mass_spec.0.1 <- readxl::read_excel(input_path ,
                                    sheet = "USE normalized abundances",
                                    col_names = TRUE,
                                    na = "NaN",
                                    skip = 6) 
# make syntactically valid names
names(mass_spec.0.1) <- 
  mass_spec.0.1 %>% 
  names() %>% 
  make.names()


# rename and reorder columns
mass_spec.0.2 <- mass_spec.0.1 %>% 
  dplyr::select(UniprotID="Accession",
                Gene.synonyms="Gene.Symbol",
                contains("SC.Donor"),
                Coverage="Coverage....",
                No.peptides="X..Peptides",
                PSMs="X..PSMs",
                Unique.Peptides="X..Unique.Peptides",
                MW.kDa="MW..kDa.") %>% 
  #separate_rows(Gene.symbol, sep = ";") %>% # put each gene symbol on its onw row
  mutate(Gene.synonyms=str_replace_all(Gene.synonyms, " ", "")) %>%  #remove spaces (replace spaces with nothing)
  mutate(Gene.symbol=word(Gene.synonyms, start = 1,sep = "\\;")) %>% #Create Gene.symbol from Gene.synonyms by extracting first word in sentence separated by semicolon
  relocate(Gene.symbol,.before = Gene.synonyms) # Put added Gene.symbol next to missing Gene.synonyms

# find ID with Gene.synonyms
symbols_with_synonyms <- mass_spec.0.2 %>% 
  pluck("Gene.synonyms") %>% 
  str_subset(";") %>% # Keep strings matching a ;
  word(start = 1,sep = "\\;")  #Create Gene.symbol from Gene.synonyms by extracting first word in sentence separated by semicolon

# See that Genes with multiple synonyms have been annotated in Gene.synonyms
mass_spec.0.2 %>% 
  select(UniprotID,
         Gene.symbol,
         Gene.synonyms) %>% 
  filter(Gene.symbol %in% symbols_with_synonyms)
##################################
# Check if Uniprot is unique
mass_spec.0.2$UniprotID[mass_spec.0.2$UniprotID %>% duplicated()] 
##################################

# shorten column names
names(mass_spec.0.2) <- str_replace(names(mass_spec.0.2), ".SC.Donor","_")
names(mass_spec.0.2)

# make longer by collecting SC_donor
mass_spec.0.3 <- mass_spec.0.2 %>% 
  pivot_longer(cols=contains("_"),
               names_to=c("SC_Donor"),
               values_to = "Norm_abundance") %>% 
  mutate(Norm_abundance=as.numeric(Norm_abundance)) # make abundances numeric

#####################################
# Fill in missing gene symbols
#####################################
# Get vector of Uniprot Ids with no matching Gene Symbol
UpIDsWithNoGS.0.0 <- mass_spec.0.3$UniprotID[mass_spec.0.3$Gene.symbol %>% is.na()] %>% # Which Uniport has missing values?
  unique()

######## Get uniprot annotation
####
# Define function uniprot_mapping
uniprot_mapping <- function(ids) {
  uri <- 'http://www.uniprot.org/uniprot/?query='
  idStr <- paste(ids, collapse="+or+")
  format <- '&format=tab'
  fullUri <- paste0(uri,idStr,format)
  dat <- read.delim(fullUri)
  dat
}

# Full anno table
Uniprot_anno.0.0 <- uniprot_mapping(UpIDsWithNoGS.0.0) %>% 
  as_tibble() # get Uniprot annotation corresponding to UPIDs

# For some reason, output is has too many Uniprot IDs??#Restrict rows
Uniprot_anno.0.0.1 <- Uniprot_anno.0.0 %>%  
  filter(Entry %in% UpIDsWithNoGS.0.0) #%>% 
#filter(Entry != "P0DOX7") # P0DOX7 Has no annotation

# Restricted anno table (to two cols)
Uniprot_anno.0.1 <- Uniprot_anno.0.0.1 %>% 
  dplyr::select(UniprotID=Entry,
                Gene.synonyms.add=Gene.names) %>% # restrict annotation to only UPIDs and GS
  #separate_rows(Gene.symbol, sep = " ") # put each GS on its own row
  mutate(Gene.synonyms.add=str_replace_all(Gene.synonyms.add,pattern = " ",replacement = ";")) %>% 
  mutate(Gene.symbol.add=word(Gene.synonyms.add,1,sep = ";")) %>% # Create Gene.symbol with only one selected Gene name from Gene.synonyms
  relocate(Gene.symbol.add,.before = Gene.synonyms.add) # Put added GS next to missing GS
## If gene name is missing (NA), replace with fake gene name "no.gene.symbol"
Uniprot_anno.0.1$Gene.symbol.add[Uniprot_anno.0.1$Gene.symbol.add==""] <- "no.gene.symbol"
Uniprot_anno.0.1$Gene.synonyms.add [Uniprot_anno.0.1$Gene.synonyms.add==""] <- "no.gene.symbol"

###############

# Merge dataframe of proteins with no GS with With Uniprot Anno
mass_spec.0.3.1 <- full_join(mass_spec.0.3,Uniprot_anno.0.1,by="UniprotID") %>% 
  relocate(Gene.symbol.add,.after = Gene.symbol) %>%  # Put added GS next to missing GS
  relocate(Gene.synonyms.add,.after = Gene.synonyms) # Put added GS next to missing GS

# When Gene.symbol is NA replace with Gene.symbol.add
mass_spec.0.3.2 <- mass_spec.0.3.1 %>% 
  mutate(Gene.symbol = case_when(is.na(Gene.symbol) ~ Gene.symbol.add,
                                 TRUE ~ Gene.symbol)) %>% 
  mutate(Gene.synonyms = case_when(is.na(Gene.synonyms) ~ Gene.symbol, # When there are no synonyms jus put Gene.symbol in synonyms
                                   TRUE ~ Gene.synonyms))


mass_spec.0.3.3 <- mass_spec.0.3.2 %>% 
  dplyr::select(-c(Gene.symbol.add,Gene.synonyms.add)) # remove Gene.symbol.y

######## Get raw data, no NA filt into wide

ms.raw.0.0 <- mass_spec.0.3.3 %>% 
  dplyr::select(UniprotID,
                Gene.symbol,
                Gene.synonyms,
                SC_Donor,
                Norm_abundance, # Here the unique column is selected
                No.peptides,
                PSMs,
                Unique.Peptides,
                MW.kDa,
                
  ) %>% 
  distinct() %>% 
  pivot_wider(names_from = SC_Donor,values_from = "Norm_abundance") %>% # make wider
  mutate_at(vars(-group_cols()), as.character) %>%  #make all columns as characters
  ungroup()

# Save raw, with some uniprot IDs filled in
# This is now saved after filtering out Keratins instead
 # fn = paste("./data/001_ms_raw.xlsx",sep="")
 # ms.raw.0.0 %>%
 #   writexl::write_xlsx(path = fn)

# initiate list
ms_numb_unique_prot_IDs <- list()
##############################
# Collect number of IDS info into list
##############################
### collect prot ID count for NPX.0.0
ms_numb_unique_prot_IDs[["input"]][["Gene.symbol"]] <- mass_spec.0.1$Gene.Symbol %>% 
  unique() %>% 
  length()
ms_numb_unique_prot_IDs[["input"]][["Uniprot"]] <- mass_spec.0.1$Accession %>% 
  unique() %>% 
  length()


### collect prot ID counts for NPX.0.3
after_proj_spec_format <- list()

after_proj_spec_format[["Gene.symbol"]] <- mass_spec.0.3.3$Gene.symbol %>% 
  unique() %>% 
  length()
after_proj_spec_format[["Uniprot"]] <- mass_spec.0.3.3$UniprotID %>% 
  unique() %>% 
  length()


ms_numb_unique_prot_IDs <- ms_numb_unique_prot_IDs %>% 
  list.append(after_proj_spec_format=after_proj_spec_format)

save(ms_numb_unique_prot_IDs,ms.raw.0.0,file = "./RData/000_ms_numb_unique_prot_IDs.RData")
