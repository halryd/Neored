library(tidyverse)
library(readxl)
library(AnnotationDbi)
library("Homo.sapiens")
library(GO.db)

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).


# remove old file



# read ms file
loaded_objs <- load("./RData/000_ms_numb_unique_prot_IDs.RData")

#### Remove Keratinization related
ms_list <- ms.raw.0.0$UniprotID

length(ms_list)

#Annotated proteisn with AnnotationDbi
ms_terms.long <- AnnotationDbi::select(Homo.sapiens,
                                       keys=ms_list,
                                       columns = c("TERM","ALIAS"),     #,"ONTOLOGY","DEFINITION"
                                       keytype = "UNIPROT") %>% 
  as_tibble()

ms_terms.long$UNIPROT %>% unique() %>% length()

# Rename variables of output and filter out only entries with BP annotation
# OBS if proteins has no BO annotation it will be lost
ms_terms.bp_only.wide <- ms_terms.long %>% 
  dplyr::rename(Gene.symbol=ALIAS) %>% 
  dplyr::rename(UniprotID=UNIPROT) %>% 
  filter(ONTOLOGY=="BP") %>% 
  dplyr::select(UniprotID,TERM) %>% 
  distinct() %>% 
  dplyr::group_by(UniprotID) %>% 
  dplyr::summarise(GO_term = paste(TERM, collapse = ",")) %>% 
  dplyr::rename(GO_BP_terms=GO_term)


# Filter in "K/eratin" associated, but not inflamm etc assoc
keratin_associated.0.0 <- ms_terms.bp_only.wide %>% 
  filter(grepl("Keratin|keratin", GO_BP_terms)) %>%  # remove any proteins with K/keratin in GO terms
  filter(!(grepl("inflamm|immune|hemopoies|hematopoie|growth factor|cytokine", GO_BP_terms))) %>% # but keep if contains any of word parts inflamm|immune|hemopoies|hematopoie|growth factor|cytokine
  filter(UniprotID!="Q02487") # we want to leave this one out since it exists in the Olink targeted platform


# # Filter in  associated to relevant words
# keratin_related_to_keep <- keratin_associated.0.0 %>% 
#   filter(!(grepl("inflamm|immune|hemopoies|hematopoie|growth factor|cytokine", GO_BP_terms)))

# Store in  variable the Keratins related that will be kept
keratin_related_to_keep.0.0 <- ms_terms.bp_only.wide %>% 
  filter(grepl("Keratin|keratin", GO_BP_terms)) %>%  # remove any proteins with K/keratin in GO terms
  filter((grepl("inflamm|immune|hemopoies|hematopoie|growth factor|cytokine", GO_BP_terms))) 


# Join keratin associated to discard  with ms.raw to get gene symobls
keratin_associated.0.1 <- ms.raw.0.0 %>% 
  dplyr::select(UniprotID,Gene.symbol,Gene.synonyms) %>% 
  inner_join(keratin_associated.0.0,by=c("UniprotID"="UniprotID"))


# Join keratin associated to discard  with ms.raw to get gene symobls
keratin_related_to_keep.0.1 <- ms.raw.0.0 %>% 
  dplyr::select(UniprotID,Gene.symbol,Gene.synonyms) %>% 
  inner_join(keratin_related_to_keep.0.0,by=c("UniprotID"="UniprotID"))

# ms_unip_without_keratin_ass  <- ms_terms.bp_only.wide %>% 
#   filter(!grepl("eratin", GO_BP_terms)) %>% 
#   pluck("UniprotID") %>% 
#   unique()
# 
# length(ms_unip_without_keratin_ass %>% unique())
# 
# 
# keratin_associated.0.0$UniprotID

# This will list proteins that do not have keratins in their name but also 
# proteins with no BP annotation
# ms.raw.0.1.alt <- ms.raw.0.0 %>% 
#   dplyr::filter(UniprotID%in%ms_unip_without_keratin_ass)

# Use keratin_associated to be discarde to filter ms.raw
ms.raw.0.1 <- ms.raw.0.0 %>% 
  dplyr::filter(!(UniprotID%in%keratin_associated.0.0$UniprotID))



# keratins <- ms.raw.0.0 %>% 
#   dplyr::filter(UniprotID%in%keratin_associated)

############################
# Remove proteins listed as contamination in "Proteomics study of human cord blood reticulocyte-derived exosomes" Varela et.al (DOI:10.1038/s41598-018-32386-2)
# Described/motivated in "The human plasma proteome: history, character, and diagnostic prospects"; (DOI: 10.1074/mcp.r200007-mcp200)
# This was done to focus on Exosome proteins. Should not be done in this analyis.
############################

# removed_in_vareal_et_al.0.0 <- readxl::read_excel("./data/potential_contaminants_diaz_varela.xlsx",sheet="potential contaminants removed",skip=3) 
# removed_uniprot_varela <- removed_in_vareal_et_al.0.0$Accession
# 
# ms.raw.0.2 <- ms.raw.0.1 %>% 
#   dplyr::filter(!(UniprotID%in%removed_uniprot_varela))

############################
# Save raw
############################
fn = paste("./data/001_ms_raw.xlsx",sep="")
ms.raw.0.1 %>%
  writexl::write_xlsx(path = fn)

# Save keratin associatied discarded, with some uniprot IDs filled in
fn = paste("./data/001_ms_filt_out_keratin_associated.xlsx",sep="")
keratin_associated.0.1  %>%
  writexl::write_xlsx(path = fn)

# save keratin_related_to_keep
fn = paste("./data/001_ms_filt_in_keratin_associated.xlsx",sep="")
keratin_related_to_keep.0.1  %>%
  writexl::write_xlsx(path = fn)


