library(tidyverse)
library(readxl)
library(AnnotationDbi)
library("Homo.sapiens")
library(GO.db)

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).



# read pea file
# loaded_objs <- load("./RData/000_pea_numb_unique_prot_IDs.RData")

001_pea_raw

pea.0.0 <- readxl::read_excel("./data/001_pea_raw.xlsx")#,sheet=2


#### Remove Keratinization related
pea_list <- pea.0.0$UniprotID

length(pea_list)

#Annotated proteisn with AnnotationDbi
pea_terms.long <- AnnotationDbi::select(Homo.sapiens,
                                       keys=pea_list,
                                       columns = c("TERM","ALIAS"),     #,"ONTOLOGY","DEFINITION"
                                       keytype = "UNIPROT") %>% 
  as_tibble()

pea_terms.long$UNIPROT %>% unique() %>% length()

# Rename variables of output and filter out only entries with BP annotation
# OBS if proteins has no BO annotation it will be lost
pea_terms.bp_only.wide <- pea_terms.long %>% 
  dplyr::rename(Gene.symbol=ALIAS) %>% 
  dplyr::rename(UniprotID=UNIPROT) %>% 
  filter(ONTOLOGY=="BP") %>% 
  dplyr::select(UniprotID,TERM) %>% 
  distinct() %>% 
  dplyr::group_by(UniprotID) %>% 
  dplyr::summarise(GO_term = paste(TERM, collapse = ",")) %>% 
  dplyr::rename(GO_BP_terms=GO_term)


# Filter in "eratin" associated
keratin_associated.0.0 <- pea_terms.bp_only.wide %>% 
  filter(grepl("Keratin|keratin", GO_BP_terms)) %>% 
  filter(!(grepl("inflamm|immune|hemopoies|hematopoie|growth factor|cytokine", GO_BP_terms))) %>% # but keep if contains any of word parts inflamm|immune|hemopoies|hematopoie|growth factor|cytokine
  filter(GO_BP_terms!="Q02487") # we want to leave this one out since it exists in the Olink targeted platform

# Filter in  associated to relevant words
keratin_related_to_keep <- keratin_associated.0.0 %>% 
  filter(grepl("inflamm|immune|hemopoies|hematopoie|growth factor|cytokine", GO_BP_terms)) 

 
# Join  with ms.raw
keratin_associated.0.1 <- pea.0.0 %>% 
  dplyr::select(UniprotID,Gene.symbol,OlinkID) %>% 
  inner_join(keratin_associated.0.0,by=c("UniprotID"="UniprotID"))

# Save keratin associatied, with some uniprot IDs filled in
fn = paste("./data/001_pea_filt_out_keratin_associated.xlsx",sep="")
keratin_associated.0.1  %>%
  writexl::write_xlsx(path = fn)

