library(AnnotationDbi)
library("Homo.sapiens")
library(GO.db)
library(tidyverse)

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).



my_loaded_objs <- load("./RData/12_formatted_out.RData")
my_loaded_objs

# Read WGCNA res
wgcna.0.0 <- readxl::read_excel("./out_r/WGCNA/pea_black_module_genes_with_go.xlsx") %>% 
  dplyr::select(UniprotID) %>% 
  mutate(inWGCNAmodule_yesIs1=1)



# Use the columns functinon of the AnnotaionDB to list the columns of the Homosapiens
columns(Homo.sapiens)

ms_list <- ms_filt_norm_test$UniprotID
ms_list_short <- ms_filt_norm_test$UniprotID[1:5]
pea_list <- pea_filt_norm_test$UniprotID

ms_terms.long <- AnnotationDbi::select(Homo.sapiens,
                      keys=ms_list,
                      columns = c("TERM","ALIAS"),     #,"ONTOLOGY","DEFINITION"
                      keytype = "UNIPROT",
                      multiVals="first") %>% 
  as_tibble()

ms_terms.bp_only.wide <- ms_terms.long %>% 
  dplyr::rename(Gene.symbol=ALIAS) %>% 
  dplyr::rename(UniprotID=UNIPROT) %>% 
  
  filter(ONTOLOGY=="BP") %>% 
  dplyr::select(UniprotID,TERM) %>% 
  distinct() %>% 
  dplyr::group_by(UniprotID) %>% 
  dplyr::summarise(GO_term = paste(TERM, collapse = ",")) %>% 
  dplyr::rename(GO_BP_terms=GO_term)

ms_filt_norm_test.go_anno <- ms_filt_norm_test %>% 
  full_join(ms_terms.bp_only.wide,by=c("UniprotID"="UniprotID"))

# PEA

pea_terms.long <- AnnotationDbi::select(Homo.sapiens,
                                  keys=pea_list,
                                  columns = c("TERM","ALIAS"),     #,"ONTOLOGY","DEFINITION"
                                  keytype = "UNIPROT") %>% 
  as_tibble()

pea_terms.bp_only.wide <- pea_terms.long %>% 
  dplyr::rename(Gene.symbol=ALIAS) %>% 
  dplyr::rename(UniprotID=UNIPROT) %>% 
  filter(ONTOLOGY=="BP") %>% 
  dplyr::select(UniprotID,TERM) %>% 
  distinct() %>% 
  dplyr::group_by(UniprotID) %>% 
  dplyr::summarise(GO_term = paste(TERM, collapse = ",")) %>% 
  dplyr::rename(GO_BP_terms=GO_term)

pea_filt_norm_test.go_anno  <- pea_filt_norm_test %>% 
  full_join(pea_terms.bp_only.wide,by=c("UniprotID"="UniprotID"))%>% 
  full_join(wgcna.0.0,by=c("UniprotID"="UniprotID")) %>% 
  mutate(inWGCNAmodule_yesIs1=tidyr::replace_na(inWGCNAmodule_yesIs1, 0))


# write formatted output to file
ms_filt_norm_test.go_anno  %>%
  distinct() %>% 
  writexl::write_xlsx(path = "./out_r/013_ms_filt_norm_test_go_anno.xlsx")

pea_filt_norm_test.go_anno  %>%
  writexl::write_xlsx(path = "./out_r/013_pea_filt_norm_test_go_anno.xlsx")

# write only GO to file



ms_terms.long  %>%
  distinct() %>% 
  writexl::write_xlsx(path = "./out_r/013_ms_go_terms.xlsx")
pea_terms.long  %>%
  writexl::write_xlsx(path = "./out_r/013_pea_go_terms.xlsx")

