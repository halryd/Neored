library(tidyverse)

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).


##################
## Read in data
##################

# # read keratins
# Made a little change
# keratins.0.0 <- readxl::read_excel("./data/001_ms_filt_out_keratin_associated.xlsx")

# read protein centric
ms.deqms.0.0 <- readxl::read_excel("./out_r/013_ms_filt_norm_test_go_anno.xlsx") 
pea.t.test.0.0 <- readxl::read_excel("./out_r/013_pea_filt_norm_test_go_anno.xlsx") 

# read_pathway_centric
ms_camera_reactome.0.0 <- readxl::read_excel("./out_r/Reactome/reactome_results.xlsx")

# read pw per prot
# No significant pathway for pea
pw_per_prot.0.0 <- readxl::read_excel("./out_r/014_pw_per_prot.xlsx")


# read correlation module centric
# No  module significantly correlatied to CD34 conncentration
pea_wgcna.0.0 <- readxl::read_excel("./out_r/WGCNA/pea_black_module_genes_with_go.xlsx")

##################
## Format
##################
sign_ms_camera_reactome <- ms_camera_reactome.0.0 %>% 
  filter(Entities.FDR<0.05)

sign_pathways.0.0 <- sign_ms_camera_reactome %>% 
  pluck("Pathway.name")

# x <- "a1~!@#$%^&*(){}_+:\"<>?,./;'[]-=" #or whatever
sign_pathways.0.1 <- str_replace_all(sign_pathways.0.0, "[[:punct:]]", "_")

ms_pea_info_per_pw <- list()

for (pw in sign_pathways.0.0 ) {
  
  # get prots of a pathway
  pw_prots <- sign_ms_camera_reactome %>% 
    select(Pathway.name,Submitted.entities.found) %>% 
    filter(Pathway.name==pw) %>% 
    pluck(2) %>%
    str_split(";")%>% 
    pluck(1) 
 
  # Get ms data for pw prots
  Plasma_lipoprotein_remodeling_ms_test <- ms.deqms.0.0 %>% 
    filter(UniprotID%in%pw_prots)%>% 
    select(UniprotID,Gene.symbol,Gene.synonyms,logFC,sca.P.Value,Possible.contamination,GO_BP_terms)
  
  # Get pea data for pw prots
  Plasma_lipoprotein_remodeling_pea_test <- pea.t.test.0.0 %>% 
    filter(UniprotID%in%pw_prots)%>% 
    select(UniprotID,OlinkID,estimate,p.value,inWGCNAmodule_yesIs1) 
  
  # join ms and pea info
  ms_pea_info_per_pw[[pw]] <- Plasma_lipoprotein_remodeling_ms_test %>% 
    full_join(Plasma_lipoprotein_remodeling_pea_test, by=c("UniprotID"="UniprotID")) %>% 
    arrange(sca.P.Value) %>% 
    relocate(GO_BP_terms,.after = inWGCNAmodule_yesIs1)
  
} 


ms_pea_info_per_pw

names(ms_pea_info_per_pw) <- sign_pathways.0.1
##################
## Write data
##################

# write formatted output to file
ms_pea_info_per_pw %>% 
  writexl::write_xlsx(path = "./out_r/024_ms_pea_info_per_pw.xlsx") 


