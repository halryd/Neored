# csv files has to be downloaded manually Reactome website before collecting into multiple sheet Excel file with this script 
# This script takes a little bit time since it collected Gene.symbols

library(writexl)
library(tidyverse)
library(AnnotationDbi)
library("Homo.sapiens")

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).


reactome_results <- list()
reactome_results_gs <- list()

# 
# ms <- read.csv("./out_r/Reactome/ms_Camera.csv",sep=",",quote="\"") %>% 
#   as_tibble() 
# 
# ms[!is.na(ms$Submitted.entities.found),]
# ms[is.na(ms$Mapped.entities),]
# 
# pea <- read.csv("./out_r/Reactome/pea_Camera.csv",sep=",",quote="\"") %>% 
#   as_tibble()
# 
# pea[!is.na(pea$Submitted.entities.found),]
# pea[is.na(pea$Mapped.entities),]
# 
# 
# names(pea) %>% cat(sep="\n")
# 
# ms  %>% 
#   filter(Pathway.identifier == "R-HSA-173623")
# klutt <- ms  %>% 
#   filter(Pathway.identifier == "R-HSA-2730905")
# 
# klutt <- pea  %>% 
#   filter(Pathway.identifier == "R-HSA-2730905")
# 
# 
# pea <- read.csv("./out_r/Reactome/pea_Camera.csv",sep=",",quote="\"") %>% 
#   as_tibble() %>% 
#   filter(Pathway.identifier == "R-HSA-173623")

# These files have to be downloaded manually form Reactome

# reactome_results[["ms_PADOG"]] <- read.csv("./out_r/Reactome/ms_PADOG.csv",sep=",",quote="\"") %>% 
#   as_tibble()
# reactome_results[["pea_PADOG"]] <- read.csv("./out_r/Reactome/pea_PADOG.csv",sep=",",quote="\"") %>% 
#   as_tibble()
reactome_results[["ms_Camera"]] <- read.csv("./out_r/Reactome/ms_Camera.csv",sep=",",quote="\"") %>% 
  as_tibble() %>% 
  filter(Pathway.identifier != "R-HSA-173623")
reactome_results[["pea_Camera"]] <- read.csv("./out_r/Reactome/pea_Camera.csv",sep=",",quote="\"") %>% 
  as_tibble() %>% 
  filter(!Pathway.identifier %in% c("R-HSA-173623","R-HSA-2730905"))
#reactome_results[["ms_ssGSEA"]] <- read.csv("./out_r/Reactome/ms_ssGSEA.csv",sep=",",quote="\"") %>% 
#  as_tibble()
#reactome_results[["pea_ssGSEA"]] <- read.csv("./out_r/Reactome/pea_ssGSEA.csv",sep=",",quote="\"") %>% 
#  as_tibble()

##############################
# Add gene symbol names
##############################

# Use the columns functinon of the AnnotaionDB to list the columns of the Homosapiens
columns(Homo.sapiens)

# check pathways where mapped entities have missing value
# hepp <- reactome_results[[1]] %>% 
#   filter(Pathway.identifier == "R-HSA-173623")

for (platformID in names(reactome_results)){
  
  # MS
  # reactome_results[[platformID]] %>% 
  #   names()
  
  # check names
  # reactome_results[[platformID]] %>% 
  #   names()
  
  # Get get column with uniprotID as text string
  #uniprotID_col <- reactome_results[[platformID]]$Submitted.entities.found %>% head()
  pwID_and_unipID <- reactome_results[[platformID]] %>%
    #head() %>% 
    dplyr::select(Pathway.identifier,Submitted.entities.found)
  # Get uniprotID string for first pathway ID
  #path1_unip_semicol_sep_string <- uniprotID_col[1]
  #uniprotID_col[1,]
  # Conv first uniprot string to vector
  #path1_unip_vec <- unlist(strsplit(path1_unip_semicol_sep_string, split=";"))
  list_of_unip_vec_per_path <- strsplit(pwID_and_unipID$Submitted.entities.found, split = ";")
  
  # Convert semicolon separated string words  to one word value per row (https://stackoverflow.com/questions/15347282/split-delimited-strings-in-a-column-and-insert-as-new-rows)
  pw_uni_long <- data.frame(Pathway.identifier = rep(pwID_and_unipID$Pathway.identifier, sapply(list_of_unip_vec_per_path, length)), 
             Submitted.entities.found = unlist(list_of_unip_vec_per_path)) %>% 
    as_tibble() %>% 
    group_by(Pathway.identifier)%>% 
    filter(Pathway.identifier != "R-HSA-2029481") %>%  # removed pw because causing trouble, move back in if possible
    filter(Pathway.identifier != "R-HSA-2029485") %>% # The uniport IDs of these pathwways genreates error 
    filter(Pathway.identifier != "R-HSA-5690714") %>% # "None of the keys entered are valid keys for 'UNIPROT'. Please use the keys method to see a listing of valid arguments." When using it with AnnotationDbi::select
    filter(Pathway.identifier != "R-HSA-9664323")
  
  # Translate each uniprot to Gene.symbol (Usually one to many synonyms)
  # Output is table format with one uniprot per row
  # uni_gs.ms <- AnnotationDbi::select(Homo.sapiens,
  #                                        keys=path1_unip_vec ,
  #                                        columns = c("UNIPROT","ALIAS"),     #,"ONTOLOGY","DEFINITION"
  #                                        keytype = "UNIPROT") %>% 
  #   as_tibble() %>% 
  #   dplyr::rename(Gene.symbol=ALIAS) %>% 
  #   dplyr::rename(UniprotID=UNIPROT) %>% 
  #   dplyr::group_by(UniprotID) %>% 
  #   dplyr::summarise(Gene.symbol = paste(Gene.symbol, collapse = ",")) 
  
  # Make function that takes a vector of uniprot and outputs df of uniprot and gene symbols
  # 
  uniprot_vec_to_uniprot_gs_df <- function(unip_vec_of_pathw) {
    uni_gs.long <- AnnotationDbi::select(Homo.sapiens,
                                          keys=unip_vec_of_pathw ,
                                          columns = c("UNIPROT","ALIAS"),     #,"ONTOLOGY","DEFINITION"
                                          keytype = "UNIPROT") %>% 
      as_tibble() %>% 
      dplyr::rename(Gene.symbol=ALIAS) %>% 
      dplyr::rename(UniprotID=UNIPROT) %>% 
      dplyr::group_by(UniprotID) %>% 
      dplyr::summarise(Gene.symbol = paste(Gene.symbol, collapse = ",")) 
    
    # Collect uniprots and Gene.symbol synonyms into strings
    uni_gs.wide <- uni_gs.long  %>% 
      dplyr::summarise(
        UniprotID = paste(UniprotID, collapse = ";"),
        Gene.symbol = paste(Gene.symbol, collapse = ";")
      )
    
    gene.symbols <- uni_gs.wide %>% 
      pluck("Gene.symbol")
    
    return(gene.symbols)
  }
  
  # test function on one unipoprt vector
  #s.ms <- uniprot_vec_to_uniprot_gs_df(path1_unip_vec)
  # Map uniprot_vec_to_uniprot_gs_df function to every uniprot set of pathways
  # Use summarise
  print(paste("summarising: ",platformID,sep=""))
  unip_gs <- pw_uni_long %>% 
    summarise(Gene.symbols=uniprot_vec_to_uniprot_gs_df(Submitted.entities.found))
  # Use function
  # pw_uni_long %>% group_map(uniprot_vec_to_uniprot_gs_df,Submitted.entities.found)
  # # Use formula
  # pw_uni_long %>% group_map(~ uniprot_vec_to_uniprot_gs_df(Submitted.entities.found))
  
  # Collect uniprots and Gene.symbol synonyms into strings
  # uni_gs %>% 
  #   dplyr::summarise(
  #     UniprotID = paste(UniprotID, collapse = ";"),
  #     Gene.symbol = paste(Gene.symbol, collapse = ";")
  #     )
  # 
  # Convert vector with Uniprotids to Gene.symbols for first pathway 
  # uni_gs
  # Joining in Gene symbols  and adding to list
  reactome_results_gs[[platformID]] <- reactome_results[[platformID]] %>% 
    full_join(unip_gs,by=c("Pathway.identifier"="Pathway.identifier"))

}
##############################
# End addd gene symbol names
##############################


write_xlsx(
  reactome_results_gs,
  path = "./out_r/Reactome/reactome_results.xlsx",
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)

# save toptables
reactome_results_head <- reactome_results %>% 
  map(head)

reactome_results_tail <- reactome_results %>% 
  map(tail)

# write_xlsx(
#   reactome_results_head,
#   path = "./out_r/Reactome/reactome_results_tt.xlsx",
#   col_names = TRUE,
#   format_headers = TRUE,
#   use_zip64 = FALSE
# )
