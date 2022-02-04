# reactomeGSA (multi-omics Gene Set Analysis)
# visualize with  ReactomeGSA

library(ReactomeGSA)
library(tidyverse)
library(DEqMS)
#library(ReactomeContentService4R)

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).



my_methods <- c("Camera") # ,"Camera""ssGSEA", "PADOG",

# Get abundances
loaded_stuff <- load("./RData/003_ab_list.RData")
loaded_stuff

# define function that creates ReactomeAnalysisRequest and submit it to Reactome

perform_reactomeGSA <- function(my_data.0.0,list_index_of_my_data){ # my_data.0.0,my_plf
  for (my_method in my_methods){
    #my_data.0.0 = .x
    my_data =  my_data.0.0 %>% # dat_abundances[[1]] %>% #my_data.0.0 %>% # 
      dplyr::select(-1) %>% 
      as.data.frame()
    row.names(my_data) = my_data.0.0 %>% #my_data %>%  #my_data.0.0 %>% 
      pluck(1) 
    #my_data
    dim(my_data)
    # specify groups
    
    cond = c(rep("High",ncol(my_data)/2),rep("Low",ncol(my_data)/2))
    
    # remove old request
    # rm(neored_ms_request)
    
    # Create a new request object using 'PADOG' for the gene set analysis
    neored_ms_request <- ReactomeAnalysisRequest(method = my_method)
    neored_ms_request
    
    # set the maximum number of allowed missing values to 100%, have already decided hwo to filter NAs
    neored_ms_request <- set_parameters(request = neored_ms_request, max_missing_values = 1)
    neored_ms_request
    
    
    neored_ms_sample_data <- as.data.frame(bind_cols(condition=as.factor(cond))) #patient.id=as.factor(names(my_data)),
    row.names(neored_ms_sample_data) <- names(my_data)
    neored_ms_sample_data %>% as_tibble()
    
    analysis_name <- paste(my_method,"_",list_index_of_my_data,"_neored",sep="")
    
    #ms_in_DEqMS_data
    neored_ms_request <- add_dataset(request = neored_ms_request, 
                                     expression_values = my_data,
                                     name = analysis_name, #"Proteomics"
                                     type = "proteomics_int",
                                     comparison_factor = "condition", 
                                     comparison_group_1 = "Low", 
                                     comparison_group_2 = "High",
                                     sample_data = neored_ms_sample_data,
                                     #additional_factors = c("patient.id")#"cell.type", 
    )
    
    
    # neored_ms_request 
    
    result <- perform_reactome_analysis(request = neored_ms_request,verbose = T)#, compress = F
    
    # save the result
    saveRDS(result, file = paste("./RData/reactome_results/ReactomeGSA_result_",my_method,"_",list_index_of_my_data,".rds",sep=""))#
    
  }
 
}

# iwalk is part of the purrr package, it applies the function perform_reactomeGSA on the two list items of ab_list[["dat_abundances"]] 

# # Remove ene.symbol from Olink
ab_list[["dat_abundances"]][[1]] <- ab_list[["dat_abundances"]][[1]] %>% 
  dplyr::select(-c(Gene.symbol))
# Remove Panel and gene.symbol from Olink
ab_list[["dat_abundances"]][[2]] <- ab_list[["dat_abundances"]][[2]] %>% 
  dplyr::select(-c(Panel,Gene.symbol))


names(ab_list[["dat_abundances"]][1:2]) <- c("ms","pea")
ab_list[["dat_abundances"]][1:2] %>% iwalk(perform_reactomeGSA)



ab_list[["dat_abundances"]][[1]]

