005\_gsea\_camera\_ReactomeGSA.R
================
xrydbh
2022-02-04

``` r
# reactomeGSA (multi-omics Gene Set Analysis)
# visualize with  ReactomeGSA

library(ReactomeGSA)
```

    ## Registered S3 method overwritten by 'gplots':
    ##   method         from     
    ##   reorder.factor DescTools

    ## 
    ## Attaching package: 'ReactomeGSA'

    ## The following object is masked from 'package:qgraph':
    ## 
    ##     pathways

``` r
library(tidyverse)
library(DEqMS)
#library(ReactomeContentService4R)

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).



my_methods <- c("Camera") # ,"Camera""ssGSEA", "PADOG",

# Get abundances
loaded_stuff <- load("./RData/003_ab_list.RData")
loaded_stuff
```

    ## [1] "ab_list"

``` r
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
```

    ## Submitting request to Reactome API...

    ## Compressing request data...

    ## Reactome Analysis submitted succesfully

    ## Converting dataset Camera_ms_abundances_eq_med_norm_neored...

    ## Mapping identifiers...

    ## Performing gene set analysis using Camera

    ## Analysing dataset 'Camera_ms_abundances_eq_med_norm_neored' using Camera

    ## Retrieving result...

    ## Submitting request to Reactome API...

    ## Compressing request data...

    ## Reactome Analysis submitted succesfully

    ## Converting dataset Camera_pea_abundances_eq_med_norm_neored...

    ## Mapping identifiers...

    ## Performing gene set analysis using Camera

    ## Creating REACTOME visualization

    ## Retrieving result...

``` r
ab_list[["dat_abundances"]][[1]]
```

    ## # A tibble: 745 × 17
    ##    UniprotID High_1 High_2 High_3 High_4 High_5 High_6 High_7 High_8 Low_1 Low_2
    ##    <chr>     <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr> <chr>
    ##  1 P04114    0.046… -0.01… 0.217… 0.801… -0.15… 0.436… -0.03… -0.50… -0.5… 0.00…
    ##  2 P01024    0.478… -0.25… -0.01… -0.21… 0.171… -0.32… 0.009… 0.066… -0.2… -2.4…
    ##  3 P21333    -1.26… -0.46… 0.177… -0.04… -0.59… 1.119… 0.237… 0.162… -1.3… -1.9…
    ##  4 Q9Y490    -0.81… -0.00… 0.240… 0.500… -0.39… 2.077… 0.576… 0.733… -0.9… -0.7…
    ##  5 P02751    -0.31… -0.21… 0.507… -0.39… -0.14… -0.65… 0.222… 1.003… 0.44… -2.1…
    ##  6 P01031    0.410… -0.07… 0.024… -0.14… 0.149… -0.23… -0.06… 0.037… -0.3… -1.4…
    ##  7 P35579    -1.27… 0.268… 0.680… 0.532… -0.62… 0.698… -0.14… -0.27… -0.9… 0.72…
    ##  8 Q9Y6R7    -0.29… -0.60… -0.55… -0.17… 0.532… 0.198… 0.640… 0.355… -0.0… -0.9…
    ##  9 P18206    -0.93… 0.158… 0.139… 0.373… -0.73… 1.453… 0.167… 0.236… -0.9… -0.7…
    ## 10 P22105    -0.15… -0.20… -0.20… -0.31… -0.00… -0.25… 0.513… 0.489… 0.29… -2.3…
    ## # … with 735 more rows, and 6 more variables: Low_3 <chr>, Low_4 <chr>,
    ## #   Low_5 <chr>, Low_6 <chr>, Low_7 <chr>, Low_8 <chr>
