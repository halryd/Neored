012\_collect\_filt\_norm\_data\_and\_testresults.R
================
xrydbh
2022-02-04

``` r
library(tidyverse)

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).



####################################
both_dat <- load("./RData/003_ab_list.RData")
# both_dat
# names(ab_list)
# names(ab_list[["dat_long"]])
# names(ab_list[["dat_meta"]])
# names(ab_list[["dat_abundances"]])
# ab_list[["dat_meta"]]$ms_meta
# ab_list[["dat_abundances"]]$ms_abundances_eq_med_norm
# ab_list[["dat_abundances"]]$ms_abundances
# ab_list[["dat_meta"]][["norm_ab_meta.eq_med_norm"]]
# ab_list[["dat_meta"]][["pea_meta.eq_med_norm"]]
###################################


# Load normalised data
# norm_ms <- load("./RData/002a_norm_ms.RData")
# norm_pea <- load("./RData/002b_norm_pea.RData")

#norm_ab_transformations.0.0


# load list (not used)
# ms_pea_data <- load("./RData/003_ab_list.RData")
# str(ms_pea_data)
# names(ab_list)

# ab_list[["dat_abundances"]]
# ab_list[["dat_meta"]]


# dim(norm_ab_transformations.0.0)
# dim(npx_transformations.0.0)

# Load test results

# read ms DEqMS output file
DEqMS_out <- readxl::read_excel("./out_r/007_DEqMS_res.xlsx" ,
                                    #sheet = "",
                                    col_names = TRUE,
                                    ) 

# read pea ttest output file
t_test_pea <- readxl::read_excel("./out_r/007_t_test_res_pea.xlsx" ,
                            #sheet = "",
                            col_names = TRUE,
) 

# join and order colunms
# ms
ms_filt_norm_test <- inner_join(ab_list[["dat_meta"]][["norm_ab_meta.eq_med_norm"]],DEqMS_out,by=c("Gene.symbol"="Gene.symbol")) %>% #norm_ab_transformations.0.0
  dplyr::select(UniprotID,Gene.symbol,Gene.synonyms=Gene.synonyms.x, contains("_"),logFC,AveExpr,t,P.Value,adj.P.Val,B,gene,count,sca.t,sca.P.Value,sca.adj.pval,Possible.contamination ) %>% 
  dplyr::select(UniprotID,Gene.symbol,Gene.synonyms,High_1,High_2,High_3,High_4,High_5,High_6,High_7,High_8,Low_1,Low_2,Low_3,Low_4,Low_5,Low_6,Low_7,Low_8,Numb_NA_High,Numb_NA_Low,logFC,AveExpr,t,P.Value,adj.P.Val,B,gene,count,sca.t,sca.P.Value,sca.adj.pval,Possible.contamination)
#cat(names(ms_filt_norm_test),sep=",")
# pea
pea_filt_norm_test <- inner_join(ab_list[["dat_meta"]][["pea_meta.eq_med_norm"]],t_test_pea ,by=c("Gene.symbol"="Gene.symbol")) %>% #npx_transformations.0.0
  dplyr::select(UniprotID=UniprotID.x,OlinkID,Gene.symbol, contains("_"),estimate,estimate1,estimate2,statistic,p.value,parameter,conf.low,conf.high,method,alternative,BH_adj,fdr_adj,Panel,Possible.contamination)
#names(pea_filt_norm_test)

save(ms_filt_norm_test,pea_filt_norm_test,file="./RData/12_formatted_out.RData")
```
