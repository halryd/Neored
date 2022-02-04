020\_compile\_supplementory\_tables.R
================
xrydbh
2022-02-04

``` r
library(tidyverse)
library(readxl)
library(writexl)

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).



#setwd("../")
# read keratins
keratins.0.0 <- readxl::read_excel("./data/001_ms_filt_out_keratin_associated.xlsx")

# read protein centric
ms.deqms.0.0 <- readxl::read_excel("./out_r/013_ms_filt_norm_test_go_anno.xlsx") 
pea.t.test.0.0 <- readxl::read_excel("./out_r/013_pea_filt_norm_test_go_anno.xlsx") 
pea.t.test.0.1 <- pea.t.test.0.0 %>% 
  dplyr::select(UniprotID,T_test_p_val=p.value,Foldchange=estimate)
# read_pathway_centric
ms_camera_reactome.0.0 <- readxl::read_excel("./out_r/Reactome/reactome_results.xlsx")

# read pw per prot
# No significant pathway for pea
pw_per_prot.0.0 <- readxl::read_excel("./out_r/014_pw_per_prot.xlsx")


# read correlation module centric
# No  module significantly correlatied to CD34 conncentration

pea_wgcna.0.0 <- readxl::read_excel("./out_r/WGCNA/pea_black_module_genes_with_go.xlsx")
pea_wgcna.0.1 <- pea_wgcna.0.0 %>% 
  left_join(pea.t.test.0.1) %>% 
  relocate(T_test_p_val,Foldchange,.before=GO_BP_terms)
```

    ## Joining, by = "UniprotID"

``` r
# Make list

supplementary.0.0 <- list()
supplementary.0.0[["S1 MS filtered out keratins"]] <- keratins.0.0
supplementary.0.0[["S2 MS data and res, per prot"]] <- ms.deqms.0.0
supplementary.0.0[["S3 PEA data and res, per prot"]] <- pea.t.test.0.0 
supplementary.0.0[["S4 MS results camera reactome"]] <- ms_camera_reactome.0.0
supplementary.0.0[["S5 MS Numb pathw per protein"]] <- pw_per_prot.0.0
supplementary.0.0[["S6 PEA results WGCNA"]] <- pea_wgcna.0.1

supplementary.0.0  %>%
  writexl::write_xlsx(path = "./out_r/020_Cord_blood_proteomics_CD34_Supplementary_tables_2020_06_15.xlsx")

#####
```
