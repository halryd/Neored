library(tidyverse)

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).



# read ms DEqMS output file
res_S4 <- readxl::read_excel("./out_r/020_Cord_blood_proteomics_CD34_Supplementary_tables_2020_06_15.xlsx" ,
                             sheet = "S4 MS results camera reactome",
                             col_names = TRUE,
)

res_S2 <- readxl::read_excel("./out_r/020_Cord_blood_proteomics_CD34_Supplementary_tables_2020_06_15.xlsx" ,
                             sheet = "S2 MS data and res, per prot",
                             col_names = TRUE,
)
# Get abundances
loaded_abundancies <- load("./RData/003_ab_list.RData")
loaded_abundancies 


#pw11 Translation
protsOfTranslation <- res_S4$Submitted.entities.found[11] %>% 
  str_replace_all(";", ",") %>% 
  str_split(",") %>% 
  unlist

pwOfprotsOfTranslation <- res_S2 %>% 
  filter(UniprotID%in%protsOfPathW )

pwOfprotsOfTranslation  %>%
  writexl::write_xlsx(path = "./014_pwOfprotsOfTranslation.xlsx")
