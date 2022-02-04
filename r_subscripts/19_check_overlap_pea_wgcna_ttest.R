library(tidyverse)
library(readxl)
library(writexl)
library(ggvenn)

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).



getwd()
pval_thresh = 0.05

# read long read protein centric
wgcna.0.0 <- readxl::read_excel("./out_r/WGCNA/pea_black_module_genes_with_go.xlsx")
ttest.0.0 <- readxl::read_excel("./out_r/007_t_test_res_pea.xlsx") %>% 
  filter(p.value<0.05) 

pea_wgcna_ttest <- wgcna.0.0 %>% inner_join(ttest.0.0, by=("UniprotID"="UniprotID")) %>% 
  dplyr::select( UniprotID,Gene.symbol,Gene.symbol.synonyms,FoldChange=estimate,GO_BP_terms)

pea_wgcna_ttest  %>%
  writexl::write_xlsx(path = "./out_r/019_pea_wgcna_ttest.xlsx")

wgcna_uni <-wgcna.0.0$UniprotID %>% 
  unique()

ttest_uni <-ttest.0.0$UniprotID %>% 
  unique()

venn_list <- list(wgcna_sign_mod_uni=wgcna_uni,ttest_p005_uni=ttest_uni)


pdf(file = paste("./out_r/019_pea_overlap_wgcna_ttest.pdf",sep=""))
ggvenn(
  venn_list , 
  fill_color = c("#0073C2FF", "#EFC000FF"), #, "#CD534CFF",, "#868686FF"
  stroke_size = 0.5, set_name_size = 4
)
dev.off()
