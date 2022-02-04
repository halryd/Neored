library(clusterProfiler)
library(pathview)
library(enrichplot)
library(annotables)
library(tidyverse)
library(readxl)
library(org.Hs.eg.db)

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).



# WGCNA Module
fn <- "./out_r/WGCNA/pea/ResultsWGCNA.xlsx"
module_of_int <- "black"
wgcna_module.0.0 <- readxl::read_excel(fn,sheet=module_of_int)


#############Get genes of module of interest and convert to Entrez 

# Get Gene.symbols
myProts <- wgcna_module.0.0 %>% 
  pluck(1)

# Get module membership
kME <- wgcna_module.0.0 %>% 
  pluck("kMEblack")

# Assigne names
names(kME) <- myProts

## Return the IDs for the gene symbols in the DE results from Genome Reference Consortium Human genome build 37; comes with package annotables
idx <- grch37$symbol %in% myProts # Generates a logical vector
ids <- grch37[idx, ] # Use logical vector to subset grch37

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
non_duplicates <- which(duplicated(ids$symbol) == FALSE)

ids <- ids[non_duplicates, ] 

## Merge the IDs with the results 
res_ids <- inner_join(wgcna_module.0.0, ids, by=c("Gene.symbol"="symbol"))

## Remove any NA values
res_entrez <- dplyr::filter(res_ids, entrez != "NA")

## Remove any Entrez duplicates
res_entrez <- res_entrez[which(duplicated(res_entrez$entrez) == F), ]

## Extract significant rows
pea_sign.0.0  <- res_entrez 

## Extract significant IDs
pea_sign_entrez <- pea_sign.0.0$entrez

pea_sign_entrez

#############Get all anlayzed genes (background) and convert to Entrez 

wgcna_background.0.0 <- readxl::read_excel(fn,sheet="AllData")

allProts <- wgcna_background.0.0  %>% 
  pluck("Gene.symbol")

## Return the IDs for the gene symbols in the DE results
idx <- grch37$symbol %in% allProts

ids <- grch37[idx, ]

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
non_duplicates <- which(duplicated(ids$symbol) == FALSE)

ids <- ids[non_duplicates, ] 

## Merge the IDs with the results 
res_ids <- inner_join(wgcna_background.0.0, ids, by=c("Gene.symbol"="symbol"))

## Remove any NA values
res_entrez <- dplyr::filter(res_ids, entrez != "NA")

## Remove any Entrez duplicates
res_entrez <- res_entrez[which(duplicated(res_entrez$entrez) == F), ]

## Extract significant rows
pea_all.0.0  <- res_entrez 

## Extract significant IDs
pea_all_entrez <- pea_all.0.0$entrez

pea_back_entrez <- res_entrez$entrez

pea_back_entrez[is.na(as.character(pea_back_entrez))]

#######Get hub genes of module of interest and convert to Entrez 
wgcna_module.0.0 %>% 
  dplyr::select(Gene.symbol,kMEblack)




ora_go_bp <- enrichGO(gene         = as.character(pea_sign_entrez),
                      universe      = as.character(pea_back_entrez),
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

head(ora_go_bp)

# pdf(file = paste("./out_r/WGCNA/pea/clusterProfiler_goplot_ora_black_module_go_bp.pdf",sep=""),
#     width = 15 , 
#     height= 10)
# goplot(ora_go_bp)
# dev.off()

ora_kegg <- enrichKEGG(gene = as.character(pea_sign_entrez),
                       universe = as.character(pea_back_entrez),
                       organism  = 'hsa',
                       minGSSize = 5,
                       maxGSSize = 500,
                       pvalueCutoff = 0.05)
ora_kegg  %>%
  as_tibble() %>% 
  writexl::write_xlsx(path = "./out_r/WGCNA/pea/ora_kegg.xlsx")

entresIDsInPw.0.0 <- c("1509","8763","5641","1200")
names(grch37)
entresIDsInPw.0.1 <- grch37$entrez %in% entresIDsInPw.0.0
entresIDsInPw.0.2 <- grch37[entresIDsInPw.0.1, ]
gsInPw <- entresIDsInPw.0.2$symbol

#temp1509/8763/5641/1200


quality_markers.0.0 <- readxl::read_excel("./data/quality_markers.xlsx")
quality_markers.0.0 %>% 
  filter(`Gene names`%in%gsInPw)# c("CAT","CA2")

# Visualize enriched KEGG pathways
browseKEGG(ora_kegg, 'hsa04142')

########################
# Get Gene.symbols
myProts <- pea_all.0.0 %>% 
  pluck("entrez")

# Get module membership
kME <- pea_all.0.0 %>% 
  pluck("kMEblack")

# Assigne names
names(kME) <- myProts


pea_all.0.0 %>% 
  dplyr::select(Gene.symbol,entrez,kMEblack)
wgcna_module.0.0
kME.2 <- kME %>% sort(decreasing = T)

gsea_go_cc <- gseGO(geneList = kME.2,
                    OrgDb        = org.Hs.eg.db,
                    ont          = "MF",
                    minGSSize    = 10,
                    maxGSSize    = 500,
                    pvalueCutoff = 0.1,
                    verbose      = FALSE)
## Extract the GSEA results
gsea_go_cc_results <- gsea_go_cc@result
gsea_go_cc_results
