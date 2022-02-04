011\_WGCNA\_ORA\_Kegg.R
================
xrydbh
2022-02-04

``` r
library(clusterProfiler)
library(pathview)
```

    ## ##############################################################################
    ## Pathview is an open source software package distributed under GNU General
    ## Public License version 3 (GPLv3). Details of GPLv3 is available at
    ## http://www.gnu.org/licenses/gpl-3.0.html. Particullary, users are required to
    ## formally cite the original Pathview paper (not just mention it) in publications
    ## or products. For details, do citation("pathview") within R.
    ## 
    ## The pathview downloads and uses KEGG data. Non-academic uses may require a KEGG
    ## license agreement (details at http://www.kegg.jp/kegg/legal.html).
    ## ##############################################################################

``` r
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
```

    ##  [1]  7064  2288  3068  9332  2811  1509  7917  8482  2167  2220 10859  8763
    ## [13]  3952  5641 30817  1200 10288

``` r
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
```

    ## integer(0)

``` r
#######Get hub genes of module of interest and convert to Entrez 
wgcna_module.0.0 %>% 
  dplyr::select(Gene.symbol,kMEblack)
```

    ## # A tibble: 23 × 2
    ##    Gene.symbol   kMEblack
    ##    <chr>            <dbl>
    ##  1 THOP1            0.928
    ##  2 FKBP4            0.912
    ##  3 BLM_hydrolase    0.879
    ##  4 HDGF             0.845
    ##  5 CD163            0.839
    ##  6 GP1BA            0.833
    ##  7 CTSD             0.826
    ##  8 ICAM-2           0.801
    ##  9 BAG6             0.784
    ## 10 SEMA7A           0.780
    ## # … with 13 more rows

``` r
ora_go_bp <- enrichGO(gene         = as.character(pea_sign_entrez),
                      universe      = as.character(pea_back_entrez),
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

head(ora_go_bp)
```

    ## [1] ID          Description GeneRatio   BgRatio     pvalue      p.adjust   
    ## [7] qvalue      geneID      Count      
    ## <0 rows> (or 0-length row.names)

``` r
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
```

    ## Reading KEGG annotation online:

    ## Reading KEGG annotation online:

``` r
ora_kegg  %>%
  as_tibble() %>% 
  writexl::write_xlsx(path = "./out_r/WGCNA/pea/ora_kegg.xlsx")

entresIDsInPw.0.0 <- c("1509","8763","5641","1200")
names(grch37)
```

    ## [1] "ensgene"     "entrez"      "symbol"      "chr"         "start"      
    ## [6] "end"         "strand"      "biotype"     "description"

``` r
entresIDsInPw.0.1 <- grch37$entrez %in% entresIDsInPw.0.0
entresIDsInPw.0.2 <- grch37[entresIDsInPw.0.1, ]
gsInPw <- entresIDsInPw.0.2$symbol

#temp1509/8763/5641/1200


quality_markers.0.0 <- readxl::read_excel("./data/quality_markers.xlsx")
quality_markers.0.0 %>% 
  filter(`Gene names`%in%gsInPw)# c("CAT","CA2")
```

    ## # A tibble: 0 × 8
    ## # … with 8 variables: Protein names <chr>, Gene names <chr>, Protein IDs <chr>,
    ## #   Majority protein IDs <chr>, Median Intensity Erythrocytes [Log10] <dbl>,
    ## #   Coefficient of Variation Erythrocytes [%] <dbl>,
    ## #   Median Intensity Plasma [Log10] <chr>,
    ## #   Ratio Erythrocytes / Plasma [Log10] <chr>

``` r
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
```

    ## # A tibble: 332 × 3
    ##    Gene.symbol entrez kMEblack
    ##    <chr>        <int>    <dbl>
    ##  1 CA1            759   0.511 
    ##  2 ICAM1         3383   0.616 
    ##  3 CHL1         10752  -0.303 
    ##  4 TGFBI         7045   0.275 
    ##  5 ENG           2022  -0.109 
    ##  6 SERPINA7      6906   0.637 
    ##  7 IGFBP3        3486  -0.140 
    ##  8 CR2           1380  -0.0971
    ##  9 SERPINA5      5104  -0.0108
    ## 10 FCGR3B        2215   0.198 
    ## # … with 322 more rows

``` r
wgcna_module.0.0
```

    ## # A tibble: 23 × 19
    ##    Gene.symbol   High_1  High_2 High_3  High_4  High_5  High_6 High_7  High_8
    ##    <chr>          <dbl>   <dbl>  <dbl>   <dbl>   <dbl>   <dbl>  <dbl>   <dbl>
    ##  1 THOP1          1.40  0.793   -0.658 -0.0796  0.182   0.330   1.06   0.345 
    ##  2 FKBP4          2.44  2.50    -1.23   0.704   1.13    1.58    1.18   0.191 
    ##  3 BLM_hydrolase  0.633 0.553   -0.311 -0.265   0.279   0.280   0.712  0.453 
    ##  4 HDGF           2.14  0.969   -1.14   0.832   0.136   1.16   -0.104  0.104 
    ##  5 CD163          0.902 0.922   -0.100 -0.257   0.309   0.249   1.08   0.438 
    ##  6 GP1BA          0.125 0.103   -0.568  0.0115 -0.0913  0.285   0.129 -0.0115
    ##  7 CTSD           1.31  1.21    -0.245 -0.408   0.679  -0.0346  0.347  0.375 
    ##  8 ICAM-2         0.285 0.00176 -0.269 -0.286   0.0702  0.250   0.422  0.115 
    ##  9 BAG6           0.995 0.746   -0.257  0.170   0.509   0.0401  0.568  0.481 
    ## 10 SEMA7A         0.558 0.110   -0.640 -0.516  -0.110   0.722   0.203  0.147 
    ## # … with 13 more rows, and 10 more variables: Low_1 <dbl>, Low_2 <dbl>,
    ## #   Low_3 <dbl>, Low_4 <dbl>, Low_5 <dbl>, Low_6 <dbl>, Low_7 <dbl>,
    ## #   Low_8 <dbl>, moduleColors <chr>, kMEblack <dbl>

``` r
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
```

    ##                    ID                                 Description setSize
    ## GO:0033218 GO:0033218                               amide binding      12
    ## GO:0003723 GO:0003723                                 RNA binding      15
    ## GO:0003676 GO:0003676                        nucleic acid binding      23
    ## GO:0042277 GO:0042277                             peptide binding      11
    ## GO:0005201 GO:0005201 extracellular matrix structural constituent      24
    ##            enrichmentScore       NES       pvalue   p.adjust    qvalues rank
    ## GO:0033218       0.7074982  1.906485 0.0008555936 0.07509174 0.07425339   86
    ## GO:0003723       0.6662499  1.901679 0.0015170048 0.07509174 0.07425339   69
    ## GO:0003676       0.5625160  1.805767 0.0025515832 0.07838646 0.07751133   69
    ## GO:0042277       0.6933931  1.824336 0.0031671296 0.07838646 0.07751133   86
    ## GO:0005201      -0.4915805 -1.843957 0.0039771780 0.07874812 0.07786896   94
    ##                              leading_edge
    ## GO:0033218 tags=75%, list=26%, signal=58%
    ## GO:0003723 tags=67%, list=21%, signal=55%
    ## GO:0003676 tags=52%, list=21%, signal=44%
    ## GO:0042277 tags=73%, list=26%, signal=56%
    ## GO:0005201 tags=58%, list=28%, signal=45%
    ##                                                                  core_enrichment
    ## GO:0033218                         7064/2288/1200/30835/5479/10288/911/8685/1471
    ## GO:0003723                    2288/3068/11164/11315/1476/328/5037/5479/2896/1399
    ## GO:0003676          2288/3068/3952/11164/5168/11315/1476/328/5037/5479/2896/1399
    ## GO:0042277                              7064/1200/30835/5479/10288/911/8685/1471
    ## GO:0005201 920/3371/4240/2202/1634/1311/5118/3910/4147/7148/7058/10418/1277/8076
