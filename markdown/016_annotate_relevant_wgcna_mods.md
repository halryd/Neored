016\_annotate\_relevant\_wgcna\_mods.R
================
xrydbh
2022-02-04

``` r
library(readxl) 
library(AnnotationDbi)
library("Homo.sapiens")
library(GO.db)
library(tidyverse)

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).


read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

ms.wgcna.mods <- read_excel_allsheets("./out_r/WGCNA/pea_unip/ResultsWGCNA.xlsx") 
ms.wgcna.mods %>% names()
```

    ## [1] "AllData"   "black"     "blue"      "brown"     "green"     "red"      
    ## [7] "turquoise" "yellow"

``` r
my_mod_df.0.0 <-ms.wgcna.mods$black %>% 
  dplyr::select(UniprotID,kMEblack)

# Extract uniprot IDs from WGCNA output
my_mod_unip <- my_mod_df.0.0$UniprotID

# my_mod_gs BLMH

# gs_unip <- AnnotationDbi::select(Homo.sapiens,
#                                  keys=my_mod_gs,
#                                  columns = c("ALIAS","UNIPROT"), 
#                                  keytype = "ALIAS")

# Annotate Uniprot with GO and genes symbol (with all synonyms)
gs_unip.0.0 <- AnnotationDbi::select(Homo.sapiens,
                                       keys=my_mod_unip,
                                       columns = c("UNIPROT","ALIAS","TERM"),     #,"ONTOLOGY","DEFINITION"
                                       keytype = "UNIPROT") %>% 
  as_tibble() %>% 
  dplyr::select(-c(EVIDENCE,ONTOLOGY)) %>% 
  dplyr::rename(Gene.symbol=ALIAS) %>% 
  dplyr::rename(UniprotID=UNIPROT) %>% 
  dplyr::rename(GO_BP_terms=TERM) %>% 
  dplyr::group_by(UniprotID) %>% 
  dplyr::mutate(GO_BP_terms= paste(GO_BP_terms, collapse = ",")) %>% 
  distinct()
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
# Collect all synonyms into one comma separated sring
gs_unip.0.1  <- gs_unip.0.0  %>% 
  dplyr::select(UniprotID,Gene.symbol) %>% 
  group_by(UniprotID) %>% 
  dplyr::summarise(Gene.symbol.synonyms=paste(Gene.symbol, collapse = ","))
  
  #GO_BP_terms = paste(GO_BP_terms, collapse = ",")

# Remove Gene.symbols gs_unip.0.0 (not summarized by uniprotID)
unip_GO.0.0 <- gs_unip.0.0 %>% 
  dplyr::select(UniprotID,GO_BP_terms) %>% 
  distinct()

list_of_unip_vec_per_path <- strsplit(unip_GO.0.0$GO_BP_terms, split = ",") # Conv first uniprot string to vector
  
# Convert semicolon separated string words  to one word value per row (https://stackoverflow.com/questions/15347282/split-delimited-strings-in-a-column-and-insert-as-new-rows)
unip_GO.0.1 <- data.frame(UniprotID = rep(unip_GO.0.0$UniprotID, sapply(list_of_unip_vec_per_path, length)), 
                          GO_BP_terms = unlist(list_of_unip_vec_per_path)) %>% 
  as_tibble() %>% 
  group_by(UniprotID) %>% 
  distinct() %>% 
  group_by(UniprotID) %>% 
  dplyr::summarise(GO_BP_terms =paste(GO_BP_terms , collapse = ","))

gs_unip_go <- gs_unip.0.1 %>% full_join(unip_GO.0.1)
```

    ## Joining, by = "UniprotID"

``` r
gs_unip_go  %>%
  writexl::write_xlsx(path = "./out_r/WGCNA/pea_black_module_genes_with_go.xlsx")

############### Anotate networkplot
library(qgraph)
library(igraph)
```

    ## 
    ## Attaching package: 'igraph'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     as_data_frame, groups, union

    ## The following objects are masked from 'package:purrr':
    ## 
    ##     compose, simplify

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     crossing

    ## The following object is masked from 'package:tibble':
    ## 
    ##     as_data_frame

    ## The following object is masked from 'package:GenomicRanges':
    ## 
    ##     union

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     union

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     union

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     normalize, path, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

``` r
networkWeightCutOff <- 30
p_val_thresh <- 0.05
module_of_interest <- "VisANTInput-TOMblack.txt"

# Get significantly differentially expressed
sign_prots <- readxl::read_excel("./out_r/007_t_test_res_pea.xlsx") %>% 
  filter(p.value<p_val_thresh) %>% 
  pluck("Gene.symbol")


# row.names(my_mat.0.0)%in%sign_prots

# cat(row.names(my_mat.0.0),sep="\",\"")
# seq(1:8)
# seq(9:16)
#my_groups <- list(klutt=c("CD163","THOP1","HDGF","BLM_hydrolase","BAG6","CTSD"),plutt=c("HO-1","LILRB1","ICAM-2","GP1BA","FCN2","Ep-CAM","FKBP4","LILRB2","SEMA7A","CD164"))
my_groups <- list(singn_in_ttest=c(1:9,11:14),not_sign=c(10,15,16))

#



edglist.0.0 <- read.table(paste("./out_r/WGCNA/pea/",module_of_interest,sep=""), 
                          sep = "" , 
                          header = F,
                          na.strings ="", 
                          stringsAsFactors= F) %>% 
  as.matrix()#, nrows = 100
edgelist.0.1 <- edglist.0.0[,c(1,2,5)]
colnames(edgelist.0.1) <- c("From","To","Weight")
edgelist.0.2 <- edgelist.0.1 %>% 
  as_tibble() %>% 
  arrange(desc(Weight)) %>% 
  head(networkWeightCutOff)
# convert edgelist to ?? matrix
g=graph.data.frame(edgelist.0.2, directed = F)
my_mat.0.0 <- get.adjacency(g, sparse = FALSE, attr='Weight')
class(my_mat.0.0) <- "numeric"

pdf(paste("./out_r/WGCNA/",module_of_interest,"_annotated.pdf",sep=""),
   # width = 20, 
   # height = 20
   ) 
qgraph(my_mat.0.0,
       graph= "default", #cor
       layout = "spring",
       edge.labels = TRUE,
       palette= "ggplot2",
       #minimum=0.14,
       legend=F,
       #bg='lightgrey',
       title.cex=1.2,
       label.cex=2,
       #label.norm="OOOOO",
       title=paste("Hubgene TOM network",module_of_interest),
       groups=my_groups,
       label.scale.equal=T,
       labels = colnames(my_mat.0.0))
```

    ## Warning in qgraph(my_mat.0.0, graph = "default", layout = "spring", edge.labels
    ## = TRUE, : Non-finite weights are omitted

``` r
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
my_mat.0.0 %>% as_tibble()
```

    ## # A tibble: 16 × 16
    ##     CD163  THOP1   HDGF BLM_hydrolase   BAG6   CTSD `HO-1` LILRB1 `ICAM-2`
    ##     <dbl>  <dbl>  <dbl>         <dbl>  <dbl>  <dbl>  <dbl>  <dbl>    <dbl>
    ##  1 NA      0.178 NA             0.188  0.119  0.121 NA     NA       NA    
    ##  2  0.178 NA      0.145         0.167  0.154  0.163  0.134  0.152   NA    
    ##  3 NA      0.145 NA            NA      0.143 NA      0.129 NA       NA    
    ##  4  0.188  0.167 NA            NA      0.121  0.137 NA     NA       NA    
    ##  5  0.119  0.154  0.143         0.121 NA     NA     NA     NA       NA    
    ##  6  0.121  0.163 NA             0.137 NA     NA     NA     NA       NA    
    ##  7 NA      0.134  0.129        NA     NA     NA     NA     NA       NA    
    ##  8 NA      0.152 NA            NA     NA     NA     NA     NA       NA    
    ##  9 NA     NA     NA            NA     NA     NA     NA     NA       NA    
    ## 10 NA      0.118 NA            NA     NA     NA     NA     NA       NA    
    ## 11 NA     NA     NA            NA      0.122 NA     NA     NA       NA    
    ## 12 NA     NA     NA            NA     NA     NA     NA     NA       NA    
    ## 13  0.153  0.185  0.177         0.144  0.166  0.142  0.162 NA       NA    
    ## 14 NA     NA     NA            NA     NA     NA     NA      0.133   NA    
    ## 15 NA     NA     NA            NA     NA     NA     NA     NA        0.131
    ## 16 NA     NA     NA            NA     NA     NA     NA     NA       NA    
    ## # … with 7 more variables: GP1BA <dbl>, FCN2 <dbl>, `Ep-CAM` <dbl>,
    ## #   FKBP4 <dbl>, LILRB2 <dbl>, SEMA7A <dbl>, CD164 <dbl>
