013\_add\_annotation.R
================
xrydbh
2022-02-04

``` r
library(AnnotationDbi)
```

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: IRanges

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

``` r
library("Homo.sapiens")
```

    ## Loading required package: OrganismDbi

    ## Loading required package: GenomicFeatures

    ## Loading required package: GenomeInfoDb

    ## Loading required package: GenomicRanges

    ## Loading required package: GO.db

    ## 

    ## Loading required package: org.Hs.eg.db

    ## 

    ## Loading required package: TxDb.Hsapiens.UCSC.hg19.knownGene

``` r
library(GO.db)
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
    ## ✓ tibble  3.1.6     ✓ dplyr   1.0.7
    ## ✓ tidyr   1.1.4     ✓ stringr 1.4.0
    ## ✓ readr   2.1.1     ✓ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::collapse()   masks IRanges::collapse()
    ## x dplyr::combine()    masks Biobase::combine(), BiocGenerics::combine()
    ## x dplyr::desc()       masks IRanges::desc()
    ## x tidyr::expand()     masks S4Vectors::expand()
    ## x dplyr::filter()     masks stats::filter()
    ## x dplyr::first()      masks S4Vectors::first()
    ## x dplyr::lag()        masks stats::lag()
    ## x ggplot2::Position() masks BiocGenerics::Position(), base::Position()
    ## x purrr::reduce()     masks GenomicRanges::reduce(), IRanges::reduce()
    ## x dplyr::rename()     masks S4Vectors::rename()
    ## x dplyr::select()     masks OrganismDbi::select(), AnnotationDbi::select()
    ## x dplyr::slice()      masks IRanges::slice()

``` r
# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).



my_loaded_objs <- load("./RData/12_formatted_out.RData")
my_loaded_objs
```

    ## [1] "ms_filt_norm_test"  "pea_filt_norm_test"

``` r
# Read WGCNA res
wgcna.0.0 <- readxl::read_excel("./out_r/WGCNA/pea_black_module_genes_with_go.xlsx") %>% 
  dplyr::select(UniprotID) %>% 
  mutate(inWGCNAmodule_yesIs1=1)



# Use the columns functinon of the AnnotaionDB to list the columns of the Homosapiens
columns(Homo.sapiens)
```

    ##  [1] "ACCNUM"       "ALIAS"        "CDSCHROM"     "CDSEND"       "CDSID"       
    ##  [6] "CDSNAME"      "CDSSTART"     "CDSSTRAND"    "DEFINITION"   "ENSEMBL"     
    ## [11] "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"    
    ## [16] "EVIDENCEALL"  "EXONCHROM"    "EXONEND"      "EXONID"       "EXONNAME"    
    ## [21] "EXONRANK"     "EXONSTART"    "EXONSTRAND"   "GENEID"       "GENENAME"    
    ## [26] "GO"           "GOALL"        "GOID"         "IPI"          "MAP"         
    ## [31] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"        
    ## [36] "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "TERM"        
    ## [41] "TXCHROM"      "TXEND"        "TXID"         "TXNAME"       "TXSTART"     
    ## [46] "TXSTRAND"     "TXTYPE"       "UCSCKG"       "UNIGENE"      "UNIPROT"

``` r
ms_list <- ms_filt_norm_test$UniprotID
ms_list_short <- ms_filt_norm_test$UniprotID[1:5]
pea_list <- pea_filt_norm_test$UniprotID

ms_terms.long <- AnnotationDbi::select(Homo.sapiens,
                      keys=ms_list,
                      columns = c("TERM","ALIAS"),     #,"ONTOLOGY","DEFINITION"
                      keytype = "UNIPROT",
                      multiVals="first") %>% 
  as_tibble()
```

    ## 'select()' returned many:many mapping between keys and columns

``` r
ms_terms.bp_only.wide <- ms_terms.long %>% 
  dplyr::rename(Gene.symbol=ALIAS) %>% 
  dplyr::rename(UniprotID=UNIPROT) %>% 
  
  filter(ONTOLOGY=="BP") %>% 
  dplyr::select(UniprotID,TERM) %>% 
  distinct() %>% 
  dplyr::group_by(UniprotID) %>% 
  dplyr::summarise(GO_term = paste(TERM, collapse = ",")) %>% 
  dplyr::rename(GO_BP_terms=GO_term)

ms_filt_norm_test.go_anno <- ms_filt_norm_test %>% 
  full_join(ms_terms.bp_only.wide,by=c("UniprotID"="UniprotID"))

# PEA

pea_terms.long <- AnnotationDbi::select(Homo.sapiens,
                                  keys=pea_list,
                                  columns = c("TERM","ALIAS"),     #,"ONTOLOGY","DEFINITION"
                                  keytype = "UNIPROT") %>% 
  as_tibble()
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
pea_terms.bp_only.wide <- pea_terms.long %>% 
  dplyr::rename(Gene.symbol=ALIAS) %>% 
  dplyr::rename(UniprotID=UNIPROT) %>% 
  filter(ONTOLOGY=="BP") %>% 
  dplyr::select(UniprotID,TERM) %>% 
  distinct() %>% 
  dplyr::group_by(UniprotID) %>% 
  dplyr::summarise(GO_term = paste(TERM, collapse = ",")) %>% 
  dplyr::rename(GO_BP_terms=GO_term)

pea_filt_norm_test.go_anno  <- pea_filt_norm_test %>% 
  full_join(pea_terms.bp_only.wide,by=c("UniprotID"="UniprotID"))%>% 
  full_join(wgcna.0.0,by=c("UniprotID"="UniprotID")) %>% 
  mutate(inWGCNAmodule_yesIs1=tidyr::replace_na(inWGCNAmodule_yesIs1, 0))


# write formatted output to file
ms_filt_norm_test.go_anno  %>%
  distinct() %>% 
  writexl::write_xlsx(path = "./out_r/013_ms_filt_norm_test_go_anno.xlsx")

pea_filt_norm_test.go_anno  %>%
  writexl::write_xlsx(path = "./out_r/013_pea_filt_norm_test_go_anno.xlsx")

# write only GO to file



ms_terms.long  %>%
  distinct() %>% 
  writexl::write_xlsx(path = "./out_r/013_ms_go_terms.xlsx")
pea_terms.long  %>%
  writexl::write_xlsx(path = "./out_r/013_pea_go_terms.xlsx")
```
