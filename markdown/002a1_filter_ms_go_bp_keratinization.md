002a1\_filter\_ms\_go\_bp\_keratinization.R
================
xrydbh
2022-02-04

``` r
library(tidyverse)
library(readxl)
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

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

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

    ## The following object is masked from 'package:Matrix':
    ## 
    ##     expand

    ## The following object is masked from 'package:rlist':
    ## 
    ##     List

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:sjmisc':
    ## 
    ##     trim

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce

    ## 
    ## Attaching package: 'AnnotationDbi'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

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

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).


# remove old file



# read ms file
loaded_objs <- load("./RData/000_ms_numb_unique_prot_IDs.RData")

#### Remove Keratinization related
ms_list <- ms.raw.0.0$UniprotID

length(ms_list)
```

    ## [1] 1419

``` r
#Annotated proteisn with AnnotationDbi
ms_terms.long <- AnnotationDbi::select(Homo.sapiens,
                                       keys=ms_list,
                                       columns = c("TERM","ALIAS"),     #,"ONTOLOGY","DEFINITION"
                                       keytype = "UNIPROT") %>% 
  as_tibble()
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
ms_terms.long$UNIPROT %>% unique() %>% length()
```

    ## [1] 1419

``` r
# Rename variables of output and filter out only entries with BP annotation
# OBS if proteins has no BO annotation it will be lost
ms_terms.bp_only.wide <- ms_terms.long %>% 
  dplyr::rename(Gene.symbol=ALIAS) %>% 
  dplyr::rename(UniprotID=UNIPROT) %>% 
  filter(ONTOLOGY=="BP") %>% 
  dplyr::select(UniprotID,TERM) %>% 
  distinct() %>% 
  dplyr::group_by(UniprotID) %>% 
  dplyr::summarise(GO_term = paste(TERM, collapse = ",")) %>% 
  dplyr::rename(GO_BP_terms=GO_term)


# Filter in "K/eratin" associated, but not inflamm etc assoc
keratin_associated.0.0 <- ms_terms.bp_only.wide %>% 
  filter(grepl("Keratin|keratin", GO_BP_terms)) %>%  # remove any proteins with K/keratin in GO terms
  filter(!(grepl("inflamm|immune|hemopoies|hematopoie|growth factor|cytokine", GO_BP_terms))) %>% # but keep if contains any of word parts inflamm|immune|hemopoies|hematopoie|growth factor|cytokine
  filter(UniprotID!="Q02487") # we want to leave this one out since it exists in the Olink targeted platform


# # Filter in  associated to relevant words
# keratin_related_to_keep <- keratin_associated.0.0 %>% 
#   filter(!(grepl("inflamm|immune|hemopoies|hematopoie|growth factor|cytokine", GO_BP_terms)))

# Store in  variable the Keratins related that will be kept
keratin_related_to_keep.0.0 <- ms_terms.bp_only.wide %>% 
  filter(grepl("Keratin|keratin", GO_BP_terms)) %>%  # remove any proteins with K/keratin in GO terms
  filter((grepl("inflamm|immune|hemopoies|hematopoie|growth factor|cytokine", GO_BP_terms))) 


# Join keratin associated to discard  with ms.raw to get gene symobls
keratin_associated.0.1 <- ms.raw.0.0 %>% 
  dplyr::select(UniprotID,Gene.symbol,Gene.synonyms) %>% 
  inner_join(keratin_associated.0.0,by=c("UniprotID"="UniprotID"))


# Join keratin associated to discard  with ms.raw to get gene symobls
keratin_related_to_keep.0.1 <- ms.raw.0.0 %>% 
  dplyr::select(UniprotID,Gene.symbol,Gene.synonyms) %>% 
  inner_join(keratin_related_to_keep.0.0,by=c("UniprotID"="UniprotID"))

# ms_unip_without_keratin_ass  <- ms_terms.bp_only.wide %>% 
#   filter(!grepl("eratin", GO_BP_terms)) %>% 
#   pluck("UniprotID") %>% 
#   unique()
# 
# length(ms_unip_without_keratin_ass %>% unique())
# 
# 
# keratin_associated.0.0$UniprotID

# This will list proteins that do not have keratins in their name but also 
# proteins with no BP annotation
# ms.raw.0.1.alt <- ms.raw.0.0 %>% 
#   dplyr::filter(UniprotID%in%ms_unip_without_keratin_ass)

# Use keratin_associated to be discarde to filter ms.raw
ms.raw.0.1 <- ms.raw.0.0 %>% 
  dplyr::filter(!(UniprotID%in%keratin_associated.0.0$UniprotID))



# keratins <- ms.raw.0.0 %>% 
#   dplyr::filter(UniprotID%in%keratin_associated)

############################
# Remove proteins listed as contamination in "Proteomics study of human cord blood reticulocyte-derived exosomes" Varela et.al (DOI:10.1038/s41598-018-32386-2)
# Described/motivated in "The human plasma proteome: history, character, and diagnostic prospects"; (DOI: 10.1074/mcp.r200007-mcp200)
# This was done to focus on Exosome proteins. Should not be done in this analyis.
############################

# removed_in_vareal_et_al.0.0 <- readxl::read_excel("./data/potential_contaminants_diaz_varela.xlsx",sheet="potential contaminants removed",skip=3) 
# removed_uniprot_varela <- removed_in_vareal_et_al.0.0$Accession
# 
# ms.raw.0.2 <- ms.raw.0.1 %>% 
#   dplyr::filter(!(UniprotID%in%removed_uniprot_varela))

############################
# Save raw
############################
fn = paste("./data/001_ms_raw.xlsx",sep="")
ms.raw.0.1 %>%
  writexl::write_xlsx(path = fn)

# Save keratin associatied discarded, with some uniprot IDs filled in
fn = paste("./data/001_ms_filt_out_keratin_associated.xlsx",sep="")
keratin_associated.0.1  %>%
  writexl::write_xlsx(path = fn)

# save keratin_related_to_keep
fn = paste("./data/001_ms_filt_in_keratin_associated.xlsx",sep="")
keratin_related_to_keep.0.1  %>%
  writexl::write_xlsx(path = fn)
```
