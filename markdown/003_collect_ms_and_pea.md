003\_collect\_ms\_and\_pea.R
================
xrydbh
2022-02-04

``` r
library(tidyverse)
library(DEqMS)
```

    ## Loading required package: limma

    ## 
    ## Attaching package: 'limma'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     plotMA

``` r
# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).



# Get abundances

ms_dat <- load("./RData/002a_norm_ms.RData")
ms_dat
```

    ## [1] "norm_ab_transformations_long" "norm_ab_transformations.0.0"

``` r
pea_dat <- load("./RData/002b_norm_pea.RData")
pea_dat
```

    ## [1] "npx_transformations_long" "npx_transformations.0.0"

``` r
ms_abundances <-  norm_ab_transformations.0.0 %>% 
  dplyr::select(-c("Gene.symbol","Gene.synonyms","Numb_NA_High","Numb_NA_Low","No.peptides","PSMs","Unique.Peptides","MW.kDa"))

pea_abundances <-  npx_transformations.0.0 %>% 
  ungroup() %>% 
  dplyr::select(-c("Gene.symbol","OlinkID","Numb_NA_High","Numb_NA_Low"))


##################
# Equal median normalization is performed in DEqMS package
##################
ms_ab_mat <- ms_abundances[,-1] %>% # Make numeric
  mutate_all(as.numeric) %>% 
  as.matrix()
rownames(ms_ab_mat) <- ms_abundances$UniprotID

boxplot(ms_ab_mat,las=2,main="MS") # check distribution
```

![](/Users/xrydbh/Personal/Projects/git_cloned/Neored_fresh/Neored/markdown/003_collect_ms_and_pea_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# MS
ms_ab_mat = equalMedianNormalization(ms_ab_mat) # normalise

ms_abundances.eq_med_norm.0.0 <- ms_ab_mat %>% # mak tibbl
  as_tibble()

ms_abundances.eq_med_norm.0.1 <- row.names(ms_ab_mat ) %>% # Put pack rownames as column 1 and also insert Gene symobl as col 2
  bind_cols(norm_ab_transformations.0.0$Gene.symbol,ms_abundances.eq_med_norm.0.0) %>% 
  dplyr::rename(UniprotID=...1,Gene.symbol=...2) %>% 
  mutate_at(vars(-c("UniprotID","Gene.symbol")),as.character)
```

    ## New names:
    ## * NA -> ...1
    ## * NA -> ...2

``` r
# Make long version of Prot_median_cent_l2_norm_ab.eq_med_norm
ms_abs_long.eq_med_norm.0.0 <- ms_abundances.eq_med_norm.0.1 %>%
  pivot_longer(cols=-c("UniprotID","Gene.symbol")) %>% 
  mutate_at(vars(c("value")),as.numeric) %>% 
  separate(
    name,
    c("Stem.cell.cont","Sample"),
    sep = "_")

# ms_for_bp.0.0 <- ms_abs_long.eq_med_norm.0.0 %>% 
#   dplyr::filter(!!as.symbol(protID)%in%top_ms)



# PEA
pea_ab_mat <- pea_abundances %>% 
  dplyr::select(-c("UniprotID","Panel")) %>% 
  mutate_all(as.numeric) %>% 
  as.matrix()
rownames(pea_ab_mat) <- pea_abundances$UniprotID
boxplot(pea_ab_mat,las=2,main="PEA")
```

![](/Users/xrydbh/Personal/Projects/git_cloned/Neored_fresh/Neored/markdown/003_collect_ms_and_pea_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

``` r
pea_ab_mat = equalMedianNormalization(pea_ab_mat) # normalise

boxplot(pea_ab_mat,las=2,main="PEA")
```

![](/Users/xrydbh/Personal/Projects/git_cloned/Neored_fresh/Neored/markdown/003_collect_ms_and_pea_files/figure-gfm/unnamed-chunk-1-3.png)<!-- -->

``` r
pea_abundances.eq_med_norm.0.0 <- pea_ab_mat %>% # make tibbl
  as_tibble()

pea_abundances.eq_med_norm.0.1 <- row.names(pea_ab_mat ) %>% # Put pack rownames as column 1 and also insert Gene symobl as col 2
  bind_cols(npx_transformations.0.0$Gene.symbol,pea_abundances.eq_med_norm.0.0,npx_transformations.0.0$Panel) %>% 
  dplyr::rename(UniprotID=...1,Gene.symbol=...2,Panel=...19) %>% 
  mutate_at(vars(-c("UniprotID","Gene.symbol")),as.character)
```

    ## New names:
    ## * NA -> ...1
    ## * NA -> ...2
    ## * NA -> ...19

``` r
# Make long version of Prot_median_cent_l2_norm_ab.eq_med_norm
pea_abs_long.eq_med_norm.0.0 <- pea_abundances.eq_med_norm.0.1 %>%
  pivot_longer(cols=-c("UniprotID","Gene.symbol","Panel")) %>% 
  mutate_at(vars(c("value")),as.numeric) %>% 
  separate(
    name,
    c("Stem.cell.cont","Sample"),
    sep = "_")
################### Fix Meta
# MS
ms_abundances.eq_med_norm.0.2 <- ms_abundances.eq_med_norm.0.1 %>% 
  dplyr::select(-Gene.symbol)
norm_ab_transformations.0.1 <- norm_ab_transformations.0.0 %>% 
  dplyr::select(-contains("High_"))%>% 
  dplyr::select(-contains("Low_"))
norm_ab_meta.eq_med_norm.0.0 <- norm_ab_transformations.0.1 %>% inner_join(ms_abundances.eq_med_norm.0.2,by=c("UniprotID"="UniprotID"))

# PEA
pea_abundances.eq_med_norm.0.2 <- pea_abundances.eq_med_norm.0.1 %>% 
  dplyr::select(-c(Gene.symbol,Panel))
npx_transformations.0.1 <- npx_transformations.0.0 %>% 
  dplyr::select(-contains("High_"))%>% 
  dplyr::select(-contains("Low_"))
npx_meta.eq_med_norm.0.0 <- npx_transformations.0.1 %>% inner_join(pea_abundances.eq_med_norm.0.2,by=c("UniprotID"="UniprotID"))
#####################


dat_long <- list(ms_abs_long_eq_med_norm=ms_abs_long.eq_med_norm.0.0,pea_abs_long_eq_med_norm=pea_abs_long.eq_med_norm.0.0,pea_long=npx_transformations_long,ms_long=norm_ab_transformations_long)
dat_meta <- list(norm_ab_meta.eq_med_norm=norm_ab_meta.eq_med_norm.0.0,ms_meta=norm_ab_transformations.0.0,pea_meta.eq_med_norm=npx_meta.eq_med_norm.0.0 ,pea_meta=npx_transformations.0.0)
dat_abundances <- list(ms_abundances_eq_med_norm=ms_abundances.eq_med_norm.0.1,pea_abundances_eq_med_norm=pea_abundances.eq_med_norm.0.1 ,pea_abundances=pea_abundances,ms_abundances=ms_abundances)
dat_abund_test <- list(ms_abundances=ms_abundances[1:15,c(1:4,10:12)],pea_abundances=pea_abundances[1:15,c(1:4,10:12)])
ab_list <- list(dat_abundances=dat_abundances,dat_long=dat_long,dat_meta=dat_meta,dat_abund_test=dat_abund_test)


save(ab_list,file="./RData/003_ab_list.RData")


(dim(dat_abund_test[[2]])[2]-1)/2
```

    ## [1] 3

``` r
rep_times <- (dim(dat_abund_test[[2]])[2]-1)/2
cond = c(rep("High",rep_times))
```
