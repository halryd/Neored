014\_doublecheck\_foldchange.R
================
xrydbh
2022-02-04

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
    ## ✓ tibble  3.1.6     ✓ dplyr   1.0.7
    ## ✓ tidyr   1.1.4     ✓ stringr 1.4.0
    ## ✓ readr   2.1.1     ✓ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::arrange()   masks plyr::arrange()
    ## x purrr::compact()   masks plyr::compact()
    ## x dplyr::count()     masks plyr::count()
    ## x dplyr::failwith()  masks plyr::failwith()
    ## x dplyr::filter()    masks stats::filter()
    ## x dplyr::id()        masks plyr::id()
    ## x dplyr::lag()       masks stats::lag()
    ## x dplyr::mutate()    masks plyr::mutate()
    ## x dplyr::rename()    masks plyr::rename()
    ## x dplyr::summarise() masks plyr::summarise()
    ## x dplyr::summarize() masks plyr::summarize()

``` r
library(gplots)
```

    ## 
    ## Attaching package: 'gplots'

    ## The following object is masked from 'package:stats':
    ## 
    ##     lowess

``` r
library(venn)
```

    ## 
    ## Attaching package: 'venn'

    ## The following object is masked from 'package:gplots':
    ## 
    ##     venn

``` r
# read ms DEqMS output file
# res_S4 <- readxl::read_excel("./out_r/020_Cord_blood_proteomics_CD34_Supplementary_tables_2020_06_15.xlsx" ,
#                                 sheet = "S4 MS results camera reactome",
#                                 col_names = TRUE,
# )
res_S4 <- readxl::read_excel("./out_r/Reactome/reactome_results.xlsx",
                             sheet = "ms_Camera",
                                col_names = TRUE,
                            
)


# Get abundances
loaded_abundancies <- load("./RData/003_ab_list.RData")
loaded_abundancies 
```

    ## [1] "ab_list"

``` r
# Get name of first pathways
res_S4$Pathway.name[1] 
```

    ## [1] "Plasma lipoprotein remodeling"

``` r
# Get gene symbols of first pathway
res_S4$Gene.symbols[1] 
```

    ## [1] "HPALP2,apo(a),APOA1;AD2,APO-E,ApoE4,LDLCQ5,LPG,APOE;Apo-AII,ApoA-II,apoAII,APOA2;APO-CII,APOC-II,APOC2;APOCIII,APOC3;HSA,PRO0883,PRO0903,PRO1341,ALB;FCHL2,FLDB,LDLCQ4,apoB-100,apoB-48,APOB;LCAT;NA;NA;BPIFF,HDLCQ10,CETP;BPIFE,HDLCQ9,PLTP"

``` r
# Get uniprotIDs of first pathway
res_S4$Submitted.entities.found[1] %>% 
  str_replace_all(";", ",")
```

    ## [1] "P04114,P55058,P04180,P02649,P06727,P08519,P02647,P02768,P02656,P02655,P02652,P11597"

``` r
# S2 MS data and res, per prot


keratinocytes <- c("Q14126","O43790","Q9Y5Y6","P01040","Q13835","P78386","P15924","Q6KB66","Q9C075","P19013","P22735","Q7Z794","Q02413","Q6A163","P07476","P22532","P35527","Q14532","P14923","P35321","Q08554","Q02487","Q86SJ6","Q92764","P78385","Q9NSB2","Q9NSB4","P05787","Q9Y446","P04259","P08779","Q15517","P23490","P13646","O76015","P13647","Q3SY84","O76013","P13645","O76011","P20930","Q15323","P04264","Q8N1N4","P31944","P35908","Q04695","Q14CN4","P02538","P07384","P02533")

dat.ms.fc <- ab_list[["dat_meta"]][["ms_meta"]] %>%
  select(-c("UniprotID","Gene.synonyms","Numb_NA_High","Numb_NA_Low","No.peptides","Unique.Peptides","MW.kDa","PSMs")) %>% 
  mutate_at(vars(-one_of(c("Gene.symbol"))),funs(as.numeric(.))) %>% 
  pivot_longer(cols=-1,
               names_to = c(".value","samp"),
               names_sep = "_") %>% 
  pivot_longer(cols = c(High,Low), names_to = "CD34", values_to = "Abundance") %>% 
  group_by(Gene.symbol,CD34) %>% 
  summarise(sum= sum(Abundance)) %>% 
  pivot_wider(names_from=CD34, values_from=sum) %>% 
  mutate(Foldchange=High-Low)
```

    ## Warning: `funs()` was deprecated in dplyr 0.8.0.
    ## Please use a list of either functions or lambdas: 
    ## 
    ##   # Simple named list: 
    ##   list(mean = mean, median = median)
    ## 
    ##   # Auto named with `tibble::lst()`: 
    ##   tibble::lst(mean, median)
    ## 
    ##   # Using lambdas
    ##   list(~ mean(., trim = .2), ~ median(., na.rm = TRUE))
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.

    ## `summarise()` has grouped output by 'Gene.symbol'. You can override using the
    ## `.groups` argument.

``` r
dat.ms.fc.unip <- ab_list[["dat_meta"]][["ms_meta"]] %>%
  select(-c("Gene.symbol","Gene.synonyms","Numb_NA_High","Numb_NA_Low","No.peptides","Unique.Peptides","MW.kDa","PSMs")) %>% 
  mutate_at(vars(-one_of(c("UniprotID"))),funs(as.numeric(.))) %>% 
  pivot_longer(cols=-1,
               names_to = c(".value","samp"),
               names_sep = "_") %>% 
  pivot_longer(cols = c(High,Low), names_to = "CD34", values_to = "Abundance") %>% 
  group_by(UniprotID,CD34) %>% 
  summarise(sum= sum(Abundance)) %>% 
  pivot_wider(names_from=CD34, values_from=sum) %>% 
  mutate(Foldchange=High-Low)
```

    ## `summarise()` has grouped output by 'UniprotID'. You can override using the
    ## `.groups` argument.

``` r
# All proteins of Plasma lipoprotein remodeling seems to have higher concentration in High CD34+
# Albumin, which has been depleeted, is one of the proteins. It shoudl problably be excluded.
dat.ms.fc %>% 
  filter(Gene.symbol%in%c("APOA1","APOE","APOA2","APOC2","ALB","APOB"))
```

    ## # A tibble: 6 × 4
    ## # Groups:   Gene.symbol [6]
    ##   Gene.symbol  High    Low Foldchange
    ##   <chr>       <dbl>  <dbl>      <dbl>
    ## 1 ALB         2.24  -0.544      2.78 
    ## 2 APOA1       1.04  -4.58       5.62 
    ## 3 APOA2       0.458 -5.67       6.13 
    ## 4 APOB        1.09   0.135      0.958
    ## 5 APOC2       7.50  -3.77      11.3  
    ## 6 APOE        3.31  -3.14       6.45

``` r
## UNIP
# Get uniprotIDs of first pathway
protsOfPathW <- res_S4$Submitted.entities.found[1] %>% 
  str_replace_all(";", ",") %>% 
  str_split(",") %>% 
  unlist
# Get FC
dat.ms.fc.unip %>% 
  filter(UniprotID %in% protsOfPathW)
```

    ## # A tibble: 12 × 4
    ## # Groups:   UniprotID [12]
    ##    UniprotID   High    Low Foldchange
    ##    <chr>      <dbl>  <dbl>      <dbl>
    ##  1 P02647     1.04  -4.58       5.62 
    ##  2 P02649     3.31  -3.14       6.45 
    ##  3 P02652     0.458 -5.67       6.13 
    ##  4 P02655     7.50  -3.77      11.3  
    ##  5 P02656     5.07  -3.76       8.83 
    ##  6 P02768     2.24  -0.544      2.78 
    ##  7 P04114     1.09   0.135      0.958
    ##  8 P04180     2.49  -1.35       3.84 
    ##  9 P06727    -0.885 -4.02       3.13 
    ## 10 P08519     6.35  NA         NA    
    ## 11 P11597     4.44  -1.03       5.47 
    ## 12 P55058     3.59  NA         NA

``` r
# Get uniprotIDs of second pathway
protsOfPathW <- res_S4$Submitted.entities.found[2] %>% 
  str_replace_all(";", ",") %>% 
  str_split(",") %>% 
  unlist
# Get FC
dat.ms.fc.unip %>% 
  filter(UniprotID %in% protsOfPathW)
```

    ## # A tibble: 64 × 4
    ## # Groups:   UniprotID [64]
    ##    UniprotID    High     Low Foldchange
    ##    <chr>       <dbl>   <dbl>      <dbl>
    ##  1 O00391     1.59   -0.789        2.38
    ##  2 O75083     1.36   -4.92         6.28
    ##  3 O94919    -1.24    0.0535      -1.29
    ##  4 P00441     0.488  -1.93         2.42
    ##  5 P00451     2.38   NA           NA   
    ##  6 P00488    -0.0595 -2.20         2.14
    ##  7 P00746    -3.37   -1.14        -2.24
    ##  8 P00747     0.391  -2.66         3.05
    ##  9 P01009     2.33   -3.97         6.29
    ## 10 P01011     2.90   -1.60         4.49
    ## # … with 54 more rows

``` r
## GS

# Get GS of first pathway
# pw1
protsOfPathW <- res_S4$Gene.symbols[1] %>% 
  str_replace_all(";", ",") %>% 
  str_split(",") %>% 
  unlist
# Get FC
dat.ms.fc %>% 
  filter(Gene.symbol %in% protsOfPathW)
```

    ## # A tibble: 10 × 4
    ## # Groups:   Gene.symbol [10]
    ##    Gene.symbol  High    Low Foldchange
    ##    <chr>       <dbl>  <dbl>      <dbl>
    ##  1 ALB         2.24  -0.544      2.78 
    ##  2 APOA1       1.04  -4.58       5.62 
    ##  3 APOA2       0.458 -5.67       6.13 
    ##  4 APOB        1.09   0.135      0.958
    ##  5 APOC2       7.50  -3.77      11.3  
    ##  6 APOC3       5.07  -3.76       8.83 
    ##  7 APOE        3.31  -3.14       6.45 
    ##  8 CETP        4.44  -1.03       5.47 
    ##  9 LCAT        2.49  -1.35       3.84 
    ## 10 PLTP        3.59  NA         NA

``` r
# pw2
protsOfPathW <- res_S4$Gene.symbols[2] %>% 
  str_replace_all(";", ",") %>% 
  str_split(",") %>% 
  unlist
# Get FC
dat.ms.fc %>% 
  filter(Gene.symbol %in% protsOfPathW)
```

    ## # A tibble: 64 × 4
    ## # Groups:   Gene.symbol [64]
    ##    Gene.symbol   High    Low Foldchange
    ##    <chr>        <dbl>  <dbl>      <dbl>
    ##  1 A1BG         0.560 -2.41       2.97 
    ##  2 A2M          1.62  -0.479      2.10 
    ##  3 ACTN1        0.531 NA         NA    
    ##  4 AHSG         1.05  -3.77       4.82 
    ##  5 ALB          2.24  -0.544      2.78 
    ##  6 ALDOA        0.715 -0.462      1.18 
    ##  7 APOA1        1.04  -4.58       5.62 
    ##  8 APOH        -0.579 -1.27       0.688
    ##  9 CAP1         1.75  -5.83       7.58 
    ## 10 CD109       -1.81  -0.819     -0.988
    ## # … with 54 more rows

``` r
#pw3
protsOfPathW <- res_S4$Gene.symbols[3] %>% 
  str_replace_all(";", ",") %>% 
  str_split(",") %>% 
  unlist
# Get FC
dat.ms.fc %>% 
  filter(Gene.symbol %in% protsOfPathW)
```

    ## # A tibble: 64 × 4
    ## # Groups:   Gene.symbol [64]
    ##    Gene.symbol   High    Low Foldchange
    ##    <chr>        <dbl>  <dbl>      <dbl>
    ##  1 A1BG         0.560 -2.41       2.97 
    ##  2 A2M          1.62  -0.479      2.10 
    ##  3 ACTN1        0.531 NA         NA    
    ##  4 AHSG         1.05  -3.77       4.82 
    ##  5 ALB          2.24  -0.544      2.78 
    ##  6 ALDOA        0.715 -0.462      1.18 
    ##  7 APOA1        1.04  -4.58       5.62 
    ##  8 APOH        -0.579 -1.27       0.688
    ##  9 CAP1         1.75  -5.83       7.58 
    ## 10 CD109       -1.81  -0.819     -0.988
    ## # … with 54 more rows

``` r
#pw4
protsOfPathW <- res_S4$Gene.symbols[4] %>% 
  str_replace_all(";", ",") %>% 
  str_split(",") %>% 
  unlist
# Get FC
dat.ms.fc %>% 
  filter(Gene.symbol %in% protsOfPathW)
```

    ## # A tibble: 12 × 4
    ## # Groups:   Gene.symbol [12]
    ##    Gene.symbol  High   Low Foldchange
    ##    <chr>       <dbl> <dbl>      <dbl>
    ##  1 ACTB        0.749 -2.47       3.22
    ##  2 FGA         1.87  -2.46       4.33
    ##  3 FGB         1.90  -2.84       4.74
    ##  4 FGG         2.40  -3.76       6.16
    ##  5 FN1         0.296 -1.92       2.21
    ##  6 ITGA2B      2.08  -2.94       5.01
    ##  7 ITGB3       1.53  -1.97       3.51
    ##  8 RAP1B       2.21  -4.40       6.61
    ##  9 TLN1        3.20  -3.39       6.60
    ## 10 VCL         1.17  -3.76       4.93
    ## 11 VWF         5.41  -4.35       9.76
    ## 12 YWHAB       0.902 -1.33       2.24

``` r
#pw5
protsOfPathW <- res_S4$Gene.symbols[5] %>% 
  str_replace_all(";", ",") %>% 
  str_split(",") %>% 
  unlist
# Get FC
dat.ms.fc %>% 
  filter(Gene.symbol %in% protsOfPathW)
```

    ## # A tibble: 13 × 4
    ## # Groups:   Gene.symbol [13]
    ##    Gene.symbol  High   Low Foldchange
    ##    <chr>       <dbl> <dbl>      <dbl>
    ##  1 ACTB        0.749 -2.47       3.22
    ##  2 FGA         1.87  -2.46       4.33
    ##  3 FGB         1.90  -2.84       4.74
    ##  4 FGG         2.40  -3.76       6.16
    ##  5 FN1         0.296 -1.92       2.21
    ##  6 ITGA2B      2.08  -2.94       5.01
    ##  7 ITGB3       1.53  -1.97       3.51
    ##  8 PEBP1       0.852 -1.68       2.53
    ##  9 RAP1B       2.21  -4.40       6.61
    ## 10 TLN1        3.20  -3.39       6.60
    ## 11 VCL         1.17  -3.76       4.93
    ## 12 VWF         5.41  -4.35       9.76
    ## 13 YWHAB       0.902 -1.33       2.24

``` r
#pw6
protsOfPathW <- res_S4$Gene.symbols[6] %>% 
  str_replace_all(";", ",") %>% 
  str_split(",") %>% 
  unlist
# Get FC
dat.ms.fc %>% 
  filter(Gene.symbol %in% protsOfPathW)
```

    ## # A tibble: 13 × 4
    ## # Groups:   Gene.symbol [13]
    ##    Gene.symbol  High   Low Foldchange
    ##    <chr>       <dbl> <dbl>      <dbl>
    ##  1 ACTB        0.749 -2.47       3.22
    ##  2 FGA         1.87  -2.46       4.33
    ##  3 FGB         1.90  -2.84       4.74
    ##  4 FGG         2.40  -3.76       6.16
    ##  5 FN1         0.296 -1.92       2.21
    ##  6 ITGA2B      2.08  -2.94       5.01
    ##  7 ITGB3       1.53  -1.97       3.51
    ##  8 PEBP1       0.852 -1.68       2.53
    ##  9 RAP1B       2.21  -4.40       6.61
    ## 10 TLN1        3.20  -3.39       6.60
    ## 11 VCL         1.17  -3.76       4.93
    ## 12 VWF         5.41  -4.35       9.76
    ## 13 YWHAB       0.902 -1.33       2.24

``` r
#pw7
protsOfPathW <- res_S4$Gene.symbols[7] %>% 
  str_replace_all(";", ",") %>% 
  str_split(",") %>% 
  unlist
# Get FC
dat.ms.fc %>% 
  filter(Gene.symbol %in% protsOfPathW)
```

    ## # A tibble: 8 × 4
    ## # Groups:   Gene.symbol [8]
    ##   Gene.symbol  High    Low Foldchange
    ##   <chr>       <dbl>  <dbl>      <dbl>
    ## 1 ALB          2.24 -0.544       2.78
    ## 2 APOA1        1.04 -4.58        5.62
    ## 3 APOC2        7.50 -3.77       11.3 
    ## 4 APOC3        5.07 -3.76        8.83
    ## 5 APOE         3.31 -3.14        6.45
    ## 6 CETP         4.44 -1.03        5.47
    ## 7 LCAT         2.49 -1.35        3.84
    ## 8 PLTP         3.59 NA          NA

``` r
#pw8
protsOfPathW <- res_S4$Gene.symbols[8] %>% 
  str_replace_all(";", ",") %>% 
  str_split(",") %>% 
  unlist
# Get FC
dat.ms.fc %>% 
  filter(Gene.symbol %in% protsOfPathW)
```

    ## # A tibble: 9 × 4
    ## # Groups:   Gene.symbol [9]
    ##   Gene.symbol  High   Low Foldchange
    ##   <chr>       <dbl> <dbl>      <dbl>
    ## 1 FGA         1.87  -2.46       4.33
    ## 2 FGB         1.90  -2.84       4.74
    ## 3 FGG         2.40  -3.76       6.16
    ## 4 FN1         0.296 -1.92       2.21
    ## 5 ITGA2B      2.08  -2.94       5.01
    ## 6 ITGB3       1.53  -1.97       3.51
    ## 7 RAP1B       2.21  -4.40       6.61
    ## 8 TLN1        3.20  -3.39       6.60
    ## 9 VWF         5.41  -4.35       9.76

``` r
#pw9
protsOfPathW <- res_S4$Gene.symbols[9] %>% 
  str_replace_all(";", ",") %>% 
  str_split(",") %>% 
  unlist
# Get FC
dat.ms.fc %>% 
  filter(Gene.symbol %in% protsOfPathW)
```

    ## # A tibble: 6 × 4
    ## # Groups:   Gene.symbol [6]
    ##   Gene.symbol   High   Low Foldchange
    ##   <chr>        <dbl> <dbl>      <dbl>
    ## 1 CAPN1       -0.733 0.123     -0.856
    ## 2 DSC2        -0.243 0.381     -0.624
    ## 3 KRT1        -6.28  2.01      -8.29 
    ## 4 KRT16       -4.05  2.00      -6.05 
    ## 5 KRT6A       -3.94  4.55      -8.49 
    ## 6 SPRR2D      -6.62  6.50     -13.1

``` r
#pw10
protsOfPathW <- res_S4$Gene.symbols[10] %>% 
  str_replace_all(";", ",") %>% 
  str_split(",") %>% 
  unlist
# Get FC
dat.ms.fc %>% 
  filter(Gene.symbol %in% protsOfPathW)
```

    ## # A tibble: 6 × 4
    ## # Groups:   Gene.symbol [6]
    ##   Gene.symbol   High   Low Foldchange
    ##   <chr>        <dbl> <dbl>      <dbl>
    ## 1 CAPN1       -0.733 0.123     -0.856
    ## 2 DSC2        -0.243 0.381     -0.624
    ## 3 KRT1        -6.28  2.01      -8.29 
    ## 4 KRT16       -4.05  2.00      -6.05 
    ## 5 KRT6A       -3.94  4.55      -8.49 
    ## 6 SPRR2D      -6.62  6.50     -13.1

``` r
## Make table with loop construct

res_S4_fdr_sign <- res_S4 %>% 
  filter(Entities.FDR<0.05)

fc_list <- list()
for (pw in 1:nrow(res_S4_fdr_sign)){
  #print("hej")
  protsOfPathW <- res_S4$Gene.symbols[pw] %>% 
    str_replace_all(";", ",") %>% 
    str_split(",") %>% 
    unlist
  # Get FC
  fc_list[[res_S4$Pathway.name[pw]]] <- dat.ms.fc %>% 
    filter(Gene.symbol %in% protsOfPathW)
}

fc_tab <- fc_list %>% 
  bind_rows(.id = "pathway") 


fc_tab_group <-  fc_tab %>% 
  ungroup() %>% 
  group_by(pathway) %>% 
  summarise(mean_fc=mean(Foldchange,na.rm =T)) %>% 
  mutate(Conc_change = if_else(mean_fc <= 0, "Lower in High CD34", "Higher in High CD34"))

fc_tab_group %>% 
  writexl::write_xlsx(path = "./out_r/014_table_with_fc_per_pw_grouped.xlsx")

fc_tab %>% 
  writexl::write_xlsx(path = "./out_r/014_table_with_fc_per_pw.xlsx")



dat.ms.fc %>% 
  slice(grep("KRT", Gene.symbol, invert = F))
```

    ## # A tibble: 5 × 4
    ## # Groups:   Gene.symbol [5]
    ##   Gene.symbol  High   Low Foldchange
    ##   <chr>       <dbl> <dbl>      <dbl>
    ## 1 KRT1        -6.28  2.01      -8.29
    ## 2 KRT16       -4.05  2.00      -6.05
    ## 3 KRT34       -1.95  8.87     -10.8 
    ## 4 KRT4        -1.63  3.38      -5.02
    ## 5 KRT6A       -3.94  4.55      -8.49

``` r
# %>% 
#   mutate_at(vars(c("High","Low")),as.numeric)



##################

prot_list <- list()
for (pw in 1:nrow(res_S4_fdr_sign)){
  #print("hej")
  protsOfPathW <- res_S4$Submitted.entities.found[pw] %>% 
    str_replace_all(";", ",") %>% 
    str_split(",") %>% 
    unlist
  # Get FC
  prot_list[[res_S4$Pathway.name[pw]]] <- protsOfPathW
}

pw_per_prot.0.0 <- prot_list %>% 
  unlist(use.names = T) %>% 
  list() %>% 
  data.frame() 

names(pw_per_prot.0.0) <- c("UniprotID")


pw_per_prot.0.1 <- pw_per_prot.0.0 %>% 
  bind_cols(row.names(pw_per_prot.0.0)) %>% 
  as_tibble()
```

    ## New names:
    ## * NA -> ...2

``` r
names(pw_per_prot.0.1) <- c("UniprotID","PW")

pw_no_digits <- pw_per_prot.0.1$PW %>% 
  str_replace_all("[:digit:]", "")

pw_per_prot.0.2 <- bind_cols(pw_per_prot.0.1,pw_no_digits) %>% 
  select(-PW)
```

    ## New names:
    ## * NA -> ...3

``` r
names(pw_per_prot.0.2) <- c("UniprotID","PW")

pw_per_prot.0.2
```

    ## # A tibble: 322 × 2
    ##    UniprotID PW                           
    ##    <chr>     <chr>                        
    ##  1 P04114    Plasma lipoprotein remodeling
    ##  2 P55058    Plasma lipoprotein remodeling
    ##  3 P04180    Plasma lipoprotein remodeling
    ##  4 P02649    Plasma lipoprotein remodeling
    ##  5 P06727    Plasma lipoprotein remodeling
    ##  6 P08519    Plasma lipoprotein remodeling
    ##  7 P02647    Plasma lipoprotein remodeling
    ##  8 P02768    Plasma lipoprotein remodeling
    ##  9 P02656    Plasma lipoprotein remodeling
    ## 10 P02655    Plasma lipoprotein remodeling
    ## # … with 312 more rows

``` r
pw_per_prot.0.3 <- pw_per_prot.0.2 %>% 
  group_by(UniprotID) %>%
  summarise(PW = toString(PW))


prots_ranked_by_pw_occ <- prot_list %>% 
  unlist() %>% 
  table() %>% 
  as_tibble() %>% 
  arrange(desc(n))

names(prots_ranked_by_pw_occ) <- c("UniprotID","numbPW")

prots_ranked_by_pw_occ
```

    ## # A tibble: 110 × 2
    ##    UniprotID numbPW
    ##    <chr>      <int>
    ##  1 P02671        12
    ##  2 P02675        12
    ##  3 P02679        12
    ##  4 P02751        12
    ##  5 P04275        12
    ##  6 P05106        12
    ##  7 P08514        12
    ##  8 Q9Y490        12
    ##  9 P61224        10
    ## 10 P18206         9
    ## # … with 100 more rows

``` r
pw_per_prot.0.4 <- prots_ranked_by_pw_occ %>% 
  full_join(pw_per_prot.0.3)
```

    ## Joining, by = "UniprotID"

``` r
my_loaded_objs <- load("./RData/12_formatted_out.RData")
my_loaded_objs
```

    ## [1] "ms_filt_norm_test"  "pea_filt_norm_test"

``` r
GS.0.0 <- ms_filt_norm_test %>% 
  select(UniprotID,Gene.symbol)

pw_per_prot.0.5 <- GS.0.0 %>% 
  right_join(pw_per_prot.0.4) %>% 
  distinct() %>% 
  arrange(desc(numbPW))
```

    ## Joining, by = "UniprotID"

``` r
pw_per_prot.0.5  %>%
  writexl::write_xlsx(path = "./out_r/014_pw_per_prot.xlsx")
 

prot_list %>% 
  unlist() %>% 
  unique()
```

    ##   [1] "P04114" "P55058" "P04180" "P02649" "P06727" "P08519" "P02647" "P02768"
    ##   [9] "P02656" "P02655" "P02652" "P11597" "P08514" "P01042" "P12814" "P12259"
    ##  [17] "P13473" "P00747" "O00391" "P00746" "P02765" "P08567" "P05452" "Q99969"
    ##  [25] "P07996" "P29622" "Q13201" "Q96C24" "P16109" "P18206" "P11021" "P02679"
    ##  [33] "Q92520" "P02675" "O94919" "Q01518" "P00451" "P02751" "P08697" "P01023"
    ##  [41] "P21333" "P02671" "P05106" "Q9NXH8" "P04217" "Q08380" "P62937" "Q6YHK3"
    ##  [49] "Q14624" "P02749" "P16284" "P23528" "P00488" "Q16610" "P02787" "P01137"
    ##  [57] "P05155" "Q9Y490" "P00441" "P01011" "Q86UX7" "Q13103" "P07737" "P37802"
    ##  [65] "Q06033" "O75083" "P49908" "P01009" "P04075" "P05121" "P02775" "P04275"
    ##  [73] "P10909" "P07225" "P31946" "P60709" "P61224" "P30086" "P08779" "Q02487"
    ##  [81] "P02538" "P07384" "P19013" "P04264" "O76011" "P22532" "P61619" "P13639"
    ##  [89] "P46782" "P36578" "P13798" "P84098" "P15880" "P68104" "Q02543" "P12081"
    ##  [97] "P63220" "P62987" "P15170" "P61247" "P62851" "P46779" "P35232" "P01130"
    ## [105] "P13497" "Q8NBP7" "Q86X29" "P02654" "P07237" "P62993"

``` r
# Response to elevated platelet cytosolic Ca2+ (Same as Platelet degranulation)
# Signaling by RAF1 mutants
# Signaling by high-kinase activity BRAF mutants
# Signaling downstream of RAS mutants
# Paradoxical activation of RAF signaling by kinase inactive BRAF
# Signaling by moderate kinase activity BRAF mutants
# Integrin signaling (same as MAP2K and MAPK activation except 1)

# Keratinization (has its own 8 prots), 
# Formation of the cornified envelope (same as Keratinization)
# Translation (has its own 16 prots)

# GRB2:SOS provides linkage to MAPK signaling for Integrins 

selected_pws1 <- c("Platelet degranulation",
                   "MAP2K and MAPK activation",
                   "Plasma lipoprotein remodeling",
                   "Signaling by RAS mutants",
                   "p130Cas linkage to MAPK signaling for integrins")
col1 <- 1:5

selected_pws2 <- c("Platelet degranulation",
                   "Response to elevated platelet cytosolic Ca2+",
                   "Signaling by RAF1 mutants",
                   "Signaling by high-kinase activity BRAF mutants",
                   "MAP2K and MAPK activation",
                   "p130Cas linkage to MAPK signaling for integrins",
                   "Signaling by RAS mutants")

selected_pws3 <- c("Paradoxical activation of RAF signaling by kinase inactive BRAF",
                   "Signaling by moderate kinase activity BRAF mutants", # Signaling downstream of RAS mutants (only the nine in center)
                   "Signaling by RAF1 mutants",
                   "Signaling by high-kinase activity BRAF mutants",
                   "MAP2K and MAPK activation",
                   "GRB2:SOS provides linkage to MAPK signaling for Integrins", # p130Cas linkage to MAPK signaling for integrins (only the ninein center)
                   "Platelet degranulation"# Signaling by RAS mutants (only the nine in center)
                   )

selected_pws4 <- c("Signaling by RAF1 mutants",
                   "Signaling downstream of RAS mutants",
                   "Signaling by moderate kinase activity BRAF mutants",
                   "p130Cas linkage to MAPK signaling for integrins",
                   "GRB2:SOS provides linkage to MAPK signaling for Integrins")

selected_pws5 <- c("Paradoxical activation of RAF signaling by kinase inactive BRAF",
                   "Signaling downstream of RAS mutants",
                   "Signaling by RAF1 mutants",
                   "Signaling by high-kinase activity BRAF mutants",
                   "MAP2K and MAPK activation",
                   "p130Cas linkage to MAPK signaling for integrins",
                   "Signaling by RAS mutants"
)

# "Signaling downstream of RAS mutants"
col2 <- 1:7

prot_list_5.1 <- prot_list[selected_pws1]
prot_list_5.2 <- prot_list[selected_pws2]
prot_list_5.3 <- prot_list[selected_pws3]
prot_list_5.4 <- prot_list[selected_pws4]

pdf(file = "./out_r/14_venn_pw_overlap_2.pdf",   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 15) # The height of the plot in inches
venn::venn(prot_list_5.2, zcolor=col2,sncs=1.2,ilcs=1.2) # , size = 1.5 , ggplot = TRUE
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
pdf(file = "./out_r/14_venn_pw_overlap_3.pdf",   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 15) # The height of the plot in inches
venn::venn(prot_list_5.3, zcolor=col2,sncs=1.2,ilcs=1.2) # , size = 1.5 , ggplot = TRUE
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
venn(prot_list_5.4, zcolor=1:5,sncs=1.2,ilcs=1.2) 
```

![](/Users/xrydbh/Personal/Projects/git_cloned/Neored_fresh/Neored/markdown/014_doublecheck_foldchange_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->
