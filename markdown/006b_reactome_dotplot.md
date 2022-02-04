006b\_reactome\_dotplot.R
================
xrydbh
2022-02-04

``` r
library(clusterProfiler)
```

    ## clusterProfiler v3.18.1  For help: https://guangchuangyu.github.io/software/clusterProfiler
    ## 
    ## If you use clusterProfiler in published research, please cite:
    ## Guangchuang Yu, Li-Gen Wang, Yanyan Han, Qing-Yu He. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology. 2012, 16(5):284-287.

    ## 
    ## Attaching package: 'clusterProfiler'

    ## The following object is masked from 'package:OrganismDbi':
    ## 
    ##     select

    ## The following object is masked from 'package:AnnotationDbi':
    ## 
    ##     select

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     slice

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     rename

    ## The following object is masked from 'package:purrr':
    ## 
    ##     simplify

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
#library(enrichplot)
#library(ggplot2)
library(tidyverse)

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).



# read long read protein centric
ms.0.0 <- readxl::read_excel("./out_r/Reactome/reactome_results.xlsx",sheet=1)
pea.0.0 <- readxl::read_excel("./out_r/Reactome/reactome_results.xlsx",sheet=2)


#pdf("./out_r/GO_BP/dotplot_ORA_ms_in_high_ora_go_bp.pdf",
#   width=10, 
# #    height=12)
# dotplot(ms_over_in_high_ora_go_bp, showCategory=30) + ggtitle("dotplot for ORA ms_over_in_high_ora_go_bp")
# #dev.off()
# 
# dotplot(ms.0.0, showCategory=30) + ggtitle("dotplot for ORA ms_over_in_high_ora_go_bp")
# 
# ms.0.0 %>% names()

ms.0.1 <- ms.0.0 %>% 
  filter(Entities.FDR<0.05) %>% 
  arrange(Entities.pValue) 

ms.0.2 <- ms.0.1 %>% 
  mutate(rank=nrow(ms.0.1):1) #%>% 
  # mutate(Pathway.name=as.factor(rank)) %>% 
  # mutate(Pathway.name=as.factor(Pathway.name))

ms.plot <- ms.0.2 %>% 
  mutate(Pathway.name = fct_reorder(Pathway.name, rank)) %>%
  ggplot(aes(x = Entities.ratio, y = Pathway.name, 
                        color = `Entities.FDR`, size = X.Entities.found)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("Entities.ratio") + 
  ggtitle("Significant pathways - MS Reactome Camera")

pdf("./out_r/GO_BP/dotplot_reactome_camera_ms.pdf",
  width=8,
   height=5)
ms.plot
dev.off()
```

    ## quartz_off_screen 
    ##                 2
