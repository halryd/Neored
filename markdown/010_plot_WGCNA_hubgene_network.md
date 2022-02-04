010\_plot\_WGCNA\_hubgene\_network.R
================
xrydbh
2022-02-04

``` r
library(qgraph)
library(igraph)
```

    ## 
    ## Attaching package: 'igraph'

    ## The following object is masked from 'package:clusterProfiler':
    ## 
    ##     simplify

    ## The following object is masked from 'package:DescTools':
    ## 
    ##     %c%

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

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

``` r
library(tidyverse)

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).



networkWeightCutOff <- 30

 for (plf in c("ms","pea")){
         input_path <- paste("./out_r/WGCNA/",plf,sep="")
         
         output_path.0.0 <- input_path 
         output_path <- paste(output_path.0.0,"/network_plots",sep="")   #"plots_wgcna_Gene.symbol/ms_Log2_norm_ab_0_na_both/softPower_12/network_plots" # coul also be set equla to plot_dir form 3_b
         dir.create(output_path)
         
         ## Plot networks based on topological overlap
         full_fns_TOM <- dir(input_path,"-TOM",full.names=T)
         full_extracts_and_matches <- full_fns_TOM %>% 
                 str_match("TOM(.*?).txt")
         modules_of_interest <- full_extracts_and_matches[,2]
         
         for (i in seq_along(modules_of_interest)){ #
                 print(modules_of_interest[i])
                 edglist.0.0 <- read.table(full_fns_TOM[i], 
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
                 pdf(paste(output_path,"/TOM_network_heatmap_",modules_of_interest[i],"_module.pdf",sep=""))#plot_dir #, 2000, 2000, res=300
                 qgraph(my_mat.0.0,
                        graph= "default", #cor
                        layout = "spring",
                        edge.labels = TRUE,
                        palette= "ggplot2",
                        #groups=sign_cyts_and_gfs.lst_0[[i]],
                        legend=F,
                        #bg='lightgrey',
                        title.cex=1.2,
                        label.cex=1,
                        #label.norm="OOOOO",
                        title=paste("Hubgene TOM network",modules_of_interest[i]),
                        label.scale.equal=T,
                        labels = colnames(my_mat.0.0))
                 dev.off()
                 
         }
         
         ## Plot networks based on adjacency
         full_fns_ADJ <- dir(input_path,"-ADJ",full.names=T)
         full_extracts_and_matches <- full_fns_ADJ %>% str_match("ADJ(.*?).txt")
         modules_of_interest <- full_extracts_and_matches[,2]
         for (i in seq_along(modules_of_interest)){
                 edglist.0.0 <- read.table(full_fns_ADJ[i], 
                                           sep = "" , 
                                           header = F,
                                           na.strings ="", 
                                           stringsAsFactors= F) %>% 
                         as.matrix()#, nrows = 100
                 edgelist.0.1 <- edglist.0.0[,c(1,2,5)]
                 colnames(edgelist.0.1) <- c("From","To","Weight")
                 edgelist.0.2 <- edgelist.0.1 %>% as_tibble() %>% 
                         arrange(desc(Weight)) %>% head(networkWeightCutOff )
                 # convert edgelist to ?? matrix
                 g=graph.data.frame(edgelist.0.2, directed = F)
                 my_mat.0.0 <- get.adjacency(g, sparse = FALSE, attr='Weight')
                 class(my_mat.0.0) <- "numeric"
                 png(paste(output_path,"/ADJ_network_heatmap_",modules_of_interest[i],"_module.png",sep=""), 2000, 2000, res=300)#plot_dir
                 qgraph(my_mat.0.0,
                        graph= "default", #cor
                        layout = "spring",
                        edge.labels = TRUE,
                        palette= "ggplot2",
                        #groups=sign_cyts_and_gfs.lst_0[[i]],
                        legend=F,
                        bg='lightgrey',
                        title.cex=1.2,
                        label.cex=1,
                        #label.norm="OOOOO",
                        title=paste("Hubgene adjacency network",modules_of_interest[i]),
                        label.scale.equal=T,
                        labels = colnames(my_mat.0.0))
                 dev.off()
                 
         }
         
         
         
 }
```

    ## Warning in dir.create(output_path): './out_r/WGCNA/ms/network_plots' already
    ## exists

    ## [1] "blue"

    ## Warning in qgraph(my_mat.0.0, graph = "default", layout = "spring", edge.labels
    ## = TRUE, : Non-finite weights are omitted

    ## [1] "brown"

    ## Warning in qgraph(my_mat.0.0, graph = "default", layout = "spring", edge.labels
    ## = TRUE, : Non-finite weights are omitted

    ## [1] "turquoise"

    ## Warning in qgraph(my_mat.0.0, graph = "default", layout = "spring", edge.labels
    ## = TRUE, : Non-finite weights are omitted

    ## Warning in qgraph(my_mat.0.0, graph = "default", layout = "spring", edge.labels
    ## = TRUE, : Non-finite weights are omitted

    ## Warning in qgraph(my_mat.0.0, graph = "default", layout = "spring", edge.labels
    ## = TRUE, : Non-finite weights are omitted

    ## Warning in qgraph(my_mat.0.0, graph = "default", layout = "spring", edge.labels
    ## = TRUE, : Non-finite weights are omitted

    ## Warning in dir.create(output_path): './out_r/WGCNA/pea/network_plots' already
    ## exists

    ## [1] "black"

    ## Warning in qgraph(my_mat.0.0, graph = "default", layout = "spring", edge.labels
    ## = TRUE, : Non-finite weights are omitted

    ## [1] "blue"

    ## Warning in qgraph(my_mat.0.0, graph = "default", layout = "spring", edge.labels
    ## = TRUE, : Non-finite weights are omitted

    ## [1] "brown"

    ## Warning in qgraph(my_mat.0.0, graph = "default", layout = "spring", edge.labels
    ## = TRUE, : Non-finite weights are omitted

    ## [1] "green"

    ## Warning in qgraph(my_mat.0.0, graph = "default", layout = "spring", edge.labels
    ## = TRUE, : Non-finite weights are omitted

    ## [1] "red"

    ## Warning in qgraph(my_mat.0.0, graph = "default", layout = "spring", edge.labels
    ## = TRUE, : Non-finite weights are omitted

    ## [1] "turquoise"

    ## Warning in qgraph(my_mat.0.0, graph = "default", layout = "spring", edge.labels
    ## = TRUE, : Non-finite weights are omitted

    ## [1] "yellow"

    ## Warning in qgraph(my_mat.0.0, graph = "default", layout = "spring", edge.labels
    ## = TRUE, : Non-finite weights are omitted

    ## Warning in qgraph(my_mat.0.0, graph = "default", layout = "spring", edge.labels
    ## = TRUE, : Non-finite weights are omitted

    ## Warning in qgraph(my_mat.0.0, graph = "default", layout = "spring", edge.labels
    ## = TRUE, : Non-finite weights are omitted

    ## Warning in qgraph(my_mat.0.0, graph = "default", layout = "spring", edge.labels
    ## = TRUE, : Non-finite weights are omitted

    ## Warning in qgraph(my_mat.0.0, graph = "default", layout = "spring", edge.labels
    ## = TRUE, : Non-finite weights are omitted

    ## Warning in qgraph(my_mat.0.0, graph = "default", layout = "spring", edge.labels
    ## = TRUE, : Non-finite weights are omitted

    ## Warning in qgraph(my_mat.0.0, graph = "default", layout = "spring", edge.labels
    ## = TRUE, : Non-finite weights are omitted

    ## Warning in qgraph(my_mat.0.0, graph = "default", layout = "spring", edge.labels
    ## = TRUE, : Non-finite weights are omitted
