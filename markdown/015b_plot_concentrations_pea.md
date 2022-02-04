015b\_plot\_concentrations\_pea.R
================
xrydbh
2022-02-04

``` r
library(tidyverse)
library(RColorBrewer)
library(readxl)
library(writexl)

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).



loaded_cyt_objs <- load("./RData/15_annoDb.cytokine.chemokine.RData")

##############################
# Format data for rank plot
########################


# read pea file
pea.raw.0.0 <- readxl::read_excel("./data/001_pea_raw.xlsx") 

missing_values_allowed_for_prot <- 0 #3 for with_at_least_5_vals (_in_high_or_low) or 
filter_missing_val_both_groups_or_either_or <- "either_or" #"either_or" # or "both"


# make longer

NPX.0.3 <- pea.raw.0.0 %>% 
  pivot_longer(cols=contains("_"),names_to="SC_Donor",values_to = "NPX") %>% 
  mutate(NPX=as.numeric(NPX)) %>% 
  mutate(LOD=as.numeric(LOD)) 

######################################
# Count numb of NAs in High and Low
# group by proteins
#################################3####
# Separate column SC_Donor into "High low Group" and "Unit 1-8"
NPX.0.3.5 <- NPX.0.3 %>%
  separate(SC_Donor,
           into = c("Group","Unit"),
           sep = "_")

# Get number of NPX below LOD (missing values) per gene and group
numb_NPX_vals_below_LOD.0.0 <- NPX.0.3.5 %>% 
  group_by(Gene.symbol,Group) %>% 
  summarise(NPX_below_LOD = sum(NPX<LOD)) %>% 
  pivot_wider(names_from = Group,values_from =NPX_below_LOD) %>% 
  dplyr::rename(Numb_NA_High=High,Numb_NA_Low=Low)
```

    ## `summarise()` has grouped output by 'Gene.symbol'. You can override using the
    ## `.groups` argument.

``` r
NPX.0.3.1 <- full_join(NPX.0.3,numb_NPX_vals_below_LOD.0.0 ,by="Gene.symbol")



# Replace UniprotID NA for Gene.symbol NT-proBNP and OlinkID  OID00131 and OID01214 by  P16860
# Also calculate average between OID00131 and OID01214

#what OlinkID have Corresponding UniprotID as NA
#NPX.0.0 %>% filter(is.na(UniProt )) %>% pluck("OlinkID") %>% unique()


# Get average for the two NT-proBNP Olink variants
NT_proBNP <- NPX.0.3.1  %>% 
  filter(Gene.symbol=="NT-proBNP") %>% 
  group_by(SC_Donor) %>% 
  mutate(NPX=(mean(NPX))) %>% # Calculate mean
  ungroup() %>% 
  filter(OlinkID=="OID00131") # keep only on of the two laternative OlinkIDs

# Replace duplicates with average
NPX.0.4 <- NPX.0.3.1  %>% 
  filter(Gene.symbol!="NT-proBNP") %>%
  bind_rows(NT_proBNP)

# Replace UniprotID NA for NT-proBNP with "P16860"
NPX.0.4$UniprotID[is.na(NPX.0.4$UniprotID)] <- "P16860"

# Remove duplicate UniprotIDs "Q8NEV9,Q14213" (IL27a/IL27b) with only one ID
NPX.0.4$UniprotID[NPX.0.4$UniprotID=="Q8NEV9,Q14213"] <- "Q8NEV9"
# Remove duplicate UniprotIDs "Q11128,P21217" (FUT3/FUT5) with only one ID
NPX.0.4$UniprotID[NPX.0.4$UniprotID=="Q11128,P21217"] <- "Q11128"


#NPX.0.4 %>% pivot_wider(names_from=SC_Donor,values_from = NPX)
####################

# summarize and arrange by NPX
pea.for.rankplot.0.0 <- NPX.0.4 %>% 
  group_by(UniprotID) %>% # group by gene symbol
  summarize(sum_npx=sum(NPX,na.rm = T)) %>% # summarize into sum per genes across samples
  filter(sum_npx!=0) %>% # remove sums that are zero
  arrange(desc(sum_npx)) 

pea.for.rankplot <- pea.for.rankplot.0.0 %>% # sort
  mutate(rank=seq(1:nrow(pea.for.rankplot.0.0))) # crate rank variable

pea.for.rankplot.gs <- pea.for.rankplot %>% inner_join(NPX.0.4,by=c("UniprotID"="UniprotID")) %>% 
  dplyr::select(Gene.symbol,sum_npx,rank) %>% 
  distinct()

###########################
# check overlap with Cytokines and chemokines
#######################

# Get cytokines as gs in ms 
cyts_in_pea.gs3 <- pea.for.rankplot.gs %>% 
  inner_join(annoDb.cytokine.chemokine,by=c("Gene.symbol"="Gene.symbol")) %>% 
  dplyr::select(-UniprotID) %>% 
  distinct() %>% 
  mutate(GO_term=replace(GO_term, GO_term=="chemokine receptor antagonist activity", "chemokine receptor binding"))
cyts_in_pea.gs3$GO_term %>% as.factor() %>% levels()
```

    ## [1] "chemokine activity"                                          
    ## [2] "chemokine receptor binding"                                  
    ## [3] "cytokine activity"                                           
    ## [4] "cytokine receptor activity"                                  
    ## [5] "cytokine receptor binding"                                   
    ## [6] "hematopoietic progenitor cell differentiation"               
    ## [7] "hematopoietic stem cell migration"                           
    ## [8] "hematopoietic stem cell migration to bone marrow"            
    ## [9] "positive regulation of hematopoietic stem cell proliferation"

``` r
# get number of ge terms matching in ms data
goTerminPEACounts <-cyts_in_pea.gs3 %>% 
  distinct %>% 
  group_by(GO_term) %>% 
  tally()

# join
cytokine_summary_table.pea.0.0 <- goTermCounts %>% 
  full_join(goTerminPEACounts,by=c("GO_term"="GO_term"))

# replace NA with 0
cytokine_summary_table.pea.0.0$n.y <- cytokine_summary_table.pea.0.0$n.y %>% 
  tidyr::replace_na(0)

# get percenttge
cytokine_summary_table.pea.0.1 <- cytokine_summary_table.pea.0.0 %>% 
  dplyr::rename(nInAnnoDB=n.x,nInMS=n.y) %>% 
  mutate(perc_in=(nInMS/nInAnnoDB)*100) %>% 
  mutate(perc_out=((nInAnnoDB-nInMS)/nInAnnoDB)*100)

# listGO_term gene members not in MS
GO_term_gene_members_not_in_pea <- annoDb.cytokine.chemokine %>% 
  dplyr::filter(!(Gene.symbol%in%cyts_in_pea.gs3$Gene.symbol)) %>% 
  dplyr::select(-UniprotID) %>% 
  distinct() %>% 
  dplyr::group_by(GO_term) %>% 
  dplyr::summarise(Gene.symbol = paste(Gene.symbol, collapse = ","))

# list  GO_term MS gene members
GO_term_pea_gene_members <- cyts_in_pea.gs3 %>% 
  dplyr::group_by(GO_term) %>% 
  dplyr::summarise(Gene.symbol = paste(Gene.symbol, collapse = ","))

cytokine_summary_table.pea.0.2  <- cytokine_summary_table.pea.0.1  %>% 
  full_join(GO_term_pea_gene_members,by=c("GO_term"="GO_term")) %>% 
  dplyr::rename(In_pea=Gene.symbol) %>% 
  full_join(GO_term_gene_members_not_in_pea,by=c("GO_term"="GO_term")) %>% 
  dplyr::rename(Not_in_pea=Gene.symbol) 

#########################
# plot
##########################

# plot with some names
rankplot.pea <- pea.for.rankplot.gs %>% 
  ggplot(aes(x=rank, y=log2(sum_npx))) + 
  geom_point(color="grey",alpha=1) +
  #geom_point(data=cyts_in_ms.gs,color="red") +
  ggrepel::geom_text_repel(aes(x = rank, 
                               y = log2(sum_npx),
                               color=GO_term, 
                               label = Gene.symbol,
                               fontface = "bold"),
                           data=cyts_in_pea.gs3,
                           max.overlaps=20,
                           nudge_x=90
  ) +  
  #scale_color_manual(values = mycolors)+
  #scale_color_brewer(palette="Paired")
  # scale_color_brewer(values = "Paired") + #    scale_color_manual(values = rainbow(13)) +  #
  
  #geom_text(aes(label=Gene.symbol,color=as.factor(GO_term)),data=cyts_in_ms.gs3) + #, check_overlap = TRUE, nudge_x = 70
  geom_text(aes(label=Gene.symbol),data=head(pea.for.rankplot.gs,n=10L), nudge_x = -90,check_overlap = TRUE,alpha=0.5)+ #, nudge_x = -90, position=position_jitter(width=1,height=1),
  geom_text(aes(label=Gene.symbol),data=tail(pea.for.rankplot.gs,n=10L),check_overlap = TRUE, nudge_x = -90,alpha=0.5)+
  ggtitle(paste("Plot of norm_ab sum and rank of ",nrow(pea.for.rankplot)," PEA proteins",sep="") )
ggtitle(paste("Demarkating ",nrow(pea.for.rankplot.gs)," MS proteins",sep="") )
```

    ## $title
    ## [1] "Demarkating 459 MS proteins"
    ## 
    ## attr(,"class")
    ## [1] "labels"

``` r
rankplot.pea
```

    ## Warning in FUN(X[[i]], ...): NaNs produced

    ## Warning in FUN(X[[i]], ...): NaNs produced

    ## Warning in FUN(X[[i]], ...): NaNs produced

    ## Warning in FUN(X[[i]], ...): NaNs produced

    ## Warning: Removed 8 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_text_repel).

    ## Warning: Removed 8 rows containing missing values (geom_text).

![](/Users/xrydbh/Personal/Projects/git_cloned/Neored_fresh/Neored/markdown/015b_plot_concentrations_pea_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
pdf(file = "./out_r/15_b_rankplot_concentrations_pea.pdf",   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 15) # The height of the plot in inches
rankplot.pea
```

    ## Warning in FUN(X[[i]], ...): NaNs produced

    ## Warning in FUN(X[[i]], ...): NaNs produced

    ## Warning in FUN(X[[i]], ...): NaNs produced

    ## Warning in FUN(X[[i]], ...): NaNs produced

    ## Warning: Removed 8 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_text_repel).

    ## Warning: Removed 8 rows containing missing values (geom_text).

``` r
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
write_xlsx(cytokine_summary_table.pea.0.2,"./out_r/015b_cytokine_summary_table_pea.xlsx")
```
