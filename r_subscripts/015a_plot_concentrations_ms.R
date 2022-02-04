##############################  
# Get all cytokines with AnnotationDB
############################## 
library(AnnotationDbi)
library("Homo.sapiens")
library(GO.db)
library(tidyverse)
library(RColorBrewer)
library(readxl)
library(writexl)

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).




columns(GO.db)
head(keys(GO.db, keytype="DEFINITION"))
head(keys(GO.db, keytype="GOID"))
head(keys(GO.db, keytype="ONTOLOGY"))

# Get all terms
all_go_terms <- keys(GO.db, keytype="TERM") 

# Get my terms with regular expression
my_cytokine_go_terms <- grep("^cytokine activity", all_go_terms, value = TRUE)
my_cytokine_receptor_go_terms <- grep("cytokine receptor", all_go_terms, value = TRUE)

my_chemokine_go_terms <- grep("^chemokine activity", all_go_terms, value = TRUE)
my_chemokine_receptor_go_terms <- grep("^chemokine receptor", all_go_terms, value = TRUE)

my_hematopoietic_cell_lineage_go_terms <- grep("hematopo", all_go_terms, value = TRUE)


my_go_terms <- c(my_cytokine_go_terms,
                 my_cytokine_receptor_go_terms,
                 my_chemokine_go_terms,
                 my_chemokine_receptor_go_terms,
                 my_hematopoietic_cell_lineage_go_terms
                 )

# get go id  of my terms (definition, and ontology)
my_go_ids <- AnnotationDbi::select(Homo.sapiens,
                                   keys=my_go_terms,
                                   columns = c("GOID"),#,"ONTOLOGY","DEFINITION"
                                   keytype = "TERM") %>% 
  pluck("GOID")

columns(Homo.sapiens)

head(keys(Homo.sapiens, keytype="TERM"))

# Get my genes symbols
annoDb.cytokine.chemokine.Gene.symb <- AnnotationDbi::select(Homo.sapiens,
                                         keys=my_go_ids,
                                         columns = c('SYMBOL'),
                                         keytype = "GOID") %>% 
  filter(GOID%in%c(my_go_ids)) %>% #"GO:0005125","GO:0008009"
  pluck("SYMBOL") %>% 
  unique()

# get my uniprot ids
annoDb.cytokine.chemokine.unipID <- AnnotationDbi::select(Homo.sapiens,
                                       keys=my_go_ids,
                                       columns = c('UNIPROT'),
                                       keytype = "GOID") %>% 
  filter(GOID%in%c("GO:0005125","GO:0008009")) %>% 
  pluck("UNIPROT") %>% 
  unique()

# Get my genes symbols and uniprot
annoDb.cytokine.chemokine <- AnnotationDbi::select(Homo.sapiens,
                                                             keys=my_go_ids,
                                                             columns = c('SYMBOL','UNIPROT','TERM'),
                                                             keytype = "GOID") %>% 
  filter(GOID%in%c(my_go_ids)) %>% #"GO:0005125","GO:0008009"
  dplyr::select(Gene.symbol="SYMBOL",UniprotID='UNIPROT',GO_term='TERM') %>% 
  distinct()%>% 
  filter(GO_term!="chemokine receptor transport within lipid bilayer") %>%  # this term had no genes in it
  filter(GO_term!="regulation of hematopoietic stem cell migration") %>% 
  filter(Gene.symbol!="NA")


annoDb.cytokine.chemokine %>% 
  filter(GO_term=="regulation of hematopoietic stem cell migration")

#############################
#
#############################

##############################
# Format data for rank plot
########################




setwd("/Users/xrydbh/OneDrive - University of Gothenburg/Neored")

# read ms file
ms.raw.0.0 <- readxl::read_excel("./data/001_ms_raw.xlsx") 


# make longer
mass_spec.0.3.3 <- ms.raw.0.0 %>% 
  pivot_longer(cols=contains("_"),names_to="SC_Donor",values_to = "Norm_abundance") %>% 
  mutate(Norm_abundance=as.numeric(Norm_abundance))

############################################
# Count numb of NAs in High and Low
# group by proteins
############################################
# Separate column SC_Donor into "High low Group" and "Unit 1-8"
mass_spec.0.3.4 <- mass_spec.0.3.3 %>%
  separate(SC_Donor,
           into = c("Group","Unit"),
           sep = "_")

# Get missing values per gene and group
missing_values_per_gene.0.0 <- mass_spec.0.3.4 %>% 
  group_by(UniprotID,Group) %>%           # group by gene symbol and High/Low
  summarise(Missing_ab_val = sum(is.na(Norm_abundance))) %>%  # get number of NAs per gene and High/Low groups
  pivot_wider(names_from = Group,values_from =Missing_ab_val) %>% 
  dplyr::rename(Numb_NA_High=High,Numb_NA_Low=Low)

mass_spec.0.3.5 <- full_join(mass_spec.0.3.3,missing_values_per_gene.0.0,by="UniprotID")



# summarize and arrange by Norm_abundance
mass_spec.for.rankplot.0.0 <- mass_spec.0.3.5 %>% 
  group_by(UniprotID) %>% # group by gene symbol
  summarize(sum_norm_ab=sum(Norm_abundance,na.rm = T)) %>% # summarize into sum per genes across samples
  filter(sum_norm_ab!=0) %>% # remove sums that are zero
  arrange(desc(sum_norm_ab)) 

mass_spec.for.rankplot <- mass_spec.for.rankplot.0.0 %>% # sort
  mutate(rank=seq(1:nrow(mass_spec.for.rankplot.0.0))) # crate rank variable

mass_spec.for.rankplot.gs.0.0 <-  mass_spec.0.3.5 %>% 
  group_by(Gene.symbol) %>% # group by gene symbol
  summarize(sum_norm_ab=sum(Norm_abundance,na.rm = T)) %>% # summarize into sum per genes across samples
  filter(sum_norm_ab!=0) %>% # remove sums that are zero
  arrange(desc(sum_norm_ab)) 

mass_spec.for.rankplot.gs <- mass_spec.for.rankplot.gs.0.0 %>% # sort
  mutate(rank=seq(1:nrow(mass_spec.for.rankplot.gs.0.0))) # crate rank variable


mass_spec.for.rankplot.gs2 <- mass_spec.for.rankplot %>% inner_join(mass_spec.0.3.5,by=c("UniprotID"="UniprotID")) %>% 
  dplyr::select(Gene.symbol,sum_norm_ab,rank) %>% 
  distinct()

# Get 5 topranked
mass_spec.for.rankplot %>% head(10)

###########################
# check overlap with Cytokines and chemokines
#######################
# Get cytokines as unip not in ms
cyts_in_ms.uni <- mass_spec.for.rankplot %>% 
  filter(UniprotID%in% annoDb.cytokine.chemokine$UniprotID,) 

# Get cytokines as gs in ms via uniprot
cyts_in_ms.gs <- mass_spec.for.rankplot %>% 
  filter(UniprotID %in% annoDb.cytokine.chemokine$UniprotID,) %>% # Filter ms data to only include cytokines and chemokines (Using uniportid as key)
  inner_join(mass_spec.0.3.5,by=c("UniprotID"="UniprotID")) %>% # replace uniprotID with Gene symbol
  dplyr::select(Gene.symbol,sum_norm_ab,rank,Numb_NA_High,Numb_NA_Low) %>% #also keep NA counts 
  distinct()

# Get cytokines as gs in ms 
cyts_in_ms.gs <- mass_spec.for.rankplot.gs %>% 
  filter(Gene.symbol %in% annoDb.cytokine.chemokine$Gene.symbol,) 

# Get cytokines as gs in ms 
cyts_in_ms.gs3 <- mass_spec.for.rankplot.gs %>% 
  inner_join(annoDb.cytokine.chemokine,by=c("Gene.symbol"="Gene.symbol")) %>% 
  dplyr::select(-UniprotID) %>% 
  distinct() %>% 
  mutate(GO_term=replace(GO_term, GO_term=="chemokine receptor antagonist activity", "chemokine receptor binding"))
cyts_in_ms.gs3$GO_term %>% as.factor() %>% levels()
    
#    filter(Gene.symbol %in% annoDb.cytokine.chemokine$Gene.symbol,) 

# get number of genes per go terms
numbGenesInGO <- annoDb.cytokine.chemokine %>% 
  pluck("Gene.symbol") %>% 
  unique() %>% 
  length()

# get number of genes per go terms
goTermCounts <- annoDb.cytokine.chemokine %>% 
  dplyr::select(-UniprotID) %>% 
  distinct() %>% 
  group_by(GO_term) %>% 
  tally()

# get number of gene matching in ms data
numbMSGenesInGO <- cyts_in_ms.gs3 %>% 
  pluck("Gene.symbol") %>% 
  unique() %>% 
  length()

percCytsInMS <- (numbMSGenesInGO/numbGenesInGO)*100


# get number of gene per go term matching in ms data
goTerminMsCounts <-cyts_in_ms.gs3 %>% 
  distinct %>% 
  group_by(GO_term) %>% 
  tally()

# join
cytokine_summary_table.0.0 <- goTermCounts %>% 
  full_join(goTerminMsCounts,by=c("GO_term"="GO_term"))

# replace NA with 0
cytokine_summary_table.0.0$n.y <- cytokine_summary_table.0.0$n.y %>% 
   tidyr::replace_na(0)

# get percenttge
cytokine_summary_table.0.1 <- cytokine_summary_table.0.0 %>% 
  dplyr::rename(nInAnnoDB=n.x,nInMS=n.y) %>% 
  mutate(perc_in=(nInMS/nInAnnoDB)*100) %>% 
  mutate(perc_out=((nInAnnoDB-nInMS)/nInAnnoDB)*100)

# listGO_term gene members not in MS
GO_term_gene_members_not_in_ms <- annoDb.cytokine.chemokine %>% 
  dplyr::filter(!(Gene.symbol%in%cyts_in_ms.gs3$Gene.symbol)) %>% 
  dplyr::select(-UniprotID) %>% 
  distinct() %>% 
  dplyr::group_by(GO_term) %>% 
  dplyr::summarise(Gene.symbol = paste(Gene.symbol, collapse = ","))

# list  GO_term MS gene members
GO_term_ms_gene_members <- cyts_in_ms.gs3 %>% 
  dplyr::group_by(GO_term) %>% 
  dplyr::summarise(Gene.symbol = paste(Gene.symbol, collapse = ","))

cytokine_summary_table.0.2 <- cytokine_summary_table.0.1 %>% 
  full_join(GO_term_ms_gene_members,by=c("GO_term"="GO_term")) %>% 
  dplyr::rename(In_ms=Gene.symbol) %>% 
  full_join(GO_term_gene_members_not_in_ms,by=c("GO_term"="GO_term")) %>% 
  dplyr::rename(Not_in_ms=Gene.symbol) 

# Get cytokines as gs not in ms
cyts_not_in_ms.gs <- annoDb.cytokine.chemokine %>% 
  filter(!c(UniprotID %in% mass_spec.for.rankplot$UniprotID)) %>% 
  pluck("Gene.symbol") %>% 
  unique()


# Get cytokines as unip not in ms
cyts_not_in_ms.uni <- annoDb.cytokine.chemokine %>% 
  filter(!c(UniprotID %in% mass_spec.for.rankplot$UniprotID)) %>% 
  pluck("UniprotID") %>% 
  unique()

numb_chemk_unip <- annoDb.cytokine.chemokine$UniprotID %>% 
  unique() %>% 
  length()
numb_chemk_gs <- annoDb.cytokine.chemokine$Gene.symbol %>% 
  unique() %>% 
  length()

numb_chemk_unip
numb_chemk_gs
cyts_in_ms.gs %>% nrow()
cyts_in_ms.uni %>% nrow()
cyts_not_in_ms.gs %>% length()
cyts_not_in_ms.uni %>% length()

cyts_in_ms.uni %>% nrow()+cyts_not_in_ms.uni %>% length()
cyts_in_ms.gs %>% nrow()+cyts_not_in_ms.gs %>% length()

  
mass_spec.for.rankplot %>% 
filter(UniprotID %in% annoDb.cytokine.chemokine$UniprotID,) %>% # Filter ms data to only include cytokines and chemokines (Using uniportid as key)
inner_join(mass_spec.0.3.5,by=c("UniprotID"="UniprotID")) %>% # replace uniprotID with Gene symbol
dplyr::select(Gene.symbol,sum_norm_ab,rank,Numb_NA_High,Numb_NA_Low) %>% #also keep NA counts 
distinct()


# Get missing cytokines in ms data (Gene symbol) 
# # Have to figure outwhat is happening here
# wtf_these_are_also_cytokines <-  annoDb.cytokine.chemokine.unipID[!c(annoDb.cytokine.chemokine$UniprotID %in% mass_spec.for.rankplot$UniprotID)] %>% # subset annoDB list, by uniprots not in ms
#   as_tibble() %>% # make tibble 
#   rename("UniprotID"="value") %>% # rename variable
#   inner_join(mass_spec.0.3.5,by=c("UniprotID"="UniprotID")) %>% 
#   dplyr::select(Gene.symbol,sum_norm_ab,rank) %>% 
#   distinct()


#####################


# 
# mssaadfs <- mass_spec.0.3.3 %>% 
#   dplyr::select(UniprotID,Gene.synonyms)
# 
# inner_join(ms_only_cyts, mssaadfs, by=c("UniprotID"="UniprotID")) %>% 
#   dplyr::select(Gene.synonyms) %>% 
#   unique()
# 
# library(viridis)

nb.cols <- 12
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
mycolors[11] <- "darkgrey"

# plot with some names
rankplot.ms <- mass_spec.for.rankplot.gs %>% 
  ggplot(aes(x=rank, y=log2(sum_norm_ab))) + 
  geom_point(color="grey",alpha=0.5) +
  #geom_point(data=cyts_in_ms.gs,color="red") +
  ggrepel::geom_text_repel(aes(x = rank, 
                      y = log2(sum_norm_ab),
                      color=GO_term, 
                      label = Gene.symbol,
                      fontface = "bold"),
                      data=cyts_in_ms.gs3,
                      max.overlaps=20,
                      nudge_x=90
                      ) +  
  scale_color_manual(values = mycolors)+
  #scale_color_brewer(palette="Paired")
 # scale_color_brewer(values = "Paired") + #    scale_color_manual(values = rainbow(13)) +  #

#geom_text(aes(label=Gene.symbol,color=as.factor(GO_term)),data=cyts_in_ms.gs3) + #, check_overlap = TRUE, nudge_x = 70
geom_text(aes(label=Gene.symbol), nudge_x = -90,data=head(mass_spec.for.rankplot.gs,n=20L),check_overlap = TRUE,alpha=0.5)+ #, nudge_x = -90,,position=position_jitter(width=1,height=1)
geom_text(aes(label=Gene.symbol),data=tail(mass_spec.for.rankplot.gs,n=20L),check_overlap = TRUE, nudge_x = -90,alpha=0.5)+
  labs(title = paste("Plot of norm_ab sum and rank of ",numbGenesInGO," MS proteins",sep=""),
       subtitle = paste("Demarkating ",numbMSGenesInGO,"(",round(percCytsInMS,2),"%) cytokines",sep=""),
       caption = "")
# ggtitle(paste("Plot of norm_ab sum and rank of ",numbGenesInGO," MS proteins",sep="") )
# ggsubtitle(paste("Demarkating ",numbMSGenesInGO,"(",percCytsInMS,") MS proteins",sep="") )

#percCytsInMS <- (numbMSGenesInGO/numbGenesInGO)*100

rankplot.ms

pdf(file = "./out_r/15_a_rankplot_concentrations_ms.pdf",   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 15) # The height of the plot in inches
rankplot.ms
dev.off()

write_xlsx(cytokine_summary_table.0.2,"./out_r/015a_cytokine_summary_table_ms.xlsx")

save(goTermCounts,annoDb.cytokine.chemokine,file="./RData/15_annoDb.cytokine.chemokine.RData")
