library(tidyverse)
library(DescTools)
library(writexl)
library(ggvenn)
library(VennDiagram) 

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).



corr_method="spearman"
p.adjust.methods=c("none") #fdr,none,BH

# Get abundances
loaded_abundancies <- load("./RData/003_ab_list.RData")
loaded_abundancies 

ms_key <- ab_list[["dat_meta"]][["ms_meta"]] %>% 
  dplyr::select(UniprotID,Gene.symbol)

pea_key <- ab_list[["dat_meta"]][["pea_meta"]] %>% 
  dplyr::select(UniprotID,Gene.symbol)

# sampName_key <- load("./RData/001_alt_samp_names.RData")
# #MS_column_head.0.1$Plf_unique_no
# #PEA_column_head.0.1

PEA_column_head.0.1 <- readxl::read_excel("./data/PEA_column_head.xlsx") 
MS_column_head.0.1 <- readxl::read_excel("./data/MS_column_head.xlsx") 




# get ms abundances
ms_ab <- ab_list %>% 
  pluck("dat_abundances") %>% 
  pluck("ms_abundances_eq_med_norm") %>%  #"ms_abundances"
  dplyr::select(-Gene.symbol)
dim(ms_ab)


# Get pea abundances
pea_ab.0.0 <- ab_list %>%
  pluck("dat_abundances") %>%
  pluck("pea_abundances_eq_med_norm") %>%
  dplyr::select(-c(Panel))

# dim(pea_ab.0.0 )

# Get pea abundances; Need to filter on number of NA ; For shared we dont want any missing values; cannot use na.omit one PEA data since we are keeping themissing values
pea_ab_filt <- ab_list %>% 
  pluck("dat_meta") %>% 
  pluck("pea_meta") %>% 
  filter(Numb_NA_High<=0&Numb_NA_Low<=0) %>% 
  dplyr::select(names(pea_ab.0.0))
  
dim(pea_ab_filt)

pea_ab <- pea_ab.0.0 %>% 
  filter(UniprotID%in%pea_ab_filt$UniprotID) %>% 
  dplyr::select(-Gene.symbol)

dim(pea_ab)

# ##########################
# # Make Venn diagram
# ##########################
# pea_meta <- pea_ab.0.0$UniprotID
# ms_meta <- ms_ab$UniprotID
# l.pea <-pea_ab.0.0$UniprotID %>% length()
# l.ms <- ms_ab$UniprotID %>% length()
# 
# venn_list <- list(MS=ms_ab$UniprotID,PEA=pea_ab$UniprotID)
# 
# 
# # pdf(file = paste("./out_r/019_pea_overlap_wgcna_ttest.pdf",sep=""))
# # ggvenn(
# #   venn_list , 
# #   fill_color = c("#0073C2FF", "#EFC000FF"), #, "#CD534CFF",, "#868686FF"
# #   stroke_size = 0.5, set_name_size = 4
# # )
# # dev.off()
# 
# venn_text=paste("Overlap of ",
#                 l.ms,
#                 " ms and \n",
#                 l.pea," Olink pea proteins",
#                 sep="")
# 
# 
# venn.diagram(venn_list, main=venn_text,fill = c("darkgreen", "darkgrey"), 
#              main.cex=3, 
#              alpha = c(0.7, 0.7), 
#              lwd =1, filename="./out_r/004_venn_overlap_ms_pea_unip.tiff",
#              cex=4,label.col=c("grey","black","darkolivegreen4"),
#              cat.cex=3,cat.dist=-0.05,cat.col=c("chartreuse4","lightgrey"))#
# 
# # Gene symbol
# ms_ab <- ab_list %>% 
#   pluck("dat_abundances") %>% 
#   pluck("ms_abundances_eq_med_norm") %>%  #"ms_abundances"
#   dplyr::select(-UniprotID)
# pea_ab <- pea_ab.0.0 %>% 
#   filter(UniprotID%in%pea_ab_filt$UniprotID) %>% 
#   dplyr::select(-UniprotID)
# venn_list <- list(MS=ms_ab$Gene.symbol,PEA=pea_ab$Gene.symbol)
# venn.diagram(main=venn_text,venn_list, fill = c("darkgreen", "darkgrey"), 
#              main.cex=3, 
#              alpha = c(0.7, 0.7), 
#              lwd =1, filename="./out_r/004_venn_overlap_ms_pea_gs.tiff",
#              cex=4,label.col=c("grey","black","darkolivegreen4"),
#              cat.cex=3,cat.dist=-0.05,cat.col=c("chartreuse4","lightgrey"))#
# 
# # end Venn


# Give platform specific sample names
names(ms_ab) <- c(names(ms_ab)[1],(MS_column_head.0.1$Plf_unique_no))
names(pea_ab) <- c(names(pea_ab)[1],(PEA_column_head.0.1$Plf_unique_no))

# Inner join ms and pea to get shared
shared.0.0 <- ms_ab %>% 
  inner_join(pea_ab,by="UniprotID") %>% 
  na.omit()

dim(shared.0.0)

# Gather all abundance values into one column
shared.long <- shared.0.0 %>% 
  pivot_longer(cols=-1) %>% 
  mutate_at(vars(c("value")),as.numeric)

# gather abundance values into MS and PEA columns respectively
shared.semilong <- shared.0.0 %>% pivot_longer(cols=-1,
                                               names_to = c(".value","samp"),
                                               names_sep = "_") %>% 
  mutate_at(vars(c("MS","PEA")),as.numeric)

# Make numeric and perform correlation test
shared.0.2 <- shared.semilong %>% 
  group_by_at("UniprotID") %>% 
  mutate_at(vars("MS","PEA"), list(as.numeric)) %>%
  summarise(Corr = cor.test(MS,PEA)$estimate,
            p = cor.test(MS,PEA, method=corr_method)$p.value,
            p_adjust= p.adjust(p, method = p.adjust.methods, n = length("UniprotID")),#dim(shared.0.2)[1]
            log_p_adjust= log2(p_adjust),
            logp = log2(p),
            Sig = p < 0.05,
            Sig_adjust = p_adjust < 0.05,
            CI_low = SpearmanRho(x=MS, y = PEA, conf.level=0.95)[2], #cor.test(MS,PEA)$conf.int[1]
            CI_high = SpearmanRho(x=MS, y = PEA, conf.level=0.95)[3])

numb_shared_prots <- length(shared.0.2$UniprotID)

# Annotate results with Gene.symbol and make Unipor-Gene.symbol united name
shared.0.3 <- inner_join(shared.0.2,ms_key,by=c("UniprotID"="UniprotID")) %>% #with ms gene.symb
  relocate(Gene.symbol,.after=UniprotID) %>% 
  dplyr::rename(Gene.symbol.ms=Gene.symbol) %>% # rename to Gene.symbol
  inner_join(pea_key,by=c("UniprotID"="UniprotID")) %>% # with pea Gene.symbol
  relocate(Gene.symbol,.after=Gene.symbol.ms) %>% 
  dplyr::rename(Gene.symbol.pea=Gene.symbol) %>%  
  # filter(Gene.symbol.pea!="FAM3C") %>% # What is going on with FAM3C ?
  # filter(Gene.symbol.ms!="FAM3C") %>% 
  unite(UnitedIDs,Gene.symbol.ms:UniprotID,sep="_")


# Replace -inf vals 
max_val<- shared.0.3$logp[!is.infinite(shared.0.3$logp)] %>% max()

min_val <- shared.0.3$logp[!is.infinite(shared.0.3$logp)] %>% min()

shared.0.3$logp[is.infinite(shared.0.3$logp)] <- min_val
shared.0.3$log_p_adjust[is.infinite(shared.0.3$log_p_adjust)] <- min_val


title_plot_anders <- paste("Plot of ",corr_method," correlation for ",numb_shared_prots," shared genes, no missing values","\nAdjust method: ",p.adjust.methods,sep="")

shared_prots_eq_med_norm_ranked_correlation <- shared.0.3 %>% ggplot(aes(x=reorder(UnitedIDs, Corr), # !!as.symbol(protID)
                          y=Corr, color=-log_p_adjust)) + #logp,log_p_adjust
  geom_hline(yintercept=0, linetype="dashed") + 
  geom_point(aes(shape=Sig_adjust), size=5) + #Sig,Sig_adjust
  labs(colour = "p") + #-logp
  theme_classic() + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) + 
  geom_linerange(aes(ymax=CI_low, ymin=CI_high)) + 
  coord_flip() + 
  xlab("Proteins measured by both MS and PEA") +
  ggtitle(title_plot_anders)

shared_prots_eq_med_norm_ranked_correlation

# remove
# remove old files
file_to_remove <- paste("004_shared_prots_eq_med_norm_ranked_correlation_",corr_method,"_corr_",sep="")
junk.0 <- dir(path="./out_r/",  pattern=file_to_remove) # ?dir
junk <- paste("./out_r/",junk.0,sep="")
file.remove(junk) # ?file.remove

# Anders plot 
pdf(file = paste("./out_r/004_shared_prots_eq_med_norm_ranked_correlation_",corr_method,"_corr_",numb_shared_prots ,"_adjust_method_",p.adjust.methods,".pdf",sep=""), #_samp10_rem The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 12) # The height of the plot in inches
print(shared_prots_eq_med_norm_ranked_correlation)

dev.off()

write_xlsx(shared.0.3,paste("./out_r/004_corr_shared.xlsx",sep=""))

# shared_old <- c("O00533","O43405","O75326","O95445","P00441","P00915","P01034","P01130","P02452","P02760","P02786","P03951","P04275","P04792","P05121","P05154","P05164","P05362","P05543","P05556","P06681","P07237","P07339","P07359","P08253","P0DOY2","P10586","P10721","P11226","P11717","P12830","P13591","P14151","P14543","P14780","P15144","P17813","P17936","P18065","P19320","P20023","P22105","P23284","P24592","P24821","P27487","P30086","P33151","P35443","P35542","P36222","P49747","P52888","P98160","Q02487","Q12860","Q13332","Q13867","Q15113","Q15166","Q15485","Q15582","Q15828","Q16698","Q16853","Q4KMG0","Q6EMK4","Q86VB7","Q8IW75","Q8NBP7","Q92520","Q99497","Q99969","Q99972","Q9UM47","Q9Y5Y7","Q9Y624","Q9Y6N7")
# 
# shared_new <- shared.0.2$UniprotID %>% unique
# 
# # Not in old: "Q12884" "Q9HCB6" "Q9NQ79"
# shared_new[!(shared_new%in%shared_old)]
# # Not in new:"P00441" "P07339" Not in raw data file
# shared_old[!(shared_old%in%shared_new)]
# 
# # P07339 in Olink input
# CTSD
# OID00622
# 
# # P00441
# SOD1
# OID01222




