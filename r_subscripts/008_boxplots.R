#` ---
#` title: "Boxplots
#` author: "Halfdan Rydbeck"
#` date: "`r Sys.Date()`"
#` output: github_document: default #html_preview: false
#` ---
#` 

library(readxl)
library(tidyverse)
library(rmarkdown)

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).




protID="UniprotID"
p_val_cut_off=0.05
fold_change_cutoff=2

#' # Get significant genes
de_q_sign_res_p005_FC_2 <- read_xlsx("./out_r/007_DEqMS_sign_res_p005_FC_2_ms.xlsx")
top_ms <- de_q_sign_res_p005_FC_2$gene

ttest_sign_res_p005_FC_2  <- read_xlsx("./out_r/007_t_test_sign_res_p005_FC_2_pea.xlsx")
top_pea <- ttest_sign_res_p005_FC_2$UniprotID

# Load normalized data
ms_pea_data <- load("./RData/003_ab_list.RData")
names(ab_list)
names(ab_list[["dat_long"]])

#ms_abs_long_eq_med_norm

# Choose long data sets and Select and split columns (for ms eq_med_norm [on top of Prot_median_cent_l2] is selected)
ms_for_bp.0.0 <- ab_list[["dat_long"]][["ms_abs_long_eq_med_norm"]] %>% 
  dplyr::select(UniprotID,Gene.symbol,Stem.cell.cont, Sample,Prot_median_cent_l2_norm_ab_eq_med_norm=value) %>% 
  dplyr::filter(!!as.symbol(protID)%in%top_ms) %>% 
  mutate_at(vars(c("Prot_median_cent_l2_norm_ab_eq_med_norm")),as.numeric) # %>% 
  # separate(
  #   SC_Donor,
  #   c("Stem.cell.cont","Sample"),
  #   sep = "_")

ms_key <- ms_for_bp.0.0 %>% 
  dplyr::select(UniprotID,Gene.symbol) %>% 
  distinct()

pea_for_bp.0.0 <- ab_list[["dat_long"]][["pea_abs_long_eq_med_norm"]] %>% 
  dplyr::select(UniprotID,Gene.symbol,Stem.cell.cont,Sample,Prot_median_cent_rel_npx_eq_med_norm=value) %>% 
  dplyr::filter(!!as.symbol(protID)%in%top_pea) %>% 
  mutate_at(vars(c("Prot_median_cent_rel_npx_eq_med_norm")),as.numeric) #%>% 
  # separate(
  #   SC_Donor,
  #   c("Stem.cell.cont","Sample"),
  #   sep = "_")

###################
# MS
###################
# There are genes where one Uniprot codes to multiple Gene.symbols; look for them
Uniprot_with_multiple_matches_in_GS_ms <- inner_join(ms_for_bp.0.0,ms_key,by=c("UniprotID"="UniprotID")) %>%
  #relocate(gene.symb,.after=UniprotID) %>%
  dplyr::rename(Gene.symbol.ms=Gene.symbol.x) %>%
  group_by(UniprotID) %>%
  dplyr::select(UniprotID,Gene.symbol.ms) %>%
  distinct() %>%
  dplyr::filter(n()>1) %>%
  arrange(UniprotID)

Uniprot_with_multiple_matches_in_GS_ms <- ms_for_bp.0.0 %>%
  group_by(UniprotID) %>%
  dplyr::select(UniprotID,Gene.symbol) %>%
  distinct() %>%
  dplyr::filter(n()>1) 

# Remove duplicatd gene.IDs
ms_for_bp <- ms_for_bp.0.0 %>% 
  dplyr::filter(Gene.symbol!="SNHG28") %>% # Removing multiple matches
  #dplyr::filter(Gene.symbol!="VSIG8") %>% #Dont remove this one RPIA (a Gene symbol synonym to corresponding to Uniprot P49247) is missing in this data set
  #dplyr::filter(Gene.symbol!="LOC101060545") %>% # or this one since since C1orf204 (a third Gene.symbol synonym ot Uniprot ID P0DPA2)
  unite(UnitedIDs,Gene.symbol:UniprotID,sep="_", remove = F)

# Check again; should have 0 rows
ms_for_bp %>%
  group_by(UniprotID) %>%
  dplyr::select(UniprotID,Gene.symbol) %>%
  distinct() %>%
  dplyr::filter(n()>1) %>%
  arrange(UniprotID)

# Find proteins with outliers
outliers_ms <- ms_for_bp %>% 
  group_by(UnitedIDs) %>% 
  dplyr::filter(Prot_median_cent_l2_norm_ab_eq_med_norm>3.1 | Prot_median_cent_l2_norm_ab_eq_med_norm < -3.1) %>% 
  pluck("UnitedIDs") %>% 
  unique()


n_prots=ms_for_bp %>% 
  dplyr::select(UnitedIDs,Gene.symbol) %>% 
  distinct() %>% 
  arrange(UnitedIDs) %>% 
  nrow()

n_prots

ms_for_bp_regular <- ms_for_bp %>%
  dplyr::filter(!(UnitedIDs%in%outliers_ms))

ms_for_bp_outliers<- ms_for_bp %>%
  dplyr::filter(UnitedIDs%in%outliers_ms)


n_prots_regular=ms_for_bp_regular %>% 
  dplyr::select(UnitedIDs,Gene.symbol) %>% 
  distinct() %>% 
  arrange(UnitedIDs) %>% 
  nrow()

n_prots_regular

n_prots_outliers=ms_for_bp_outliers %>% 
  dplyr::select(UnitedIDs,Gene.symbol) %>% 
  distinct() %>% 
  arrange(UnitedIDs) %>% 
  nrow()

n_prots_outliers


ms_for_bp_regular_plot <- ms_for_bp_regular  %>% ggplot(aes(x=Stem.cell.cont, 
                                  y=Prot_median_cent_l2_norm_ab_eq_med_norm,
                                  fill=Stem.cell.cont)) + 
  facet_wrap("UnitedIDs",ncol=6) + #vars(Gene.symbol)
  #theme() +
  geom_boxplot(outlier.size = 0) +
  geom_point(position=position_jitterdodge(),alpha=0.5, size=1) +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_classic() +
  theme(legend.position = "none",strip.text.x = element_text(size = 6)) + # , colour = "orange"
  ggtitle(paste("Mass spectrometry\n",n_prots_regular," with narrow distribution out of ",n_prots," proteins with p < ",p_val_cut_off,"; FC > ",fold_change_cutoff,sep="")) +
  labs(y= "Median centred log2 normalized abundance", x = "Top table protein")

ms_for_bp_regular_plot 

# Narrow distributinons
pdf(file = paste("./out_r/008_boxplot_p005_FC_2_sign_DEqMS_narrow_dist_eq_med_norm.pdf",sep=""))
ms_for_bp_regular_plot 
dev.off()

# Wider distributions
pdf(file = paste("./out_r/008_boxplot_p005_FC_2_sign_DEqMS_wider_dist_eq_med_norm.pdf",sep=""))
ms_for_bp_outliers %>% ggplot(aes(x=Stem.cell.cont, 
                                  y=Prot_median_cent_l2_norm_ab_eq_med_norm,
                                  fill=Stem.cell.cont)) + 
  facet_wrap("UnitedIDs",ncol=6) + #vars(Gene.symbol)
  #theme() +
  geom_boxplot(outlier.size = 0) +
  geom_point(position=position_jitterdodge(),alpha=0.5, size=1) +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_classic() +
  theme(legend.position = "none",strip.text.x = element_text(size = 6)) + # , colour = "orange"
  ggtitle(paste("Mass spectrometry\n",n_prots_outliers," with wider distribution out of ",n_prots," proteins with p < ",p_val_cut_off,"; FC > ",fold_change_cutoff,sep="")) +
  labs(y= "Median centred log2 normalized abundance", x = "Top table protein")
dev.off()

###################
# PEA
###################


Uniprot_with_multiple_matches_in_GS_pea <- pea_for_bp.0.0 %>% 
  group_by(UniprotID) %>% 
  dplyr::select(UniprotID,Gene.symbol) %>% 
  distinct() %>% 
  filter(n()>1) %>% 
  arrange(UniprotID)

pea_for_bp <- pea_for_bp.0.0 %>% 
  unite(UnitedIDs,Gene.symbol:UniprotID,sep="_", remove = F)

n_prots_pea=pea_for_bp %>% 
  dplyr::select(UnitedIDs,Gene.symbol) %>% 
  distinct() %>% 
  arrange(UnitedIDs) %>% 
  nrow()

# Make boxplot
pea_bp <- pea_for_bp %>% ggplot(aes(x=Stem.cell.cont, 
                          y=Prot_median_cent_rel_npx_eq_med_norm,
                          fill=Stem.cell.cont
                          )
                          ) + 
  facet_wrap("UnitedIDs",ncol=4) + #vars(Gene.symbol)
  geom_boxplot(outlier.size = 0) +
  geom_point(position=position_jitterdodge(),alpha=0.5, size=1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Olink") +
  ggtitle(paste("Olink\n",n_prots_pea," proteins with p < ",p_val_cut_off,"; FC > ",fold_change_cutoff,sep="")) +
  
  labs(y= "Median centred relative NPX (to sample sum)", x = "Top table protein")

# Boxplot for report
pea_bp

# Boxplot to file
pdf(file = paste("./out_r/008_boxplot_p005_FC_2_sign_ttest_pea_eq_med_norm.pdf",sep=""))
pea_bp
dev.off()

# ##########################
# Boxplot of shared significant
# ##########################
# Get significant genes
signMS_unip <- read_xlsx("./out_r/007_DEqMS_res.xlsx") %>% 
  dplyr::filter(sca.P.Value<p_val_cut_off) %>% 
  dplyr::select(UniprotID=gene)

sigPEA_unip <- read_xlsx("./out_r/007_t_test_res_pea.xlsx") %>% 
  dplyr::filter(p.value<p_val_cut_off) %>% 
  dplyr::select( UniprotID)


# get significant in both platforms
sign_in_both_plfs <- signMS_unip %>% inner_join(sigPEA_unip) %>% 
  pluck("UniprotID")


ms <- ab_list[["dat_long"]][["ms_abs_long_eq_med_norm"]] %>% 
  dplyr::select(UniprotID,Gene.symbol,Stem.cell.cont, Sample,Prot_median_cent_l2_norm_ab_eq_med_norm=value) %>% 
  mutate_at(vars(c("Prot_median_cent_l2_norm_ab_eq_med_norm")),as.numeric)%>% 
  filter(UniprotID%in%sign_in_both_plfs) %>% 
  mutate(platform="ms") %>% 
  mutate(eq_med_norm_ab=Prot_median_cent_l2_norm_ab_eq_med_norm,.keep = "unused")

ms_for_bp.0.0 
  

pea <- ab_list[["dat_long"]][["pea_abs_long_eq_med_norm"]] %>% 
  dplyr::select(UniprotID,Gene.symbol,Stem.cell.cont,Sample,Prot_median_cent_rel_npx_eq_med_norm=value) %>% 
  mutate_at(vars(c("Prot_median_cent_rel_npx_eq_med_norm")),as.numeric) %>% 
  filter(UniprotID%in%sign_in_both_plfs) %>% 
  mutate(platform="pea")%>% 
  mutate(eq_med_norm_ab=Prot_median_cent_rel_npx_eq_med_norm,.keep = "unused")

shared.long <- ms %>% bind_rows(pea)%>% 
  unite(UnitedIDs,Gene.symbol:UniprotID,sep="_", remove = F)



# png(paste("/Users/xrydbh/OneDrive/OneDrive\ -\ Göteborgs\ Universitet/Diff_protExp_in_highVsLow_Cd34_cordPlasma/plots_r/both/shared_prots/",abund_trans,"_boxplot_",used_pars_string,".png",sep=""), 3000,2000,res=300)
# print(ggplot(shared.long, aes(x = Gene.symbol,y= eq_med_norm_ab, fill=as.factor(Gene.symbol))) +
#         geom_boxplot(alpha = .3) +
#         theme(axis.text.x=element_text(angle=60, hjust=1,)) +
#         theme(legend.position = "none") + 
#         ggtitle("my_title"))
# 
# dev.off()


for (prot in sign_in_both_plfs){
  ms_pea_for_bp_one_gene <- shared.long %>% filter(UniprotID==prot)
  unified_prot_name <- ms_pea_for_bp_one_gene[1,]  %>% pluck("UnitedIDs")
  pdf(file = paste("./out_r/008_boxplot_shared_ms_pea_UnitedIDs_",prot,".pdf",sep=""), 
      width=3, 
      height=3)
  print(
    ms_pea_for_bp_one_gene %>% ggplot(aes(x=platform, 
                                          y=eq_med_norm_ab,
                                          fill=Stem.cell.cont,alpha=platform)) + 
      facet_wrap("Stem.cell.cont") + #vars(Gene.symbol),ncol=4
      geom_boxplot(outlier.size = 0) +
      geom_point(position=position_jitterdodge(),alpha=0.5, size=1) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      theme_classic() +
      theme(legend.position = "none") +
      ggtitle(paste("Significant for both platforms \n", unified_prot_name ,sep="")) +
      
      labs(y= "Eq_med_norm_ab", x = "Top table protein")
  )
  
  dev.off()
}  

# shared.0.0 <- sel_transform_alts.l[["ms"]] %>% 
#   inner_join(sel_transform_alts.l[["pea"]],by=protID)  %>% 
#   na.omit()

# ##########################
# # end  Boxplot of shared significant
# ##########################