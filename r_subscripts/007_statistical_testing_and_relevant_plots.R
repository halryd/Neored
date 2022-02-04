library(tidyverse)
library(broom)
library(writexl)
library(DEqMS)
library(ggrepel)
library('EnhancedVolcano')
library(readxl)
library(sjmisc)

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).



p_val_cut_off=0.05
fold_change_cutoff= 2


# Get abundances
loaded_objs <- load("./RData/003_ab_list.RData")
loaded_objs


########################################################
# DEqMS; MS
######################################################
dat.ms <- ab_list[["dat_abundances"]][["ms_abundances"]] %>% 
  dplyr::select(-UniprotID) %>% 
  mutate_all(funs(as.numeric(.))) %>% 
  as.data.frame()

row.names(dat.ms) <- ab_list[["dat_abundances"]][["ms_abundances"]] %>% 
  pluck("UniprotID")

# boxplot(dat.ms,las=2,main="MS")

dat.ms = equalMedianNormalization(dat.ms)

# boxplot(dat.ms,las=2,main="MS")

cond = as.factor(c(rep("High",8),rep("Low",8)))

# The function model.matrix is used to generate the design matrix
design = model.matrix(~0+cond) # 0 means no intercept for the linear model
colnames(design) = gsub("cond","",colnames(design))

# you can define one or multiple contrasts here
x <- c("High-Low")
contrast =  makeContrasts(contrasts=x,levels=design)
fit1 <- lmFit(dat.ms, design)
fit2 <- contrasts.fit(fit1,contrasts = contrast)
fit3 <- eBayes(fit2)


# Get PSMs from ms metadata
psm.0.1 <- ab_list[["dat_meta"]][["ms_meta"]] %>% 
  dplyr::select(UniprotID,PSMs) %>% 
  mutate(PSMs=as.numeric(PSMs)) %>% 
  as.data.frame()
row.names(psm.0.1) <- ab_list[["dat_meta"]][["ms_meta"]]$UniprotID


fit3$count = psm.0.1[rownames(fit3$coefficients),"PSMs"] 
fit4 = spectraCounteBayes(fit3)

# Visualize the fit curve - variance dependence on quantified PSM
# n=30 limits the boxplot to show only proteins quantified by <= 30 PSMs.
VarianceBoxplot(fit4,n=30,main="MS",xlab="PSM count")

VarianceScatterplot(fit4,main="MS")

# Extract the results as a data frame
DEqMS.results = outputResult(fit4,coef_col = 1) 

# DEqMS.results.uni <- DEqMS.results %>% 
#   as_tibble()
 
# if you are not sure which coef_col refers to the specific contrast,type
head(fit4$coefficients)

prot_ID_key_ms <- ab_list[["dat_meta"]]$ms_meta %>% 
  dplyr::select( UniprotID,Gene.symbol,Gene.synonyms)

DEqMS.results <- DEqMS.results %>% 
  left_join(prot_ID_key_ms,by= c("gene"="UniprotID")) 

rownames(DEqMS.results) <- DEqMS.results$Gene.symbol
  

# write_xlsx(DEqMS.results,
#            path = paste("./out_r/007_DEqMS_res.xlsx",sep="")) # write with contamination instead

DEqMS.sign.results.5perc.ms <- DEqMS.results %>%
  dplyr::filter(sca.P.Value<p_val_cut_off)

DEqMS.sign.results.5perc.FC.5.ms <- DEqMS.results %>%
  dplyr::filter(sca.P.Value<p_val_cut_off) %>%
  dplyr::filter(abs(logFC) > log2(fold_change_cutoff))

dim(DEqMS.sign.results.5perc.ms)
dim(DEqMS.sign.results.5perc.FC.5.ms)

write_xlsx(DEqMS.sign.results.5perc.FC.5.ms,
           path = paste("./out_r/007_DEqMS_sign_res_p005_FC_",fold_change_cutoff,"_ms.xlsx",sep=""))


DEqMS.results$log.sca.pval = -log10(DEqMS.results$sca.P.Value)

fit4$p.value = fit4$sca.p


uniprotIDs <- rownames(fit4$coefficients)
Gene.Symbols <- prot_ID_key_ms %>% 
  filter(UniprotID%in%uniprotIDs) %>% 
  pluck("Gene.symbol")

pdf(file = paste("./out_r/007_DEqMS_volcano.pdf",sep=""))
print(EnhancedVolcano(DEqMS.results,
                      lab = rownames(DEqMS.results),
                      x = 'logFC',
                      y = 'sca.P.Value',
                      xlim = c(-3,3),
                      ylim = c(0,4.5),
                      labSize = 4.0,
                      FCcutoff = log2(fold_change_cutoff),
                      pCutoff = p_val_cut_off, 
                      title = 'Volcano plot MS DEqMS',
                      subtitle = paste('P-val cut off: ',p_val_cut_off,', FC cut off: ',fold_change_cutoff,sep="")
)
)# ,selectLab=top10_norm_ab
dev.off()

################################ Demarcates possible contaminant
# read file with clinical data
# Crete function to rad multiple Excel sheets
read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

# Read quality marker categories into list
contaminants.0.0 <- read_excel_allsheets("./data/quality_markers.xlsx") #samples as rows
contaminants.0.0 %>% names()

# contaminants.0.0$erythr %>% 
#   #select("Protein IDs",`Gene names`) %>% 
#   as_tibble()

# Convert list into data frame and select columns
contaminants.0.1 <- contaminants.0.0 %>% 
  bind_rows(.id = "Contamination.type") %>% 
  dplyr::select(UniprotID="Protein IDs",Possible.contamination="Contamination.type") %>% 
  as_tibble() %>% 
  mutate(UniprotID=strsplit(as.character(UniprotID), ";")) %>%
  unnest(UniprotID) 
  
# Annotate DEqMS with potential contaminants
DEqMS.results.cont.0.0 <- DEqMS.results %>% 
  left_join(contaminants.0.1, by=c("gene"="UniprotID")) %>% 
  mutate(Possible.contamination = replace(Possible.contamination, str_detect(Gene.symbol,"KRT"), "keratin")) %>% 
  mutate(Possible.contamination = replace(Possible.contamination, str_detect(Gene.symbol,"LPA"), "Known genetic high variabilty")) %>% 
  mutate(Possible.contamination = replace(Possible.contamination, str_detect(Gene.symbol,"HLA"), "Known genetic high variabilty")) 

#Easier solution  
# DEqMS.results.cont.0.0 %>% tidyr::replace_na(list(Possible.contamination = "Not a contaminant"))

DEqMS.results.cont.0.1 <- DEqMS.results.cont.0.0 %>% 
  dplyr::mutate(Possible.contamination=tidyr::replace_na(Possible.contamination,"Not a contaminant")) #$Contamination %>% replace_na("Not a contaminant")

dim(DEqMS.results.cont.0.1)

# Collect contamination anotatin inot coma separated list
DEqMS.results.cont.0.2 <- DEqMS.results.cont.0.1 %>% 
  group_by(gene) %>% 
  dplyr::mutate(Possible.contamination = paste(Possible.contamination, collapse = ","))

dim(DEqMS.results.cont.0.2)

DEqMS.results.cont.0.0$Possible.contamination %>% 
  is.na() %>% 
  head()

df <- DEqMS.results.cont.0.0 %>% 
  mutate_if(is.character, 
                str_replace_all, pattern = "<NA>", replacement = "newVar")


write_xlsx(DEqMS.results.cont.0.2,
           path = paste("./out_r/007_DEqMS_res.xlsx",sep=""))


############################### End DEqMS
###############################
# ttest PEA
###############################


pea_long.0.0 <- ab_list[["dat_long"]][["pea_abs_long_eq_med_norm"]] %>% 
  dplyr::select(UniprotID,Stem.cell.cont,Sample,Prot_median_cent_rel_npx=value) #%>% 
  # separate(
  #   SC_Donor,
  #   c("Stem.cell.cont","Sample"),
  #   sep = "_")

group_variables <- c("UniprotID","Stem.cell.cont")

t_test_res_pea <- pea_long.0.0  %>% 
  group_by_at(group_variables) %>% 
  nest() %>% 
  spread(key = Stem.cell.cont, value = data) %>% 
  mutate(
    t_test = map2(High, Low, ~{t.test(.x[,"Prot_median_cent_rel_npx"], .y[,"Prot_median_cent_rel_npx"]) %>% tidy()}),
    High = map(High, nrow),
    Low = map(Low, nrow)
  ) %>%
  unnest(cols = c(High, Low, t_test)) %>% #,foldchange
  arrange(p.value) %>% 
  as.data.frame()

t_test_res_pea <- t_test_res_pea %>% 
  mutate(BH_adj=p.adjust(p.value, method = "BH", n = dim(t_test_res_pea)[1])) %>% 
  mutate(fdr_adj=p.adjust(p.value, method = "fdr", n = dim(t_test_res_pea)[1]))

prot_ID_key_pea <- ab_list[["dat_meta"]]$pea_meta %>% 
  dplyr::select( UniprotID,Gene.symbol)

t_test_res_pea <- t_test_res_pea %>% 
  left_join(prot_ID_key_pea,by= c("UniprotID"="UniprotID")) 

row.names(t_test_res_pea) <- t_test_res_pea$Gene.symbol

pdf(file = paste("./out_r/007_ttest_pea_volcano.pdf",sep=""))
print(EnhancedVolcano(t_test_res_pea,
                      lab = row.names(t_test_res_pea),
                      x = 'estimate',
                      y = 'p.value',
                      xlim = c(-2.5,2.5),
                      ylim = c(0,4.5),
                      labSize = 4.0,
                      FCcutoff = log2(fold_change_cutoff),
                      pCutoff = p_val_cut_off, title = 'Volcano plot PEA', 
                      subtitle = paste('P-val cut off: ',p_val_cut_off,', FC cut off: ',fold_change_cutoff,sep=""))
)# ,selectLab=top10_norm_ab
dev.off()

# write_xlsx(t_test_res_pea,
#            path = paste("./out_r/007_t_test_res_pea.xlsx",sep=""))


ttest.sign.results.5perc.pea <- t_test_res_pea %>%
  dplyr::filter(p.value<p_val_cut_off)

ttest.sign.results.5perc.FC.5.pea <- t_test_res_pea %>%
  dplyr::filter(p.value<p_val_cut_off) %>%
  dplyr::filter(abs(estimate) > log2(fold_change_cutoff))


dim(ttest.sign.results.5perc.pea)
dim(ttest.sign.results.5perc.FC.5.pea)

write_xlsx(ttest.sign.results.5perc.FC.5.pea,
           path = paste("./out_r/007_t_test_sign_res_p005_FC_",fold_change_cutoff,"_pea.xlsx",sep=""))

# Make count table
numbMS <- dim(dat.ms)[1]
numbPEA <- pea_long.0.0$UniprotID %>% unique() %>% length()
numbSign <- c(dim(DEqMS.sign.results.5perc.ms)[1],dim(ttest.sign.results.5perc.pea)[1])
percSign <- c(dim(DEqMS.sign.results.5perc.ms)[1]/numbMS,dim(ttest.sign.results.5perc.pea)[1]/numbPEA)
numbSignFC05 <- c(dim(DEqMS.sign.results.5perc.FC.5.ms)[1],dim(ttest.sign.results.5perc.FC.5.pea)[1])
percSignFC05 <- c(dim(DEqMS.sign.results.5perc.FC.5.ms)[1]/numbMS,dim(ttest.sign.results.5perc.FC.5.pea)[1]/numbPEA)
platform <- c("MS","PEA")

countsOfSignificant <- tibble(platform=platform,numbSign=numbSign,percSign=percSign,numbSignFC05=numbSignFC05,percSignFC05=percSignFC05 )


write_xlsx(countsOfSignificant,
           path = paste("./out_r/007_countsOfSignificant.xlsx",sep=""))

# Annotae PEA with potential contaminants
t_test_res_pea.cont.0.0 <- t_test_res_pea %>% 
  left_join(contaminants.0.1, by=c("Gene.symbol"="UniprotID")) %>% 
  mutate(Possible.contamination = replace(Possible.contamination, str_detect(Gene.symbol,"KRT"), "keratin")) %>% 
  mutate(Possible.contamination = replace(Possible.contamination, str_detect(Gene.symbol,"LPA"), "Known genetic high variabilty")) %>% 
  mutate(Possible.contamination = replace(Possible.contamination, str_detect(Gene.symbol,"HLA"), "Known genetic high variabilty")) 



t_test_res_pea.cont.0.1 <- t_test_res_pea.cont.0.0 %>% 
  mutate(Possible.contamination=tidyr::replace_na(Possible.contamination,"Not a contaminant")) #$Contamination %>% replace_na("Not a contaminant")

dim(t_test_res_pea.cont.0.1)

# Collect contamination anotatin inot coma separated list
# Not needed for PEA
# t_test_res_pea.cont.0.2 <- t_test_res_pea.cont.0.1 %>% 
#   group_by(UniprotID) %>% 
#   dplyr::mutate(Possible.contamination = paste(Possible.contamination, collapse = ","))
# 
# dim(t_test_res_pea.cont.0.2)

write_xlsx(t_test_res_pea.cont.0.1,
           path = paste("./out_r/007_t_test_res_pea.xlsx",sep=""))



#########################
signOFRegulationMS <- DEqMS.results %>%
  dplyr::filter(sca.P.Value<p_val_cut_off) %>% 
  mutate(property=case_when(logFC < 0 ~ "Negative", 
                            logFC > 0 ~ "Positive", 
                            TRUE  ~ "Zero")) %>% 
  dplyr::select(gene,logFC,property) %>% 
  group_by(property) %>% 
  tally() 

ms_tot_sign.0.0 <- signOFRegulationMS %>%  
  summarise(across(n, ~ sum(., is.na(.), 0))) %>% 
  t()

ms_tot_sign.0.1 <- c("ms_tot_sign",ms_tot_sign.0.0 )
#names(ms_tot_sign.0.0) <- c("signOfRegulation","n")

ms_tot_sign <- ms_tot_sign.0.1 %>% 
  t() %>% 
  as_tibble() %>% 
  mutate(V2=as.numeric(V2)) %>% 
  dplyr::select(property=V1,n=V2)



ms_tot <- c("ms_tot",row.names(dat.ms) %>% unique() %>% length()) %>% 
  t() %>% 
  as_tibble() %>% 
  mutate(V2=as.numeric(V2)) %>% 
  dplyr::select(property=V1,n=V2)

ms_prot_counts <- signOFRegulationMS %>% bind_rows(ms_tot_sign,ms_tot)


ms_prot_counts
perc <- c(ms_prot_counts$n[1]/ms_prot_counts$n[3],ms_prot_counts$n[2]/ms_prot_counts$n[3],ms_prot_counts$n[3]/ms_prot_counts$n[4],ms_prot_counts$n[4]/ms_prot_counts$n[4])

ms_prot_counts_and_perc <-  bind_cols(platform=rep("ms",nrow(ms_prot_counts)),ms_prot_counts,perc=perc)
################

# get nuber of up and down regualted
signOFRegulationPEA <- t_test_res_pea %>%
  dplyr::filter(p.value<p_val_cut_off) %>% 
  mutate(property=case_when(estimate  < 0 ~ "Negative", 
                            estimate  > 0 ~ "Positive", 
                            TRUE  ~ "Zero")) %>% 
  dplyr::select(UniprotID,p.value,property) %>% 
  group_by(property) %>% 
  tally()


# Get number of significant proteins
pea_tot_sign.0.0 <- signOFRegulationPEA %>%  
  summarise(across(n, ~ sum(., is.na(.), 0))) %>% 
  t()

pea_tot_sign.0.1 <- c("pea_tot_sign",pea_tot_sign.0.0 )
#names(pea_tot_sign.0.0) <- c("signOfRegulation","n")

pea_tot_sign <- pea_tot_sign.0.1 %>% 
  t() %>% 
  as_tibble() %>% 
  mutate(V2=as.numeric(V2)) %>% 
  dplyr::select(property=V1,n=V2)


pea_tot <- c("pea_tot",pea_long.0.0$UniprotID %>% unique() %>% length()) %>% 
  t() %>% 
  as_tibble() %>% 
  mutate(V2=as.numeric(V2)) %>% 
  dplyr::select(property=V1,n=V2)

pea_prot_counts <- signOFRegulationPEA %>% bind_rows(pea_tot_sign,pea_tot)


pea_prot_counts
perc <- c(pea_prot_counts$n[1]/pea_prot_counts$n[3],pea_prot_counts$n[2]/pea_prot_counts$n[3],pea_prot_counts$n[3]/pea_prot_counts$n[4],pea_prot_counts$n[4]/pea_prot_counts$n[4])

pea_prot_counts_and_perc <-  bind_cols(platform=rep("pea",nrow(pea_prot_counts)),pea_prot_counts,perc=perc)

pea_prot_counts_and_perc 

signOFRegualtionCountsMS_PEA <- ms_prot_counts_and_perc  %>% bind_rows(pea_prot_counts_and_perc)
names(signOFRegualtionCountsMS_PEA) <- c("platform","protein_property","n","perc")

write_xlsx(signOFRegualtionCountsMS_PEA,"./out_r/007_signOFRegualtionCountsMS_PEA.xlsx")




