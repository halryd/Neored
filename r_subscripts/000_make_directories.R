library(tidyverse)
library(readxl)

# Three version of sample names are used: 
# Unit_no like: 19900088
# Plf_unique_no like: MS_01 or PEA_01
# SC_Donor High_1
# A conversion key is made with this script

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).



# Create output directory structure, if not allready existing
dir.create("RData")
dir.create("RData/reactome_results")
dir.create("out_r/Reactome", recursive = TRUE)
dir.create("out_r/GO_BP", recursive = TRUE)
dir.create("out_r/WGCNA/ms", recursive = TRUE)
dir.create("out_r/WGCNA/pea", recursive = TRUE)


