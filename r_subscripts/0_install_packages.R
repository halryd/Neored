###########################
# Install form The Comprehensive R Archive Network (CRAN)
###########################

# First specify the packages of interest
baseR_packages <- c("tidyverse","writexl","broom","heatmap3","openxlsx","qgraph","rlist","sjmisc")

## Then load or install&load all
package.check <- lapply(
  baseR_packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE, repos="https://ftpmirror1.infania.net/mirror/CRAN/",source="mac.binary")
      library(x, character.only = TRUE)
    }
  }
)

###########################
# install from Bioconductor
###########################

# First install BiocManager

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12",update = F)

# Then specify the packages of interest
bioconductor_packages <- c("ReactomeGSA", "WGCNA", "EnhancedVolcano", "clusterProfiler", "pathview", "enrichplot","org.Hs.eg.db")

# Then install/update bioconductor_packages. This one has, as far as I know, to be commented out 
# when packages has been installed once to avoid spending time on re-installing
BiocManager::install(bioconductor_packages,update=F)


###########################
# install with Devtools
###########################

install.packages("devtools", repos="https://ftpmirror1.infania.net/mirror/CRAN/")
devtools::install_github("stephenturner/annotables")

#igraph??
#ggrepel??
baseR_packages
knitr::write_bib(c(baseR_packages,bioconductor_packages,"annotables"), "./data/neored_packages.bib")
  
  
  




