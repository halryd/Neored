# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).



rmarkdown::render("./r_subscripts/021_compile_supplementory_figures.Rmd", c("html_document", "pdf_document","word_document"),output_dir="./out_r",output_file=c("021_Cord_blood_proteomics_CD34_Supplementary_figures_2021_08_16","021_Cord_blood_proteomics_CD34_Supplementary_figures_2021_08_16","021_Cord_blood_proteomics_CD34_Supplementary_figures_2021_08_16"))


