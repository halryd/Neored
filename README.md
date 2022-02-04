R scripts for analysis of mass spectrometry and Olink PEA proteomics
data from cord blood
================

-   [What this workflow does](#what-this-workflow-does)
-   [The purpose of this repository](#the-purpose-of-this-repository)
-   [How to run the analysis (the
    r-scripts)](#how-to-run-the-analysis-the-r-scripts)
-   [Take a look the code](#take-a-look-the-code)
-   [Links to Reactome webserver to browse
    results](#links-to-reactome-webserver-to-browse-results)
-   [Generating the README.md file for Github (the file that you are
    reading
    now)](#generating-the-readmemd-file-for-github-the-file-that-you-are-reading-now)
-   [Generating the README.html for local
    browsing](#generating-the-readmehtml-for-local-browsing)
-   [References to all used R
    packages](#references-to-all-used-r-packages)

Supporting repository for the manuscript “The proteome signature of cord
blood plasma with high hematopoietic stem and progenitor cell count”

## What this workflow does

<!-- ```{r sunburst, echo=F, message=FALSE, fig.cap ="Three main types of analyses"} -->
<!-- source("./r_subscripts/make_graphcal_abstract/sunburst_plot.R") -->
<!-- fig -->
<!-- ``` -->
<!-- ![](./out_r/023_neored_graphical_ab_sunburst.png) -->

![](./data/023_neored_graphical_ab_sunburst.png)

This repository contains a set of R scripts that performs protein,
pathway and correlation module centric analysis of proteomics data.

More specifically it uses data from two different technologies that to a
large extent complement each other in terms of covered proteins but also
to a small extent overlaps in terms of that coverage. It takes two Excel
files as input. One with Olink NPX protein expression values and one
with mass spectrometry normalized abundance (formated output from
Proteome Discoverer).

The specific data analyzed in this project is collected from cord blood
plasma. 8 samples with high CDC34 concentration and 8 samples with low
such concentration. The project is trying to find proteomics biomarkers
for CD34 concentration. It looks for biomarkers on the single protein
level but also looks more widely into pathways with differing expression
patterns.

It looks for differentially expressed genes using the packages t-test
and DEqMS, pathways using ReactomeGSA and correlation modules using
WGCNA.

## The purpose of this repository

The main purpose of this repository is to make the published analysis
reproducible. A secondary usage would be for anyone to run the code for
another project with similar datasets. Hope fully the repository will
also develop beyond its state at publication. It applies a set for
recently published algorithms and software proven to outperform their
predecessors in a workflow to identify biomarkers.

## How to run the analysis (the r-scripts)

There is a README.Rmd. It can be used to run all or a subset of the R
scripts of the analysis. It also generates the README.md file for github
site.

To reproduce the published analysis:

1.  clone this repository
2.  Checkout commit 778d19a5f01fe38871da7a947e42700041e4886b from the
    master branch
3.  Move or delete the out\_r folder and make a new empty out\_r
4.  Edit “path\_to\_my\_project” in README.Rmd
5.  Rerun the scripts in the “r\_subscripts” folder. This is easiest
    done by running the README.Rmd file. It runs the R scripts and
    render markdown files and generates a README.md file with links to
    the markdown reports markdown files. (This was the only way I could
    render markdown form the scripts and keep the main direcory as
    working directory. I am sure there is a tidyer way to make README
    fils) Note: There is also a README\_local\_html.Rmd that can be used
    to generate a corresponding README file in html that can be used
    locally, not at github.)
6.  For each subscript to be run from the README file(s) it might have
    to be “uncommented” by removing the hash in front of it.

## Take a look the code

If you would like to see what is going on in the scripts the code and
output can be accessed with the linkes below. The links are in order of
executions.

-   [install packages](markdown/0_install_packages.md)
-   [ms\_project\_specific\_formatting](markdown/000a_ms_project_specific_formatting.md)
-   [pea\_project\_specific\_formatting](markdown/000b_pea_project_specific_formatting.md)
-   [make\_per\_plf\_sample\_name\_keys](markdown/001_make_per_plf_sample_name_keys.md)
-   [filter\_ms\_go\_bp\_keratinization.](markdown/002a1_filter_ms_go_bp_keratinization.md)
-   [ms\_format\_normalize\_naFilter](markdown/002a2_ms_format_normalize_naFilter.md)
-   [pea\_format\_normalize\_naFilter](markdown/002b_pea_format_normalize_naFilter.md)
-   [collect\_ms\_and\_pea](markdown/003_collect_ms_and_pea.md)
-   [correlation\_platform\_shared](markdown/004_correlation_platform_shared.md)
-   [gsea\_padog\_ReactomeGSA](markdown/005_gsea_padog_ReactomeGSA.md)
-   [collect\_reactome\_results\_tables](markdown/006a_collect_reactome_results_tables.md)
-   [visualize\_reactome\_results\_in\_dot\_plot](markdown/006b_reactome_dotplot.md)
-   [statistical\_testing\_and\_relevant\_plots](markdown/007_statistical_testing_and_relevant_plots.md)
-   [boxplots](markdown/008_boxplots.md)
-   [WGCNA using gene symbol gene ID](markdown/009_WGCNA.md)
-   [WGCNA using uniprot gene ID](markdown/009b_WGCNA_uniprot.md)
-   [plot\_WGCNA\_hubgene\_network](markdown/010_plot_WGCNA_hubgene_network.md)
-   [WGCNA\_ORA\_Kegg](markdown/011_WGCNA_ORA_Kegg.md)
-   [collect\_filt\_norm\_data\_and\_testresults](markdown/012_collect_filt_norm_data_and_testresults.md)
-   [add\_annotation](markdown/013_add_annotation.md)
-   [plot\_concentrations\_ms](markdown/015a_plot_concentrations_ms.md)
-   [plot\_concentrations\_pea](markdown/015b_plot_concentrations_pea.md)
-   [annotate\_relevant\_wgcna\_mods](markdown/016_annotate_relevant_wgcna_mods.md)
-   [compile\_supplementory\_file](markdown/020_compile_supplementory_tables.md)

<!-- When you want to create local html report uncommen this and comment the markdown references above -->
<!-- * [install packages](rendered_html/0_install_packages.html) -->
<!-- * [ms_project_specific_formatting](rendered_html/000a_ms_project_specific_formatting.html) -->
<!-- * [pea_project_specific_formatting](rendered_html/000b_pea_project_specific_formatting.html) -->
<!-- * [make_per_plf_sample_name_keys](rendered_html/001_make_per_plf_sample_name_keys.html) -->
<!-- * [ms_format_normalize_naFilter](rendered_html/002a_ms_format_normalize_naFilter.html) -->
<!-- * [pea_format_normalize_naFilter](rendered_html/002b_pea_format_normalize_naFilter.html) -->
<!-- * [collect_ms_and_pea](rendered_html/003_collect_ms_and_pea.html) -->
<!-- * [correlation_platform_shared](rendered_html/004_correlation_platform_shared.html) -->
<!-- * [gsea_padog_ReactomeGSA](rendered_html/005_gsea_padog_ReactomeGSA.html) -->
<!-- * [collect_reactome_results_tables](rendered_html/006_collect_reactome_results_tables.html) -->
<!-- * [statistical_testing_and_relevant_plots](rendered_html/007_statistical_testing_and_relevant_plots.html) -->
<!-- * [boxplots](rendered_html/008_boxplots.html) -->
<!-- * [WGCNA](rendered_html/009_WGCNA.html) -->
<!-- * [plot_WGCNA_hubgene_network](rendered_html/010_plot_WGCNA_hubgene_network.html) -->
<!-- * [WGCNA_ORA_Kegg](rendered_html/011_WGCNA_ORA_Kegg.html) -->
<!-- * [collect_filt_norm_data_and_testresults](rendered_html/012_collect_filt_norm_data_and_testresults.html) -->
<!-- * [add_annotation](rendered_html/013_add_annotation.html) -->
<!-- * [plot_concentrations_ms](rendered_html/015a_plot_concentrations_ms.html) -->
<!-- * [plot_concentrations_pea](rendered_html/015b_plot_concentrations_pea.html) -->
<!-- * [annotate_relevant_wgcna_mods](rendered_html/016_annotate_relevant_wgcna_mods.html) -->
<!-- * [compile_supplementory_file](rendered_html/020_compile_supplementory_tables.html) -->

## Links to Reactome webserver to browse results

Only active at Reactome server for seven days. Then this analysis (or
script 006) has to be rerun.

### Correlation Adjusted MEan RAnk (CAMERA)

-   [MS](https://www.reactome.org/PathwayBrowser/#/DTAB=AN&ANALYSIS=MjAyMjAyMDQyMDA3NDNfMjYzMg%3D%3D)
-   [PEA](https://www.reactome.org/PathwayBrowser/#/DTAB=AN&ANALYSIS=MjAyMjAyMDQyMDA3NTNfMjYzNA%3D%3D)

## Generating the README.md file for Github (the file that you are reading now)

Run/Knit “README.Rmd”

## Generating the README.html for local browsing

Run “render\_html.R”

## References to all used R packages

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-R-EnhancedVolcano" class="csl-entry">

Blighe, Kevin, Sharmila Rana, and Myles Lewis. 2020. *EnhancedVolcano:
Publication-Ready Volcano Plots with Enhanced Colouring and Labeling*.
<https://github.com/kevinblighe/EnhancedVolcano>.

</div>

<div id="ref-R-org.Hs.eg.db" class="csl-entry">

Carlson, Marc. 2020. *Org.hs.eg.db: Genome Wide Annotation for Human*.

</div>

<div id="ref-R-qgraph" class="csl-entry">

Epskamp, Sacha, Giulio Costantini, Jonas Haslbeck, and Adela Isvoranu.
2021. *Qgraph: Graph Plotting Methods, Psychometric Data Visualization
and Graphical Model Estimation*.
<https://CRAN.R-project.org/package=qgraph>.

</div>

<div id="ref-qgraph2012" class="csl-entry">

Epskamp, Sacha, Angélique O. J. Cramer, Lourens J. Waldorp, Verena D.
Schmittmann, and Denny Borsboom. 2012. “<span
class="nocase">qgraph</span>: Network Visualizations of Relationships in
Psychometric Data.” *Journal of Statistical Software* 48 (4): 1–18.

</div>

<div id="ref-R-ReactomeGSA" class="csl-entry">

Griss, Johannes. 2021. *ReactomeGSA: Client for the Reactome Analysis
Service for Comparative Multi-Omics Gene Set Analysis*.
<https://github.com/reactome/ReactomeGSA>.

</div>

<div id="ref-ReactomeGSA2020" class="csl-entry">

Griss, Johannes, Guilherme Viteri, Konstantinos Sidiropoulos, Vy Nguyen,
Antonio Fabregat, and Henning Hermjakob. 2020. “ReactomeGSA - Efficient
Multi-Omics Comparative Pathway Analysis.” *bioRxiv*.
<https://doi.org/10.1101/2020.04.16.044958>.

</div>

<div id="ref-R-WGCNA" class="csl-entry">

Langfelder, Peter, Steve Horvath with contributions by Chaochao Cai, Jun
Dong, Jeremy Miller, Lin Song, Andy Yip, and Bin Zhang. 2021. *WGCNA:
Weighted Correlation Network Analysis*.
<http://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/>.

</div>

<div id="ref-WGCNA2008" class="csl-entry">

Langfelder, Peter, and Steve Horvath. 2008. “WGCNA: An r Package for
Weighted Correlation Network Analysis.” *BMC Bioinformatics*, no. 1:
559.
<https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559>.

</div>

<div id="ref-WGCNA2012" class="csl-entry">

———. 2012. “Fast R Functions for Robust Correlations and Hierarchical
Clustering.” *Journal of Statistical Software* 46 (11): 1–17.
<https://www.jstatsoft.org/v46/i11/>.

</div>

<div id="ref-R-pathview" class="csl-entry">

Luo, Weijun. 2020. *Pathview: A Tool Set for Pathway Based Data
Integration and Visualization*. <https://pathview.uncc.edu/>.

</div>

<div id="ref-pathview2013" class="csl-entry">

Luo, Weijun, Brouwer, and Cory. 2013. “Pathview: An r/Bioconductor
Package for Pathway-Based Data Integration and Visualization.”
*Bioinformatics* 29 (14): 1830–31.
<https://doi.org/10.1093/bioinformatics/btt285>.

</div>

<div id="ref-sjmisc2018" class="csl-entry">

Lüdecke, Daniel. 2018. “Sjmisc: Data and Variable Transformation
Functions.” *Journal of Open Source Software* 3 (26): 754.
<https://doi.org/10.21105/joss.00754>.

</div>

<div id="ref-R-sjmisc" class="csl-entry">

———. 2021. *Sjmisc: Data and Variable Transformation Functions*.
<https://strengejacke.github.io/sjmisc/>.

</div>

<div id="ref-R-writexl" class="csl-entry">

Ooms, Jeroen. 2021. *Writexl: Export Data Frames to Excel Xlsx Format*.
<https://CRAN.R-project.org/package=writexl>.

</div>

<div id="ref-R-rlist" class="csl-entry">

Ren, Kun. 2021. *Rlist: A Toolbox for Non-Tabular Data Manipulation*.
<https://CRAN.R-project.org/package=rlist>.

</div>

<div id="ref-R-broom" class="csl-entry">

Robinson, David, Alex Hayes, and Simon Couch. 2022. *Broom: Convert
Statistical Objects into Tidy Tibbles*.
<https://CRAN.R-project.org/package=broom>.

</div>

<div id="ref-R-openxlsx" class="csl-entry">

Schauberger, Philipp, and Alexander Walker. 2021. *Openxlsx: Read, Write
and Edit Xlsx Files*. <https://CRAN.R-project.org/package=openxlsx>.

</div>

<div id="ref-R-annotables" class="csl-entry">

Turner, Stephen. 2022. *Annotables: Ensembl Annotation Tables*.
<https://github.com/stephenturner/annotables>.

</div>

<div id="ref-R-tidyverse" class="csl-entry">

Wickham, Hadley. 2021. *Tidyverse: Easily Install and Load the
Tidyverse*. <https://CRAN.R-project.org/package=tidyverse>.

</div>

<div id="ref-tidyverse2019" class="csl-entry">

Wickham, Hadley, Mara Averick, Jennifer Bryan, Winston Chang, Lucy
D’Agostino McGowan, Romain François, Garrett Grolemund, et al. 2019.
“Welcome to the <span class="nocase">tidyverse</span>.” *Journal of Open
Source Software* 4 (43): 1686. <https://doi.org/10.21105/joss.01686>.

</div>

<div id="ref-R-clusterProfiler" class="csl-entry">

Yu, Guangchuang. 2021a. *clusterProfiler: Statistical Analysis and
Visualization of Functional Profiles for Genes and Gene Clusters*.
<https://yulab-smu.top/biomedical-knowledge-mining-book/>.

</div>

<div id="ref-R-enrichplot" class="csl-entry">

———. 2021b. *Enrichplot: Visualization of Functional Enrichment Result*.
<https://yulab-smu.top/biomedical-knowledge-mining-book/>.

</div>

<div id="ref-clusterProfiler2012" class="csl-entry">

Yu, Guangchuang, Li-Gen Wang, Yanyan Han, and Qing-Yu He. 2012.
“clusterProfiler: An r Package for Comparing Biological Themes Among
Gene Clusters.” *OMICS: A Journal of Integrative Biology* 16 (5):
284–87. <https://doi.org/10.1089/omi.2011.0118>.

</div>

<div id="ref-R-heatmap3" class="csl-entry">

Zhao, Shilin, Linlin Yin, Yan Guo, Quanhu Sheng, and Yu Shyr. 2021.
*Heatmap3: An Improved Heatmap Package*.
<https://CRAN.R-project.org/package=heatmap3>.

</div>

</div>
