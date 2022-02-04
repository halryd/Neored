006c\_reactome\_results\_in\_hierarch\_tree.R
================
xrydbh
2022-02-04

``` r
# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).



path <- c(
  "Reactome/Transport of small molecules/16 Plasma lipoprotein assembly, remodeling, and clearance (20-98)", 
  "Reactome/Transport of small molecules/16 Plasma lipoprotein assembly, remodeling, and clearance (20-98)/1 Plasma lipoprotein remodeling (13-54)", 
  "Reactome/Transport of small molecules/16 Plasma lipoprotein assembly, remodeling, and clearance (20-98)/7 HDL remodeling (9-24)",
  "Reactome/Hemostasis/3 Response to elevated platelet cytosolic Ca2+ (64-146)", 
  "Reactome/Hemostasis/3 Response to elevated platelet cytosolic Ca2+ (64-146)/2 Platelet degranulation (64-139)", 
  "Reactome/Disease/Diseases of signal transduction by growth factor receptors and second messengers/Oncogenic MAPK signaling", 
  "Reactome/Disease/Diseases of signal transduction by growth factor receptors and second messengers/Oncogenic MAPK signaling/4 Signaling by RAF1 mutants (12-49)", 
  "Reactome/Disease/Diseases of signal transduction by growth factor receptors and second messengers/Oncogenic MAPK signaling/Signaling by high-kinase activity BRAF mutants (13-44)",
  "Reactome/Disease/Diseases of signal transduction by growth factor receptors and second messengers/Oncogenic MAPK signaling/15 Signaling by RAS mutants (14-54)",
  "Reactome/Disease/Diseases of signal transduction by growth factor receptors and second messengers/Oncogenic MAPK signaling/12 Signaling downstream of RAS mutants (14-54)",
  "Reactome/Disease/Diseases of signal transduction by growth factor receptors and second messengers/Oncogenic MAPK signaling/14 Paradoxical activation of RAF signaling by kinase inactive BRAF (14-54)",
  "Reactome/Disease/Diseases of signal transduction by growth factor receptors and second messengers/Oncogenic MAPK signaling/13 Signaling by moderate kinase activity BRAF mutants (14-54)",
  "Reactome/Signal transduction/MAPK family signaling cascades/MAPK1-MAPK3 signaling/RAF-MAP kinase cascade/6 MAP2K and MAPK activation (13-49)",
  "Reactome/18 Integrin signaling (10-39)",
  "Reactome/18 Integrin signaling (10-39)/8 p130Cas linkage to MAPK signaling for integrins (9-22)",
  "Reactome/18 Integrin signaling (10-39)/17 GRB2:SOS provides linkage to MAPK signaling for Integrins (10-20)",
  "Reactome/Developmental biology",
  "Reactome/Developmental biology/9 Keratinization (8-226)",
  "Reactome/Developmental biology/9 Keratinization (8-226)/10 Formation of the cornified envelope (8-13)",
  "Reactome/Metabolism of proteins",
  "Reactome/Metabolism of proteins/11 Translation (16-339)"
)


library(data.tree); library(plyr)

png(file = "./out_r/Reactome/Significant_pathways_in_pw_hierarchy.png",   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 15) # The height of the plot in inches
(mytree <- data.tree::as.Node(data.frame(pathString = path)))
```

    ##                                                                                     levelName
    ## 1  Reactome                                                                                  
    ## 2   ¦--Transport of small molecules                                                          
    ## 3   ¦   °--16 Plasma lipoprotein assembly, remodeling, and clearance (20-98)                 
    ## 4   ¦       ¦--1 Plasma lipoprotein remodeling (13-54)                                       
    ## 5   ¦       °--7 HDL remodeling (9-24)                                                       
    ## 6   ¦--Hemostasis                                                                            
    ## 7   ¦   °--3 Response to elevated platelet cytosolic Ca2+ (64-146)                           
    ## 8   ¦       °--2 Platelet degranulation (64-139)                                             
    ## 9   ¦--Disease                                                                               
    ## 10  ¦   °--Diseases of signal transduction by growth factor receptors and second messengers  
    ## 11  ¦       °--Oncogenic MAPK signaling                                                      
    ## 12  ¦           ¦--4 Signaling by RAF1 mutants (12-49)                                       
    ## 13  ¦           ¦--Signaling by high-kinase activity BRAF mutants (13-44)                    
    ## 14  ¦           ¦--15 Signaling by RAS mutants (14-54)                                       
    ## 15  ¦           ¦--12 Signaling downstream of RAS mutants (14-54)                            
    ## 16  ¦           ¦--14 Paradoxical activation of RAF signaling by kinase inactive BRAF (14-54)
    ## 17  ¦           °--13 Signaling by moderate kinase activity BRAF mutants (14-54)             
    ## 18  ¦--Signal transduction                                                                   
    ## 19  ¦   °--MAPK family signaling cascades                                                    
    ## 20  ¦       °--MAPK1-MAPK3 signaling                                                         
    ## 21  ¦           °--RAF-MAP kinase cascade                                                    
    ## 22  ¦               °--6 MAP2K and MAPK activation (13-49)                                   
    ## 23  ¦--18 Integrin signaling (10-39)                                                         
    ## 24  ¦   ¦--8 p130Cas linkage to MAPK signaling for integrins (9-22)                          
    ## 25  ¦   °--17 GRB2:SOS provides linkage to MAPK signaling for Integrins (10-20)              
    ## 26  ¦--Developmental biology                                                                 
    ## 27  ¦   °--9 Keratinization (8-226)                                                          
    ## 28  ¦       °--10 Formation of the cornified envelope (8-13)                                 
    ## 29  °--Metabolism of proteins                                                                
    ## 30      °--11 Translation (16-339)

``` r
dev.off()
```

    ## quartz_off_screen 
    ##                 2
