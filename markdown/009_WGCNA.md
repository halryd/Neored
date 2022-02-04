009\_WGCNA.R
================
xrydbh
2022-02-04

``` r
library(WGCNA)
```

    ## Loading required package: dynamicTreeCut

    ## Loading required package: fastcluster

    ## 
    ## Attaching package: 'fastcluster'

    ## The following object is masked from 'package:stats':
    ## 
    ##     hclust

    ## 
    ## Attaching package: 'WGCNA'

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     cor

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     cor

    ## The following object is masked from 'package:stats':
    ## 
    ##     cor

``` r
library(tidyverse)
library(heatmap3)
library(openxlsx)

# Set path
# If script is running separately move up one step with setwd(".."). Only needed once (for first script run separately).



############################### 
# Set parameters
############################### 
technology_platforms <- c("pea","ms")
RCutoff = .85
MCutheight = .1
PowerUpper = 20
minModuleSize = 20
softPower_ms = 5
softPower_pea = 10
exportForVisant=T
allowWGCNAThreads()
```

    ## Allowing multi-threading with up to 8 threads.

``` r
############################### 

# Get abundances
loaded_abundances <- load("./RData/003_ab_list.RData")
loaded_abundances
```

    ## [1] "ab_list"

``` r
##################### Define datExpr and allDat (names used throughout WGCNA tutorial on input dataframes)

for(technology_platform in technology_platforms){
  out_dir <- paste("./out_r/WGCNA/",technology_platform,sep="")
  dir.create(out_dir, showWarnings = F, recursive = T)
  
  tech_plf <- paste(technology_platform,"_meta",sep="")
 
  if (technology_platform=="ms"){
    datExpr <- ab_list[["dat_meta"]][[tech_plf]] %>% 
      dplyr::select(-c("UniprotID","Gene.symbol","Gene.synonyms","Numb_NA_High","Numb_NA_Low","No.peptides","PSMs","Unique.Peptides","MW.kDa")) %>%
      t() %>% 
      as_tibble() %>% 
      mutate_all(list(~as.numeric(.))) %>% 
      as.data.frame()
    
    names(datExpr) <- ab_list[["dat_meta"]][[tech_plf]] %>% 
      pluck("Gene.symbol") 
    
    allDat <- ab_list[["dat_meta"]][[tech_plf]] %>% 
      dplyr::select(-c("UniprotID","Gene.synonyms","Numb_NA_High","Numb_NA_Low","No.peptides","PSMs","Unique.Peptides","MW.kDa")) %>% 
      mutate_at(vars(-one_of(c("Gene.symbol"))),list(~as.numeric(.))) 
  } else if(technology_platform=="pea"){
    datExpr <- ab_list[["dat_meta"]][[tech_plf]] %>% 
      dplyr::select(-c("UniprotID","Gene.symbol","OlinkID","Numb_NA_High","Numb_NA_Low","Panel")) %>%
      t() %>% 
      as_tibble() %>% 
      mutate_all(list(~as.numeric(.))) %>% 
      as.data.frame()
    
    names(datExpr) <- ab_list[["dat_meta"]][[tech_plf]] %>% 
      pluck("Gene.symbol") 
    
    allDat <- ab_list[["dat_meta"]][[tech_plf]] %>% 
      dplyr::select(-c("UniprotID","OlinkID","Numb_NA_High","Numb_NA_Low","Panel")) %>% 
      mutate_at(vars(-one_of(c("Gene.symbol"))),list(~as.numeric(.))) 
    
  }
  prot_count <- allDat %>% 
    pluck("Gene.symbol") %>% 
    unique() %>% 
    length()
  
  # Create a grouping vector
  Group <- allDat %>% 
    dplyr::select(-1) %>% 
    names() %>% 
    str_remove("_") %>% str_remove("[:digit:]")
  ######################  

  #########################
  # Build WGCNA network
  #########################
  
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=PowerUpper, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, RsquaredCut = RCutoff, verbose = 5)
  
  
  # Plot the results to pdf :
  png(paste("./out_r/WGCNA/",technology_platform,"/ScaleFreeTopology_",technology_platform,".png",sep=""), 3000,2000,res=300)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  h = 0.85
  abline(h=h,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")    
  
  dev.off()
  
  
  # Get the soft power
  if(technology_platform=="ms"){
    softPower = softPower_ms# sft$powerEstimate;
  }else if(technology_platform=="pea"){
    softPower = softPower_pea# sft$powerEstimate;
  }
 
  
  # Build the adjacency table - use "signed" for proteomics data
  adjacency = adjacency(datExpr, power = softPower, type="signed");
  
  # Turn adjacency into topological overlap distance
  TOM = TOMsimilarity(adjacency)
  dissTOM = 1-TOM
  
  
  # Clustering using TOM-based dissimilarity
  
  proTree = hclust(as.dist(dissTOM), method = "average");
  
  
  # Module identification using dynamic tree cut
  dynamicMods = cutreeDynamic(dendro = proTree, distM = dissTOM, 
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize);
  
  print("Dynamic tree cut results:")
  print(table(dynamicMods))
  
  
  # Convert numeric labels into colors
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  # Plot the dendrogram and colors underneath
  sizeGrWindow(8,6)
  png(paste("./out_r/WGCNA/",technology_platform,"/DendroColor_and_colors_underneath.png",sep=""), 2000, 2000, res=300)
  plotDendroAndColors(proTree, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Protein dendrogram and module colors")
  dev.off()
  # Merge clusters
  
  mergedClust = mergeCloseModules(datExpr, dynamicColors, cutHeight = MCutheight, verbose = 3)
  
  mergedColors = mergedClust$colors;
  
  mergedMEs = mergedClust$newMEs;
  
  
  #############################################
  # Calculate and plot module eigenproteins
  # - Dendrogram
  # - Heatmap
  # - Boxplot
  # - KME
  #############################################
  
  
  # Rename to moduleColors
  moduleColors = mergedColors
  
  
  print("Modules after merging:")
  print(table(moduleColors))
  
  # Plot dendrogram
  
  png(paste("./out_r/WGCNA/",technology_platform,"/DendroColorMergedClust.png",sep=""), 2000, 2000, res=300)
  plotDendroAndColors(proTree, cbind(dynamicColors, mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  
  # Get the module eigenproteins
  MEs = mergedMEs
  
  rownames(MEs) = rownames(datExpr)
  
  
  # reorder MEs by color names of modules
  MEs = MEs[,order(colnames(MEs))]
  
  # Plot module profiles with eigenproteins overlaid
  WGCNAClusterID = moduleColors
  
  ############################
  #Plogo2 functions
  #############################
  plotErrorBarsLines <- function (v, barSizes, lines, labels = NULL, col = "blue", 
                                  ylim = c(min(lines), max(lines)), ...) 
  {
    barSizes[is.na(barSizes)] <- 0
    topBars <- v + 0.5 * barSizes
    bottomBars <- v - 0.5 * barSizes
    N <- length(v)
    if (is.null(labels)) 
      labels <- 1:N
    ylims <- c(min(bottomBars, ylim[1], min(lines)), max(topBars, 
                                                         ylim[2], max(lines)))
    par(pch = 19, xaxt = "n")
    plot(as.numeric(labels), v, ylim = ylims, col = col, type = "b", 
         lwd = 3, ...)
    par(xaxt = "s")
    
    for (i in 1:N) {
      lines(c(i, i), c(topBars[i], bottomBars[i]))
    }
    for (i in 1:ncol(lines)) {
      lines(as.numeric(labels), lines[, i], lwd = 0.5, lty = "dotted", 
            col = "gray")
    }
  }
  
  plotClusterProfileWGCNA <- function(cluster.data, moduleColors, group, MEs=NULL, 
                                      ylab="Abundance", 
                                      file="ClusterPatterns.png", ...) {
    
    gp = group
    noClusters <- nlevels(as.factor(moduleColors))
    
    r.temp <- aggregate(t(cluster.data), by=list(gp=gp), FUN=mean)
    ag.sample <- r.temp[,-1]
    rownames(ag.sample) <- r.temp[,1]
    ag.genes <- aggregate(t(ag.sample), by=list(Cluster=moduleColors), FUN=mean)
    ag.sd <- aggregate(t(ag.sample), by=list(Cluster=moduleColors), FUN=sd)
    ag.matrix <- as.matrix(ag.genes[,-1])
    
    if(!is.null(MEs) ) {        
      r.temp <- aggregate(MEs, by=list(gp=gp), FUN=mean)
      ag.matrix <- t(r.temp[,-1])
      colnames(ag.matrix) <- r.temp[,1]
    }
    ag.counts <- summary(as.factor(moduleColors))
    ag.bars <- as.matrix(ag.sd[,-1])
    
    fScale = max(8,noClusters)/8
    
    
    png(file, 2000, 3000*fScale, res=300)
    par(bg=gray(.95), fg=gray(0.3), mar= c(8, 6, 2, 1) + 0.1, col.main="black", col.sub="black", col.lab="black", col.axis="black")
    layout(matrix(1:(ceiling(noClusters/2)*2), ncol=2, byrow=TRUE))
    NSig <- noClusters
    cols = levels(as.factor(moduleColors) )
    for(i in 1:NSig) {
      gname <-  paste(levels(as.factor(moduleColors))[i], "(", ag.counts[i], "proteins )")
      lines <- ag.sample[, moduleColors==levels(as.factor(moduleColors))[i], drop=FALSE]
      plotErrorBarsLines(ag.matrix[i,], 2*ag.bars[i,], lines, 
                         labels=1:ncol(ag.matrix), 
                         col=cols[i],  main=gname, # bgcol="gray", split=split,
                         ylab=ylab, xlab="",
                         ylim=c(min(ag.matrix), max(ag.matrix)), ...)
      axis(1,at=1:ncol(ag.matrix), las=2, labels=colnames(ag.matrix), col="black", ...)
      abline(h=0, lty="dotted")
    }
    
    dev.off()
    
  }
  
  ########################
  
  # Maybe need to filter her first with : t(datExpr) %>% as_tibble() %>% # na.omit()
  t(datExpr) %>% 
    as_tibble()%>% 
    select_if(function(x) any(is.na(x))) %>% 
    summarise_each(funs(sum(is.na(.)))) -> t_datExpr_NA
  if (dim(t_datExpr_NA)[2]==0){
    #sum_NA <- sum(t_datExpr_NA[1,])
    ## plotClusterProfileWGCNA cannot handle missing values
    #if (sum_NA == 0){
    print(paste("./out_r/WGCNA/",technology_platform,"/WGCNAClusterPattenAve.png",sep=""))
    plotClusterProfileWGCNA(t(datExpr), WGCNAClusterID, Group,  MEs= MEs,
                            ylab="Average log ratio", file=paste("./out_r/WGCNA/",technology_platform,"/WGCNAClusterPattenME.png",sep=""),
                            cex.main=1.8, cex.lab=1.7, cex.axis=1.5)
    dim(t(datExpr))
    
    plotClusterProfileWGCNA(t(datExpr), WGCNAClusterID, Group,  
                            ylab="Average log ratio", file=paste("./out_r/WGCNA/",technology_platform,"/WGCNAClusterPattenAve.png",sep=""),
                            cex.main=1.8, cex.lab=1.7, cex.axis=1.5)
    # }
    
  } 
  ########################### 
  # Module trait relationships
  #####################
  Group_numeric=c(rep(1,8),rep(0,8))
  
  nGenes = ncol(datExpr);
  stemcell.conc.groups <- as.data.frame(Group_numeric)
  names(stemcell.conc.groups)="StemCellConc"
  row.names(stemcell.conc.groups) <- row.names(datExpr)
  nSamples = nrow(stemcell.conc.groups);
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr, dynamicColors)$eigengenes
  MEs1 = orderMEs(MEs0)
  #moduleTraitCor is a df with correlations of eigengenes and traits
  moduleTraitCor = WGCNA::cor(MEs1, stemcell.conc.groups, use = "p");
  #moduleTraitCor = biserial.cor(MEs, clin_out)
  #p-values of eigengene trait correlation
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  
  # Will display correlations and their p-values
  textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 10, 3, 3));
  # Display the correlation values within a heatmap plot
  
  
  pdf(paste("./out_r/WGCNA/",technology_platform,"/Eigenprotein_stemCellGroup_correlation_heatmap.pdf",sep=""), 
      width=4, 
      height=4

  )     
  par(mar = c(6, 10, 3, 3));
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(stemcell.conc.groups),
                 yLabels = names(MEs1),
                 ySymbols = names(MEs1),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.7,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
  ###############################################
  
  # dendrogram and heatmap for eigenproteins
  png(paste("./out_r/WGCNA/",technology_platform,"/Dendrogram eigenproteins.png",sep=""), 2000,2000,res=300)                    
  plotEigengeneNetworks(MEs, "Eigenprotein Network", marHeatmap = c(3,4,2,2), marDendro = c(3,4,2,5),
                        plotDendrograms = TRUE, xLabelsAngle = 90,heatmapColors=blueWhiteRed(50))   
  dev.off()
  
  
  png(paste("./out_r/WGCNA/",technology_platform,"/Heatmap eigenproteins.png",sep=""), 550,500,res=100)
  heatmap3(t(MEs), #distfun = function(x) dist(x, method="euclidean"), 
           ColSideColors=rainbow(nlevels(as.factor(Group)))[as.factor(Group)],
           method = "average", 
           main="Module eigenproteins")
  legend("topleft", fill=rainbow(nlevels(as.factor(Group)))[1:nlevels(as.factor(Group))],
         legend=levels(as.factor(Group)), cex=.6, xpd=TRUE, inset=-.1 )
  dev.off()
  
  
  # Network heatmap
  
  # Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
  plotTOM = dissTOM^7;
  # Set diagonal to NA for a nicer plot
  diag(plotTOM) = NA;
  
  png(paste("./out_r/WGCNA/",technology_platform,"/Network heatmap.png",sep=""), 2000, 2000, res=300)
  TOMplot(plotTOM, proTree, moduleColors, main = "Network heatmap plot, all proteins")
  dev.off()
  
  
  # Export networks for visulisation with VisANT - can take a while!
  if(exportForVisant) { 
    fn_list_visAnt_input_TOM <- list()
    fn_list_visAnt_input_ADJ <- list()
    for(module in moduleColors) {
      #module="brown"
      # Select module probe
      probes = names(datExpr)
      inModule = (moduleColors==module);
      modProbes = probes[inModule];
      if (length(modProbes)>2){ ## halfdan added to avoid inclusion of grey module with only 1 protein
        # Select the corresponding Topological Overlap
        
        for(NmodTOM in 1:2) {           
          # NmodTOM=1
          modTOM = TOM[inModule, inModule];
          visFile = paste("./out_r/WGCNA/",technology_platform,"/VisANTInput-TOM", module, ".txt", sep="")
          fn_list_visAnt_input_TOM[[module]] <- visFile
          
          if(NmodTOM == 2) {
            modTOM = adjacency[inModule, inModule]
            visFile = paste("./out_r/WGCNA/",technology_platform,"/VisANTInput-ADJ", module, ".txt", sep="")
            fn_list_visAnt_input_ADJ[[module]] <- visFile
          } 
          
          dimnames(modTOM) = list(modProbes, modProbes)
          # Export the network into an edge list file VisANT can read
          vis = exportNetworkToVisANT(modTOM,
                                      file = visFile,
                                      weighted = TRUE,
                                      threshold = 0
          )         
        }
      }
      
    }
    save(fn_list_visAnt_input_TOM,fn_list_visAnt_input_ADJ,file=paste("./out_r/WGCNA/",technology_platform,"/full_paths_visAnt_input.RData",sep=""))
  }
  
  # Get KME - module membership - correlation between proteins and eigenproteins
  
  kmes = signedKME(datExpr, MEs)
  # rownames(datExpr)
  # colnames(datExpr)
  
  # separate results by modules, order by kME, hub proteins on top
  
  dat.res = data.frame(allDat, moduleColors , kmes)
  
  list.cluster.dat = lapply(levels(as.factor(moduleColors)), 
                            function(x) {dtemp = dat.res[dat.res$moduleColors == x,];
                            dtemp[order(dtemp[,paste0('kME',x)==colnames(dtemp)], decreasing=TRUE),
                                  -setdiff(grep("^kME", colnames(dtemp)), which(paste0('kME',x)==colnames(dtemp)))]} )
  
  names(list.cluster.dat) =     levels(as.factor(moduleColors))
  
  
  
  # Boxplot for eigenproteins
  
  ag.temp = aggregate(MEs, by=list(Group=Group), FUN=mean)
  ag.eigengenes = t(ag.temp[,-1])
  colnames(ag.eigengenes) = ag.temp[,1]
  
  fScale = max(8,nlevels(as.factor(moduleColors)))/8
  
  
  #MEs <- MEs[,-4]
  
  png(paste("./out_r/WGCNA/",technology_platform,"/Boxplot eigenproteins.png", sep=""), 2000, 3000*fScale, res=300)
  
  par(mar= c(7, 4, 2, 1) + 0.1)
  layout(matrix(1:(ceiling(nlevels(as.factor(moduleColors))/2)*2), ncol=2, byrow=TRUE))
  cols = levels(as.factor(moduleColors))
  for(ii in 1:ncol(MEs))    
    boxplot(MEs[,ii] ~ Group, las=2, col=cols[ii], ylab = "log ratio",
            main=paste(colnames(MEs)[ii], table(moduleColors)[ii] ), cex.main=1.7, cex.lab=1.7, cex.axis=1.5 )
  
  dev.off()
  
  # Boxplot for top 6 hub proteins
  
  # list.cluster.dat[["grey"]] <- NULL
  
  
  for(ii in 1:length(list.cluster.dat)) {
    
    png(paste0("./out_r/WGCNA/",technology_platform,"/Boxplot hub proteins - ", names(list.cluster.dat)[ii], ".png"), 2000, 2500, res=300)
    par(oma= c(5, 2, 2, 1) + 0.1)
    ## Check what is the max number of genes in subcluster/module
    if (dim(list.cluster.dat[[ii]])[1]< 6){
      numb_hub_prots_to_plot <- dim(list.cluster.dat[[ii]])[1]
    } else{
      numb_hub_prots_to_plot <- 6
    }
    ##  split plotting area
    layout(matrix(1:6, ncol=2))
    for(jj in 1:numb_hub_prots_to_plot){
      x.0 <- list.cluster.dat[[ii]][jj,2:17]
      x.1 <- sapply(x.0, as.numeric)
      x.2 <- t(log(x.1)) %>% as.vector()
      #x <- sapply(x, as.numeric
      boxplot(x.2 ~ Group, 
              main=paste(list.cluster.dat[[ii]][jj,1],"\nkME=", round(list.cluster.dat[[ii]][jj,ncol(list.cluster.dat[[ii]])],2)), 
              col=rainbow(nlevels(as.factor(Group))), ylab="Log ratio", cex.main=1.5, las=2,cex.lab=1.2, )
    } 
    
    dev.off()
    
  }
  
  
  # Output results
  
  wb = createWorkbook()
  
  addWorksheet(wb, "AllData")
  
  writeData(wb, "AllData", dat.res)
  
  # write modules only tabs
  
  for(ii in 1:length(list.cluster.dat)) {
    addWorksheet(wb, names(list.cluster.dat)[ii])
    writeData(wb, names(list.cluster.dat)[ii], list.cluster.dat[[ii]])
    
  }
  
  saveWorkbook(wb, paste("./out_r/WGCNA/",technology_platform,"/ResultsWGCNA.xlsx",sep=""), overwrite=TRUE)
  
  # output the eigenprotein 
  write.csv(MEs,paste("./out_r/WGCNA/",technology_platform,"/Module eigenprotein.csv",sep=""))
  
  ## Collect used parameters
  Parameter=c("softPower","RCutoff","MCutheight","PowerUpper","minModuleSize","prot_count")
  Value=c(softPower,RCutoff,MCutheight,PowerUpper,minModuleSize,prot_count)
  used_parameters <- cbind(Parameter,Value) %>% as_tibble()
  
  # Save used parametersas Excel workbook 
  used_parameters  %>%
    writexl::write_xlsx(path = paste("./out_r/WGCNA/",technology_platform,"/Used_parameters.xlsx",sep=""))
}
```

    ## pickSoftThreshold: will use block size 417.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 417 of 417
    ##    Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
    ## 1      1   0.3110  2.150         0.9470 115.000  116.0000 157.00
    ## 2      2   0.0107  0.181         0.9360  46.800   47.0000  82.80
    ## 3      3   0.0464 -0.321         0.9320  23.200   22.9000  50.00
    ## 4      4   0.3300 -0.842         0.9520  12.900   12.5000  32.60
    ## 5      5   0.4890 -0.989         0.9690   7.850    7.3500  22.30
    ## 6      6   0.5900 -1.060         0.9380   5.070    4.4700  15.90
    ## 7      7   0.6660 -1.090         0.9340   3.430    3.0000  11.50
    ## 8      8   0.7760 -1.050         0.9730   2.420    1.9900   8.60
    ## 9      9   0.8250 -1.070         0.9470   1.750    1.3500   6.74
    ## 10    10   0.8470 -1.230         0.9680   1.310    0.9790   5.83
    ## 11    12   0.8120 -1.470         0.9310   0.775    0.5080   4.54
    ## 12    14   0.8460 -1.630         0.8980   0.492    0.2950   3.72
    ## 13    16   0.9040 -1.660         0.9120   0.329    0.1690   3.12
    ## 14    18   0.2370 -2.540         0.0201   0.229    0.0979   2.64
    ## 15    20   0.2770 -2.610         0.0993   0.165    0.0584   2.26

    ## ..connectivity..
    ## ..matrix multiplication (system BLAS)..
    ## ..normalization..
    ## ..done.
    ##  ..cutHeight not given, setting it to 0.986  ===>  99% of the (truncated) height range in dendro.
    ##  ..done.
    ## [1] "Dynamic tree cut results:"
    ## dynamicMods
    ##  1  2  3  4  5  6  7 
    ## 96 92 77 57 39 33 23

    ##  mergeCloseModules: Merging modules whose distance is less than 0.1
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 7 module eigengenes in given set.
    ##    Calculating new MEs...
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 7 module eigengenes in given set.
    ## [1] "Modules after merging:"
    ## moduleColors
    ##     black      blue     brown     green       red turquoise    yellow 
    ##        23        92        77        39        33        96        57

    ## Warning: `summarise_each_()` was deprecated in dplyr 0.7.0.
    ## Please use `across()` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.

    ## [1] "./out_r/WGCNA/pea/WGCNAClusterPattenAve.png"

    ## Warning in greenWhiteRed(50): WGCNA::greenWhiteRed: this palette is not suitable for people
    ## with green-red color blindness (the most common kind of color blindness).
    ## Consider using the function blueWhiteRed instead.

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## pickSoftThreshold: will use block size 745.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 745 of 745
    ##    Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
    ## 1      1    0.577  0.731          0.648  287.00   297.000  425.0
    ## 2      2    0.118 -0.128          0.283  155.00   154.000  297.0
    ## 3      3    0.648 -0.446          0.706   96.70    88.800  226.0
    ## 4      4    0.808 -0.656          0.857   65.70    55.500  180.0
    ## 5      5    0.862 -0.762          0.880   47.20    35.600  147.0
    ## 6      6    0.871 -0.835          0.913   35.30    24.500  124.0
    ## 7      7    0.902 -0.900          0.935   27.30    17.100  106.0
    ## 8      8    0.889 -0.979          0.929   21.60    12.200   92.7
    ## 9      9    0.852 -1.050          0.908   17.40     9.130   82.1
    ## 10    10    0.814 -1.100          0.890   14.30     6.710   73.4
    ## 11    12    0.813 -1.190          0.897   10.10     3.820   60.2
    ## 12    14    0.839 -1.220          0.929    7.38     2.260   50.5
    ## 13    16    0.809 -1.310          0.922    5.59     1.430   43.2
    ## 14    18    0.828 -1.340          0.936    4.35     0.952   37.5
    ## 15    20    0.791 -1.410          0.909    3.46     0.683   32.9

    ## ..connectivity..
    ## ..matrix multiplication (system BLAS)..
    ## ..normalization..
    ## ..done.
    ##  ..cutHeight not given, setting it to 0.95  ===>  99% of the (truncated) height range in dendro.
    ##  ..done.
    ## [1] "Dynamic tree cut results:"
    ## dynamicMods
    ##   1   2   3 
    ## 417 231  97

    ##  mergeCloseModules: Merging modules whose distance is less than 0.1
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 3 module eigengenes in given set.
    ##    Calculating new MEs...
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 3 module eigengenes in given set.
    ## [1] "Modules after merging:"
    ## moduleColors
    ##      blue     brown turquoise 
    ##       231        97       417

    ## Warning in greenWhiteRed(50): WGCNA::greenWhiteRed: this palette is not suitable for people
    ## with green-red color blindness (the most common kind of color blindness).
    ## Consider using the function blueWhiteRed instead.

    ## Warning in greenWhiteRed(50): NaNs produced

    ## Warning in greenWhiteRed(50): NaNs produced

    ## Warning in greenWhiteRed(50): NaNs produced

    ## Warning in greenWhiteRed(50): NaNs produced

    ## Warning in greenWhiteRed(50): NaNs produced

    ## Warning in greenWhiteRed(50): NaNs produced

    ## Warning in greenWhiteRed(50): NaNs produced

    ## Warning in greenWhiteRed(50): NaNs produced

    ## Warning in bplt(at[i], wid = width[i], stats = z$stats[, i], out = z$out[z$group
    ## == : Outlier (-Inf) in boxplot 1 is not drawn

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced

    ## Warning in log(x.1): NaNs produced
