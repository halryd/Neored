library(WGCNA)
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
############################### 

# Get abundances
loaded_abundances <- load("./RData/003_ab_list.RData")
loaded_abundances

##################### Define datExpr and allDat (names used throughout WGCNA tutorial on input dataframes)

for(technology_platform in technology_platforms){
  out_dir <- paste("./out_r/WGCNA/",technology_platform,"_unip",sep="")
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
      pluck("UniprotID") 
    
    allDat <- ab_list[["dat_meta"]][[tech_plf]] %>% 
      dplyr::select(-c("Gene.symbol","Gene.synonyms","Numb_NA_High","Numb_NA_Low","No.peptides","PSMs","Unique.Peptides","MW.kDa")) %>% 
      mutate_at(vars(-one_of(c("UniprotID"))),list(~as.numeric(.))) 
  } else if(technology_platform=="pea"){
    datExpr <- ab_list[["dat_meta"]][[tech_plf]] %>% 
      dplyr::select(-c("UniprotID","Gene.symbol","OlinkID","Numb_NA_High","Numb_NA_Low","Panel")) %>%
      t() %>% 
      as_tibble() %>% 
      mutate_all(list(~as.numeric(.))) %>% 
      as.data.frame()
    
    names(datExpr) <- ab_list[["dat_meta"]][[tech_plf]] %>% 
      pluck("UniprotID") 
    
    allDat <- ab_list[["dat_meta"]][[tech_plf]] %>% 
      dplyr::select(-c("Gene.symbol","OlinkID","Numb_NA_High","Numb_NA_Low","Panel")) %>% 
      mutate_at(vars(-one_of(c("UniprotID"))),list(~as.numeric(.))) 
    
  }
  prot_count <- allDat %>% 
    pluck("UniprotID") %>% 
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
  png(paste("./out_r/WGCNA/",technology_platform,"_unip/ScaleFreeTopology_",technology_platform,".png",sep=""), 3000,2000,res=300)
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
  png(paste("./out_r/WGCNA/",technology_platform,"_unip/DendroColor_and_colors_underneath.png",sep=""), 2000, 2000, res=300)
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
  
  png(paste("./out_r/WGCNA/",technology_platform,"_unip/DendroColorMergedClust.png",sep=""), 2000, 2000, res=300)
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
    print(paste("./out_r/WGCNA/",technology_platform,"_unip/WGCNAClusterPattenAve.png",sep=""))
    plotClusterProfileWGCNA(t(datExpr), WGCNAClusterID, Group,  MEs= MEs,
                            ylab="Average log ratio", file=paste("./out_r/WGCNA/",technology_platform,"_unip/WGCNAClusterPattenME.png",sep=""),
                            cex.main=1.8, cex.lab=1.7, cex.axis=1.5)
    dim(t(datExpr))
    
    plotClusterProfileWGCNA(t(datExpr), WGCNAClusterID, Group,  
                            ylab="Average log ratio", file=paste("./out_r/WGCNA/",technology_platform,"_unip/WGCNAClusterPattenAve.png",sep=""),
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
  moduleTraitCor = cor(MEs1, stemcell.conc.groups, use = "p");
  #moduleTraitCor = biserial.cor(MEs, clin_out)
  #p-values of eigengene trait correlation
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  
  # Will display correlations and their p-values
  textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 10, 3, 3));
  # Display the correlation values within a heatmap plot
  
  
  pdf(paste("./out_r/WGCNA/",technology_platform,"_unip/Eigenprotein_stemCellGroup_correlation_heatmap.pdf",sep=""), 
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
  png(paste("./out_r/WGCNA/",technology_platform,"_unip/Dendrogram eigenproteins.png",sep=""), 2000,2000,res=300)					
  plotEigengeneNetworks(MEs, "Eigenprotein Network", marHeatmap = c(3,4,2,2), marDendro = c(3,4,2,5),
                        plotDendrograms = TRUE, xLabelsAngle = 90,heatmapColors=blueWhiteRed(50))	
  dev.off()
  
  
  png(paste("./out_r/WGCNA/",technology_platform,"_unip/Heatmap eigenproteins.png",sep=""), 550,500,res=100)
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
  
  png(paste("./out_r/WGCNA/",technology_platform,"_unip/Network heatmap.png",sep=""), 2000, 2000, res=300)
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
          visFile = paste("./out_r/WGCNA/",technology_platform,"_unip/VisANTInput-TOM", module, ".txt", sep="")
          fn_list_visAnt_input_TOM[[module]] <- visFile
          
          if(NmodTOM == 2) {
            modTOM = adjacency[inModule, inModule]
            visFile = paste("./out_r/WGCNA/",technology_platform,"_unip/VisANTInput-ADJ", module, ".txt", sep="")
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
    save(fn_list_visAnt_input_TOM,fn_list_visAnt_input_ADJ,file=paste("./out_r/WGCNA/",technology_platform,"_unip/full_paths_visAnt_input.RData",sep=""))
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
  
  names(list.cluster.dat) = 	levels(as.factor(moduleColors))
  
  
  
  # Boxplot for eigenproteins
  
  ag.temp = aggregate(MEs, by=list(Group=Group), FUN=mean)
  ag.eigengenes = t(ag.temp[,-1])
  colnames(ag.eigengenes) = ag.temp[,1]
  
  fScale = max(8,nlevels(as.factor(moduleColors)))/8
  
  
  #MEs <- MEs[,-4]
  
  png(paste("./out_r/WGCNA/",technology_platform,"_unip/Boxplot eigenproteins.png", sep=""), 2000, 3000*fScale, res=300)
  
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
    
    png(paste0("./out_r/WGCNA/",technology_platform,"_unip/Boxplot hub proteins - ", names(list.cluster.dat)[ii], ".png"), 2000, 2500, res=300)
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
  
  saveWorkbook(wb, paste("./out_r/WGCNA/",technology_platform,"_unip/ResultsWGCNA.xlsx",sep=""), overwrite=TRUE)
  
  # output the eigenprotein	
  write.csv(MEs,paste("./out_r/WGCNA/",technology_platform,"_unip/Module eigenprotein.csv",sep=""))
  
  ## Collect used parameters
  Parameter=c("softPower","RCutoff","MCutheight","PowerUpper","minModuleSize","prot_count")
  Value=c(softPower,RCutoff,MCutheight,PowerUpper,minModuleSize,prot_count)
  used_parameters <- cbind(Parameter,Value) %>% as_tibble()
  
  # Save used parametersas Excel workbook 
  used_parameters  %>%
    writexl::write_xlsx(path = paste("./out_r/WGCNA/",technology_platform,"_unip/Used_parameters.xlsx",sep=""))
}
