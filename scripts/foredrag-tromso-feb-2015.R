## -------- FUNCTION ------------
## Make plots, function:
## Pick out gene to investigate:
set.seed(1)
gene <- featureNames(eset.noNorm)[sample(1:nrow(eset.noNorm),1)] ## Randomly chosen number..
## "cl36CSOeh6ule_H8uE"

geneWiseMDSplots <- function(aeset, name, gene=gene){
  pdf(file.path(resultPath,paste0("geneWisePlot_MDS_all_3_runs_",name,".pdf")),width=15)
  par(mfrow=c(2,2));
  plotGeneWiseBoxPlot(aeset, colLabel="Case_ctrl", batchLabel="New_plate", gene=gene, legend=TRUE, main=paste(name,"Run 1-3, \n Plate (=Day), gene",gene));
  abline(h=0)
  plotGeneWiseBoxPlot(aeset, colLabel="Case_ctrl", batchLabel="Chip", gene=gene, legend=TRUE, main=paste(name,"Run 1-3, \n Chip type, gene,",gene));
  abline(h=0)
  plotGeneWiseBoxPlot(aeset, colLabel="Case_ctrl", batchLabel="Run", gene=gene, legend=TRUE, main=paste(name,"Run 1-3, \n Run, gene",gene));
  abline(h=0)
  
  ## MDS Plot:
  par(mfrow=c(2,2))
  plotMDS(aeset, colLabel="New_plate", symLabel="Case_ctrl", main=paste(name,"Run 1-3, \n Plate (=Day)"))
  ## PCA plot colored with run:
  plotMDS(aeset, colLabel="Run", symLabel="Case_ctrl", main=paste(name,"Run 1-3, \n Chip type"));
  ## Color by Chip:
  plotMDS(aeset, colLabel="Chip", symLabel="Case_ctrl", main=paste(name,"Run 1-3, \n Run"));
  dev.off()
}


## -------- END FUNCTION ------------


## Make plots: 
geneWiseMDSplots(eset.noNorm, "noNorm")
geneWiseMDSplots(d.eset.noNorm, "CCdiff_noNorm")





############################
### NORMALIZE DATA #########
############################
## Perform quantile normalization on each run separetly and merge, 
## or merge first and then normalize?

## Remember that the data eset.noNorm are already log2-transformed.
normData <- lumiN(eset.noNorm,method="quantile")
## Take difference between case and control:
d.normData<- normData[,case_labnr]
exprs(d.normData) <- exprs(normData)[,case_labnr] - exprs(normData)[,ctrl_labnr]

##############################
## Make new plots with 
## normalized data
##############################
geneWiseMDSplots(normData, "qNorm")
geneWiseMDSplots(d.normData, "CC-diff_qNorm") ## Case-control DIFFERENCE


## Density plots for all runs:
## plot density for all individuals in one plot, color by run: 
## HAve to first plot run 1, then run 3, then run 2 to see them all: 
densityAllRuns <- function(aeset,name){
  d.uu <- split(sampleNames(aeset),pData(aeset)$Run)
  order.run <- c(3,1,2)
  col.run <- c("blue","green","red")
  pdf(file.path(resultPath,paste0("density_all_3_runs_CCdiff_",name,".pdf")),width=9)
  #   par(mfrow=c(2,2))
  #   i <- 1 ## plot first density outside loop for run1, then add the remaineder densities as "lines" (see code below).
  #   #plot(density(exprs(aeset)[,d.uu[[i]][1]]),xlab="",main = paste("log2(case) - log2(ctrl),",name),lwd=.5,xlim=c(-3,3),ylim=c(0,4))
  #   for (i in order.run){
  #     plot(density(exprs(aeset)[,d.uu[[i]][1]]),xlab="",main = paste("Run",i,"\n log2(case) - log2(ctrl),",name),lwd=.5,xlim=c(-3,3),ylim=c(0,4))
  #     for (j in d.uu[[i]]){
  #       lines(density(exprs(aeset)[,j]),lwd=.5,col=col.run[i])
  #     }
  #   }
  #   
  ## All plots on top of each other:
  i <- 1 ## plot first density outside loop for run1, then add the remaineder densities as "lines" (see code below).
  plot(density(exprs(aeset)[,d.uu[[i]][1]]),xlab="",main = paste("log2(case) - log2(ctrl),",name),lwd=.5,xlim=c(-3,3),ylim=c(0,4))
  for (i in order.run){
    #plot(density(exprs(aeset)[,d.uu[[i]][1]]),xlab="",main = paste("log2(case) - log2(ctrl),",name),lwd=.5,ylim=c(0,1.3))#,xlim=c(-3,3),ylim=c(0,4))
    for (j in d.uu[[i]]){
      lines(density(exprs(aeset)[,j]),lwd=.5,col=col.run[i])
    }
  }
  
  
  dev.off()
}

densityAllRuns(d.normData,"quantile_normalized2")
densityAllRuns(d.eset.noNorm, "no_normalization")
densityAllRuns(d.combat.normData, "combat_adjusted_quantile_normalized")


densityAllRuns(eset.noNorm, "Test_ind_no_normalization")
densityAllRuns(normData, "Test_ind_quantile_normalization")

###############################
## Normalize before we merge:
###############################
normData1 <- lumiN(eset1,method="quantile")
normData2 <- lumiN(eset2,method="quantile")
normData3 <- lumiN(eset3,method="quantile")
esets.preNorm <- list(normData1, normData2, normData3)

## Merge:
preNormData <- merge(esets.preNorm)
pData(preNormData)$New_plate <- factor(pData(preNormData)$New_plate,as.character(c(1:9,10:18)))
levels(pData(preNormData)$New_plate)
## Take difference between case and control:
d.preNormData<- preNormData[,case_labnr]
exprs(d.preNormData) <- exprs(preNormData)[,case_labnr] - exprs(preNormData)[,ctrl_labnr]

## Make plots:
geneWiseMDSplots(preNormData,gene,"pre_qNorm")
geneWiseMDSplots(d.preNormData, gene,"CC-diff_pre_qNorm")

#############################
## Try Combat
############################
## Merge the three prenormalized runs and apply combat: 
batch = pData(normData)$New_plate
modcombat = model.matrix(~1, data=pData(normData))
combat_edata = ComBat(dat=exprs(normData), batch=batch, mod=modcombat)

combat.normData <- normData
exprs(combat.normData) <- combat_edata ## Add combat adjusted data.
geneWiseMDSplots(combat.normData,gene,"combat_qNorm")

## CC-Difference:
aeset <- combat.normData
tmp <- aeset[,case_labnr]
exprs(tmp) <- exprs(aeset)[,case_labnr] - exprs(aeset)[,ctrl_labnr]
d.combat.normData <- tmp

geneWiseMDSplots(d.combat.normData,gene,"CCdiff_combat_qNorm")

## Combat on the difference:

d.batch = pData(d.normData)$New_plate
d.modcombat = model.matrix(~1, data=pData(d.normData))
d.combat_edata = ComBat(dat=exprs(d.normData), batch=d.batch, mod=d.modcombat)

combat.D.normData <- d.normData
exprs(combat.D.normData) <- d.combat_edata ## Add combat adjusted data.



############################
## Mean centering
###########################
bmc.normData <- bmcFunc(normData,"New_plate")
## Finc case-control difference:
d.bmc.normData<- bmc.normData[,case_labnr]
exprs(d.bmc.normData) <- exprs(bmc.normData)[,case_labnr] - exprs(bmc.normData)[,ctrl_labnr]


######### PLOT 

## Want to plot genewise plot for one gene and four different normalization/batch adjustment methods
## in on page:

## Change the plotGeneWiseBoxPlot
plotGeneWiseBoxPlot2 <- plotGeneWiseBoxPlot
as.list (body(plotGeneWiseBoxPlot))
body(plotGeneWiseBoxPlot2)[[14]] <- substitute(min_y <- 4)
body(plotGeneWiseBoxPlot2)[[15]] <- substitute(max_y <- 9)

plotGeneWiseBoxPlot3 <- plotGeneWiseBoxPlot
as.list (body(plotGeneWiseBoxPlot))
body(plotGeneWiseBoxPlot3)[[14]] <- substitute(min_y <- -3)
body(plotGeneWiseBoxPlot3)[[15]] <- substitute(max_y <- 2)

#listEsets <- list(eset.noNorm, normData, bmc.normData, combat.normData)
set.seed(16)
genes <- c(gene, featureNames(eset.noNorm)[sample(1:nrow(eset.noNorm),3)])
name <- "All-adjust-methods-Plate-adjusted"
pdf(file.path(resultPath,paste0("geneWisePlot_PCA_all_3_runs_",name,".pdf")),width=15)
par(mfrow=c(2,2));
for (gg in genes){
  plotGeneWiseBoxPlot2(eset.noNorm, colLabel="Case_ctrl", batchLabel="New_plate", gene=gg, legend=TRUE, main=paste("No normalization \n gene",gg));
  abline(h=0)
  plotGeneWiseBoxPlot2(normData, colLabel="Case_ctrl", batchLabel="New_plate", gene=gg, legend=TRUE, main=paste("Quantile normalized"));
  abline(h=0)
  plotGeneWiseBoxPlot3(bmc.normData, colLabel="Case_ctrl", batchLabel="New_plate", gene=gg, legend=TRUE, main=paste("Quantile normalized \n Mean-centering plate batch adjusted"));
  abline(h=0)
  plotGeneWiseBoxPlot2(combat.normData, colLabel="Case_ctrl", batchLabel="New_plate", gene=gg, legend=TRUE, main=paste("Quantile normalized \n Combat plate batch adjusted"));
  abline(h=0)
}
dev.off()



## CAse-control difference:
plotGeneWiseBoxPlot4 <- plotGeneWiseBoxPlot
as.list (body(plotGeneWiseBoxPlot))
body(plotGeneWiseBoxPlot4)[[14]] <- substitute(min_y <- -3)
body(plotGeneWiseBoxPlot4)[[15]] <- substitute(max_y <- 3)

name <- "All-adjust-methods-Plate-adjusted-CCdiff"
pdf(file.path(resultPath,paste0("geneWisePlot_PCA_all_3_runs_",name,".pdf")),width=15)
par(mfrow=c(2,2));
for (gg in genes){
  plotGeneWiseBoxPlot4(d.eset.noNorm, colLabel="Case_ctrl", batchLabel="New_plate", gene=gg, legend=TRUE, main=paste("CC-difference, No normalization \n gene",gg));
  abline(h=0)
  plotGeneWiseBoxPlot4(d.normData, colLabel="Case_ctrl", batchLabel="New_plate", gene=gg, legend=TRUE, main=paste("CC-difference, Quantile normalized"));
  abline(h=0)
  plotGeneWiseBoxPlot4(d.bmc.normData, colLabel="Case_ctrl", batchLabel="New_plate", gene=gg, legend=TRUE, main=paste("CC-difference, Quantile normalized \n Mean-centering plate batch adjusted"));
  abline(h=0)
  plotGeneWiseBoxPlot4(d.combat.normData, colLabel="Case_ctrl", batchLabel="New_plate", gene=gg, legend=TRUE, main=paste("CC-difference, Quantile normalized \n Combat plate batch adjusted"));
  abline(h=0)
}
dev.off()


########## PCA PLots: 

## Run pca to investigate the first two principal components:
pcaPlots <- function(aeset, name){
  #require(ggplot)
  ex <- t(exprs(aeset))
  n <- nrow(ex)
  pcaResult <- prcomp(ex)
  pcData <- data.frame(pcaResult$x)
  ## look at how much the first few pc-components explain of the variation in the data:
  pcVar <- (pcaResult$sdev)^2 / sum(pcaResult$sdev^2) ## percent per component.
  
  y.max<- x.max <-  max(abs(c(max(pcData[1:2]), min(pcData[1:2]))))
  y.min <- x.min <- (- y.max)
  
  #pdf(file.path(resultPath,paste("pca-plot-n=",n,"all-3-runs-breastcancer",name,".pdf",sep="")))
  pl1 <- ggplot(pcData,aes(x=PC1,y=PC2)) + geom_point(aes(color=as.factor(pData(aeset)$New_plate)))+ ggtitle(paste("PCA plot of",n,"individuals. Colored by Plate. \n",name))  + xlim(x.min,x.max) + ylim(y.min,y.max) + xlab(paste("PC1 (",round(pcVar[1]*100)," %)")) + ylab(paste("PC2 (",round(pcVar[2]*100)," %)"))
  pl2 <- ggplot(pcData,aes(x=PC1,y=PC2)) + geom_point(aes(color=as.factor(pData(aeset)$Run)))+ ggtitle(paste("PCA plot of",n,"individuals.Colored by Run. \n",name))  + xlim(x.min,x.max) + ylim(y.min,y.max) + xlab(paste("PC1 (",round(pcVar[1]*100)," %)")) + ylab(paste("PC2 (",round(pcVar[2]*100)," %)"))
  pl3 <- ggplot(pcData,aes(x=PC1,y=PC2)) + geom_point(aes(color=as.factor(pData(aeset)$Chip)))+ ggtitle(paste("PCA plot of",n,"individuals. Colored by Chip. \n",name))  + xlim(x.min,x.max) + ylim(y.min,y.max) + xlab(paste("PC1 (",round(pcVar[1]*100)," %)")) + ylab(paste("PC2 (",round(pcVar[2]*100)," %)"))
  pl4 <- ggplot(pcData,aes(x=PC1,y=PC2)) + geom_point(aes(color=as.factor(pData(aeset)$User_ID)))+ ggtitle(paste("PCA plot of",n,"individuals. Colored by lab personnel. \n",name))  + xlim(x.min,x.max) + ylim(y.min,y.max) + xlab(paste("PC1 (",round(pcVar[1]*100)," %)")) + ylab(paste("PC2 (",round(pcVar[2]*100)," %)"))
  ppp <- list(pl1, pl2, pl3, pl4)
  return(ppp)
  #dev.off()
}

pl1 <- pcaPlots(eset.noNorm,"no normalization")
pl2 <- pcaPlots(normData,"quantile normalized")
pl3 <- pcaPlots(bmc.normData,"quantile normalized and batch adjusted with mean centering")
pl4 <- pcaPlots(combat.normData,"quantile normalized and batch adjusted with ComBat")

pdf(file.path(resultPath,paste("pca-plot-all-3-runs-breastcancer",name,".pdf",sep="")),width=15)
grid.arrange(pl1[[1]], pl1[[2]], pl1[[3]],pl1[[4]], ncol=2)
grid.arrange(pl2[[1]], pl2[[2]], pl2[[3]],pl2[[4]], ncol=2)
grid.arrange(pl3[[1]], pl3[[2]], pl3[[3]],pl3[[4]], ncol=2)
grid.arrange(pl4[[1]], pl4[[2]], pl4[[3]],pl4[[4]], ncol=2)
dev.off()

## Case control difference: 
pl1 <- pcaPlots(d.eset.noNorm,"no normalization")
pl2 <- pcaPlots(d.normData,"quantile normalized")
pl3 <- pcaPlots(d.bmc.normData,"quantile normalized and batch adjusted with mean centering")
pl4 <- pcaPlots(d.combat.normData,"quantile normalized and batch adjusted with ComBat")

pdf(file.path(resultPath,paste("pca-plot-all-3-runs-breastcancer-CCdiff",name,".pdf",sep="")),width=15)
grid.arrange(pl1[[1]], pl1[[2]], pl1[[3]],pl1[[4]], ncol=2)
grid.arrange(pl2[[1]], pl2[[2]], pl2[[3]],pl2[[4]], ncol=2)
grid.arrange(pl3[[1]], pl3[[2]], pl3[[3]],pl3[[4]], ncol=2)
grid.arrange(pl4[[1]], pl4[[2]], pl4[[3]],pl4[[4]], ncol=2)
dev.off()

pdf(file.path(resultPath,"pca-cc-diff-all-3-runs.pdf"),width=12)
ggplot(pcData,aes(x=PC1,y=PC2)) + geom_point(aes(color=as.factor(pData(aeset)$New_plate)))+ ggtitle(paste("PCA plot of",n,"individuals.\n Colored by Plate."))  + xlim(x.min,x.max) + ylim(y.min,y.max) + xlab(paste("PC1 (",round(pcVar[1]*100)," %)")) + ylab(paste("PC2 (",round(pcVar[2]*100)," %)"))
dev.off()