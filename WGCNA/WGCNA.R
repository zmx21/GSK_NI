library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads(nThreads = 40)

load('../../Count_Data/CV_Filtered/MicrogliaGeneCVFiltered.rda')
exprMatCoding <- as.data.frame(t(MicrogliaGeneCVFiltered$coding))
exprMatAll <- as.data.frame(t(rbind(MicrogliaGeneCVFiltered$coding,MicrogliaGeneCVFiltered$noncoding)))

##############################################SOFT THRESHOLD POWER#####################################3
PlotThresholdingIndex <- function(sft){
  # Plot the results:
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sftCoding = pickSoftThreshold(exprMatCoding, powerVector = powers, verbose = 5)
sftAll = pickSoftThreshold(exprMatAll, powerVector = powers, verbose = 5)
PlotThresholdingIndex(sftCoding)
PlotThresholdingIndex(sftAll)

#From above, we see that threshold of 6 reaches the Scale indepdence topology cutoff nearly.
softPower = 6;
adjacencyCoding = adjacency(exprMatCoding, power = softPower)
adjacencyAll= adjacency(exprMatAll, power = softPower)

dissTOMCodingSigned = 1-TOMsimilarity(adjacencyCoding,TOMType = 'signed')
dissTOMCodingUnsigned = 1-TOMsimilarity(adjacencyCoding,TOMType = 'unsigned')
dissTOMAllSigned = 1-TOMsimilarity(adjacencyAll,TOMType = 'signed')
dissTOMAllUnsigned = 1-TOMsimilarity(adjacencyAll,TOMType = 'unsigned')

# Call the hierarchical clustering function
geneTreeCodingSigned = hclust(as.dist(dissTOMCodingSigned), method = "average")
geneTreeCodingUnsigned = hclust(as.dist(dissTOMCodingUnsigned), method = "average")
geneTreeAllSigned = hclust(as.dist(dissTOMAllSigned), method = "average")
geneTreeAllUnsigned = hclust(as.dist(dissTOMAllUnsigned), method = "average")

# # Plot the resulting clustering tree (dendrogram)
# sizeGrWindow(12,9)
# plot(geneTreeCodingSigned, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
#      labels = FALSE, hang = 0.04)


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 10;
# Module identification using dynamic tree cut:
modulesCodingSigned = cutreeDynamic(dendro = geneTreeCodingSigned, distM = dissTOMCodingSigned,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)

modulesCodingUnsigned = cutreeDynamic(dendro = geneTreeCodingUnsigned, distM = dissTOMCodingUnsigned,
                                    deepSplit = 2, pamRespectsDendro = FALSE,
                                    minClusterSize = minModuleSize)

modulesAllSigned = cutreeDynamic(dendro = geneTreeAllSigned, distM = dissTOMAllSigned,
                                      deepSplit = 2, pamRespectsDendro = FALSE,
                                      minClusterSize = minModuleSize)

modulesAllUnsigned = cutreeDynamic(dendro = geneTreeAllUnsigned, distM = dissTOMAllUnsigned,
                                      deepSplit = 2, pamRespectsDendro = FALSE,
                                      minClusterSize = minModuleSize)

WriteGMT <- function(moduleVector,geneNames,path,type){
  uniqueClusters <- unique(names(moduleVector))
  outFile <- file(path,open = 'w')
  for(i in 1:length(uniqueClusters)){
    currentClusterName <- uniqueClusters[i]
    currentGenes <- geneNames[which(names(moduleVector) == currentClusterName)]
    currentGenesAsString <- paste0(currentGenes,collapse = '\t')
    strToWrite <- paste(paste0(type,'_level_1_cluster_',i),
                        length(currentGenes),
                        currentGenesAsString,
                        sep = '\t')
    write(strToWrite,file = outFile,sep = "",append = T)
  }
  close.connection(outFile)
}
WriteGMT(moduleVector = modulesCodingSigned,
         rownames(adjacencyCoding),
         path = '../../WGCNA_clusters/CodingWGCNASigned.gmt',
         type = 'coding')
WriteGMT(moduleVector = modulesCodingUnsigned,
         rownames(adjacencyCoding),
         path = '../../WGCNA_clusters/CodingWGCNAUnsigned.gmt',
         type = 'coding')
WriteGMT(moduleVector = modulesAllSigned,
         rownames(adjacencyAll),
         path = '../../WGCNA_clusters/AllWGCNASigned.gmt',
         type='all')
WriteGMT(moduleVector = modulesAllUnsigned,
         rownames(adjacencyAll),
         path = '../../WGCNA_clusters/AllWGCNAUnsigned.gmt',
         type='all')
