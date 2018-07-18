library(dpl)
invisible(lapply(dir('../../Count_Data/PASCAL_Results/sigClusters/'),function(x) load(paste0('../../Count_Data/PASCAL_Results/sigClusters/',x))))
clusterNames <- c(paste0('CP',1:nrow(sigClustersCodingPearson)),
                  paste0('AP',1:nrow(sigClustersAllPearson)),
                  paste0('CW',1:nrow(sigClustersCodingWGCNAUnsigned)),
                  paste0('AW',1:nrow(sigClustersAllWGCNAUnsigned)))
allSigClusters <- rbind(sigClustersCodingPearson %>% dplyr::arr,
                        sigClustersAllPearson,
                        sigClustersCodingWGCNAUnsigned,
                        sigClustersAllWGCNAUnsigned)
rownames(allSigClusters) <- clusterNames


allSigClustersGenes <- allSigClusters$Genes
numSigClusters <- nrow(allSigClusters)
similiarityMatrix <- matrix(NA,nrow = numSigClusters,ncol = numSigClusters)
rownames(similiarityMatrix) <- clusterNames;colnames(similiarityMatrix) <- clusterNames
jaccardMatrix <- matrix(NA,nrow = numSigClusters,ncol = numSigClusters)
rownames(jaccardMatrix) <- clusterNames;colnames(jaccardMatrix) <- clusterNames

for(i in 1:nrow(similiarityMatrix)){
  for(j in 1:ncol(similiarityMatrix)){
    iGenes <- allSigClustersGenes[[i]]
    jGenes <- allSigClustersGenes[[j]]
    similiarityMatrix[i,j] <- length(intersect(iGenes,jGenes)) / min(c(length(iGenes),length(jGenes)))
    jaccardMatrix[i,j] <- length(intersect(iGenes,jGenes)) / length(union(iGenes,jGenes))
  }
}
library(corrplot)
corrplot(similiarityMatrix, method = "circle",is.corr = F)
simClust <- hclust(d=as.dist(1-similiarityMatrix))
plot(simClust,cex=0.8)

corrplot(jaccardMatrix, method = "circle",is.corr = F)
jaccardClust <- hclust(d=as.dist(1-jaccardMatrix))
plot(jaccardClust,cex=0.8)


#Monogenic genes
monoGenicFamilialFiles <- c('familial_ALS','familial_MS','monogenic_AD','monogenic_recessive_PD')
monoGenicFamilialGenes <- lapply(monoGenicFamilialFiles,function(x) data.table::fread(paste0('/local/data/public/zmx21/zmx21_private/GSK/Monogenic_disease/',x,'.txt'),header = F) %>% {.$V1})
names(monoGenicFamilialGenes) <- monoGenicFamilialFiles
allGeneNames <- lapply(allSigClusters$GeneNames,function(x) unlist(strsplit(split = " ",x = x)))

for(i in 1:length(monoGenicFamilialFiles)){
  overlap <- sapply(allSigClusters$GeneNames,function(x) length(intersect(monoGenicFamilialGenes[[i]],unlist(strsplit(split = " ",x = x)))))
  curDf <- data.frame(overlap)
  colnames(curDf) <- monoGenicFamilialFiles[i]
  allSigClusters <- cbind(allSigClusters,curDf)
}
