##############################################################################################
#Fig 3.3.2.1
##############################################################################################
library(dplyr)
library(dendextend)
#Dendrograme and heatmap of subcluster measures between significant modules.
PlotSimilarity <- function(allSigClusters,type){
  allSigClustersGenes <- allSigClusters$Genes
  numSigClusters <- nrow(allSigClusters)
  similiarityMatrix <- matrix(NA,nrow = numSigClusters,ncol = numSigClusters)
  rownames(similiarityMatrix) <- rownames(allSigClusters);colnames(similiarityMatrix) <- rownames(allSigClusters)
  jaccardMatrix <- matrix(NA,nrow = numSigClusters,ncol = numSigClusters)
  rownames(jaccardMatrix) <- rownames(allSigClusters);colnames(jaccardMatrix) <- rownames(allSigClusters)
  #For each pair of clusters, calculate overlap
  for(i in 1:nrow(similiarityMatrix)){
    for(j in 1:ncol(similiarityMatrix)){
      iGenes <- allSigClustersGenes[[i]]
      jGenes <- allSigClustersGenes[[j]]
      #Similiarity is based on length of intersect / length of smallest cluster within the pair.
      similiarityMatrix[i,j] <- length(intersect(iGenes,jGenes)) / min(c(length(iGenes),length(jGenes)))
      jaccardMatrix[i,j] <- length(intersect(iGenes,jGenes)) / length(union(iGenes,jGenes))
    }
  }
  #Heat Map
  if(type=='corr'){
    library(corrplot)
    rownames(similiarityMatrix) <- paste0(rownames(similiarityMatrix),' (',sapply(allSigClustersGenes,length),')')
    colnames(similiarityMatrix) <- paste0(colnames(similiarityMatrix),' (',sapply(allSigClustersGenes,length),')')
    
    # png(file =  '../../FinalFigures/ClustSimCorr.png',width = 600,height = 790)
    corrplot(similiarityMatrix, method = "circle",
             is.corr = F,type = 'lower',tl.srt = 45,
             tl.col = sapply(colnames(similiarityMatrix),function(x) ifelse(substr(start = 2,stop = 2,x=x)=='P','blue','red')))
    text(17, y = -5, labels = 'Subcluster Measure')
    # dev.off()

  }#Dendrogram
  else{
    simClust <- hclust(d=as.dist(1-similiarityMatrix),method = 'average')
    simClust <- as.dendrogram(simClust,hang = 0.1)
    plot(simClust %>% 
           set("leaves_pch", sapply(labels(simClust),function(x) ifelse(substr(x = x,start = 1,stop = 1)=="C", 16, 17))) %>%
           set("leaves_col", sapply(labels(simClust),function(x) ifelse(substr(x = x,start = 2,stop = 2)=="W", "red", "blue"))),cex=0.8,ylab = '1 - Subcluster Measure')  
    legend("top",inset = c(0,-0.26),legend = c("WGCNA Coding","WGCNA Coding\n + lncRNA","Pearson Coding",'Pearson Coding\n + lncRNA'),
           col=c('red','red','blue','blue'), pch= c(16,17,16,17),cex = 1,ncol = 2,y.intersp=1.4)
    
  }
}

###############################SIGNIFICANT CLUSTERS#####################################3
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingPearson.rda')
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersNonNeurologicalPearsonCoding.rda')
sigClustersCodingPearson <- rbind(sigClustersCodingPearson,sigClustersNonNeurologicalPearsonCoding)
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersAllPearson.rda')
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersNonNeurologicalPearsonAll.rda')
sigClustersAllPearson <- rbind(sigClustersAllPearson,sigClustersNonNeurologicalPearsonAll)
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingWGCNAUnsigned.rda')
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersAllWGCNAUnsigned.rda')
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingPearsonTranscripts.rda')
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingNonNeurologicalPearsonTranscripts.rda')
sigClustersCodingPearsonTranscripts <- rbind(sigClustersCodingPearsonTranscripts,
                                             sigClustersCodingNonNeurologicalPearsonTranscripts)
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingWGCNATranscripts.rda')

#Get overlapping cluster information
clusterNames <-   c(paste0('CPG',1:nrow(sigClustersCodingPearson)),
                    paste0('APG',1:nrow(sigClustersAllPearson)),
                    paste0('CWG',1:nrow(sigClustersCodingWGCNAUnsigned)),
                    paste0('AWG',1:nrow(sigClustersAllWGCNAUnsigned)),
                    paste0('CPT',1:nrow(sigClustersCodingPearsonTranscripts)),
                    paste0('CWT',1:nrow(sigClustersCodingWGCNATranscripts)))
allSigClusters <- rbind(sigClustersCodingPearson,
                      sigClustersAllPearson,
                      sigClustersCodingWGCNAUnsigned,
                      sigClustersAllWGCNAUnsigned,
                      sigClustersCodingPearsonTranscripts,
                      sigClustersCodingWGCNATranscripts) 
rownames(allSigClusters) <- clusterNames
PlotSimilarity(allSigClusters,type='dendro')
PlotSimilarity(allSigClusters,type='corr')
