library(dplyr)
library(dendextend)
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
  if(type=='corr'){
    library(corrplot)
    corrplot(similiarityMatrix, method = "circle",is.corr = F)
  }else{
    # ## function to set label color
    # labelCol <- function(x) {
    #   if (is.leaf(x)) {
    #     ## fetch label
    #     label <- attr(x, "label") 
    #     ## set label color to red for A and B, to blue otherwise
    #     attr(x, "nodePar") <- list(lab.col=ifelse(substr(x = label,start = 2,stop = 2)=="W", "red", "blue"))
    #   }
    #   return(x)
    # }
    
    simClust <- hclust(d=as.dist(1-similiarityMatrix))
    simClust <- as.dendrogram(simClust,hang = 0.1)
    par(mar = c(5, 5, 5, 13),
        xpd = NA) # allow content to go into outer margin 
    
    plot(simClust %>% 
           set("leaves_pch", sapply(labels(simClust),function(x) ifelse(substr(x = x,start = 1,stop = 1)=="C", 16, 17))) %>%
           set("leaves_col", sapply(labels(simClust),function(x) ifelse(substr(x = x,start = 2,stop = 2)=="W", "red", "blue"))),cex=0.8)  
    legend("topright", inset = c(-0.3, 0),legend = c("WGCNA Coding","WGCNA Coding\n + lncRNA","Pearson Coding",'Pearson Coding\n + lncRNA'),
           col=c('red','red','blue','blue'), pch= c(16,17,16,17),cex = 1)
    
    }
}

###############################SIGNIFICANT CLUSTERS#####################################3
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingPearson.rda')
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersAllPearson.rda')
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingWGCNAUnsigned.rda')
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersAllWGCNAUnsigned.rda')

clusterNames <- c(paste0('CP',1:nrow(sigClustersCodingPearson)),
                  paste0('AP',1:nrow(sigClustersAllPearson)),
                  paste0('CW',1:nrow(sigClustersCodingWGCNAUnsigned)),
                  paste0('AW',1:nrow(sigClustersAllWGCNAUnsigned)))
allSigClusters <- rbind(sigClustersCodingPearson,
                        sigClustersAllPearson,
                        sigClustersCodingWGCNAUnsigned,
                        sigClustersAllWGCNAUnsigned)
rownames(allSigClusters) <- clusterNames
PlotSimilarity(allSigClusters,type='dendro')

###############################TOP TEN Percent CLUSTERS#####################################3
load('../../Count_Data/PASCAL_Results/Microglia_Pearson_cor0p2_abs.rda')
load('../../Count_Data/PASCAL_Results/JoinedDfMicrogliaWGCNAUnsigned.rda')
top20PearsonCoding <- JoinedDfMicrogliaPearson %>% dplyr::filter(Biotype=='coding') %>% dplyr::arrange(adjPvalue) %>%{.[1:20,]}
top20PearsonAll <- JoinedDfMicrogliaPearson %>% dplyr::filter(Biotype=='all') %>% dplyr::arrange(adjPvalue) %>%{.[1:20,]}

top20WGCNACoding <- JoinedDfMicrogliaWGCNAUnsigned %>% dplyr::filter(Biotype=='coding') %>% dplyr::arrange(adjPvalue) %>%{.[1:20,]}
top20WGCNAAll <- JoinedDfMicrogliaWGCNAUnsigned %>% dplyr::filter(Biotype=='all') %>% dplyr::arrange(adjPvalue) %>%{.[1:20,]}

clusterNames <- c(paste0('CP',1:nrow(top20PearsonCoding)),
                  paste0('AP',1:nrow(top20PearsonAll)),
                  paste0('CW',1:nrow(top20WGCNACoding)),
                  paste0('AW',1:nrow(top20WGCNAAll)))
top20Clusters <- rbind(top20PearsonCoding,
                        top20PearsonAll,
                        top20WGCNACoding,
                        top20WGCNAAll)
rownames(top20Clusters) <- clusterNames
PlotSimilarity(top20Clusters,type='dendro')


#Monogenic genes
# monoGenicFamilialFiles <- c('familial_ALS','familial_MS','monogenic_AD','monogenic_recessive_PD')
# monoGenicFamilialGenes <- lapply(monoGenicFamilialFiles,function(x) data.table::fread(paste0('/local/data/public/zmx21/zmx21_private/GSK/Monogenic_disease/',x,'.txt'),header = F) %>% {.$V1})
# names(monoGenicFamilialGenes) <- monoGenicFamilialFiles
# allGeneNames <- lapply(allSigClusters$GeneNames,function(x) unlist(strsplit(split = " ",x = x)))
# 
# for(i in 1:length(monoGenicFamilialFiles)){
#   overlap <- sapply(allSigClusters$GeneNames,function(x) length(intersect(monoGenicFamilialGenes[[i]],unlist(strsplit(split = " ",x = x)))))
#   curDf <- data.frame(overlap)
#   colnames(curDf) <- monoGenicFamilialFiles[i]
#   allSigClusters <- cbind(allSigClusters,curDf)
# }
