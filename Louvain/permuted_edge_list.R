library(qusage)
library(dplyr)
library(data.table)
PermutateNodeLabels <- function(trueGMT,outPath,permNum){
  clusterNames <- names(trueGMT)
  clusterLevels <- sapply(clusterNames,function(x){
    splitStr <- unlist(strsplit(x,split = '_'))
    level <- as.numeric(splitStr[grep(splitStr,pattern = 'level')+1])
    return(level)
  })
  GMTDf <- data_frame(Name=clusterNames,Level=clusterLevels,Size=sapply(trueGMT,length),Genes=trueGMT)
  uniqueLevels <- unique(clusterLevels)
  PermutedGMTDf <- data_frame(Name=character(),Level=numeric(),Genes=list())
  for(i in 1:length(uniqueLevels)){
    curLevel <- uniqueLevels[i]
    currentClusters <- dplyr::filter(GMTDf,Level==curLevel)
    currentGenes <- unlist(currentClusters$Genes)
    permutedGenes <- currentGenes[sample(1:length(currentGenes),size = length(currentGenes),replace = F)]
    
    permutedGeneList <- vector(mode = 'list',length = nrow(currentClusters))
    counter <- 1
    for(j in 1:length(permutedGeneList)){
      curSize <- currentClusters$Size[j]
      permutedGeneList[[j]] <- permutedGenes[counter:(counter+curSize-1)]
      counter <- counter + curSize
    }
    PermutedGMTDf <- rbind(PermutedGMTDf,cbind(currentClusters %>% dplyr::select(Name,Level,Size),data_frame(Genes=permutedGeneList)))
  }
  PermutedGMTDf$Genes <- lapply(PermutedGMTDf$Genes,function(x) paste0(x,collapse = '\t'))
  PermutedGMTDf$Name <- sapply(PermutedGMTDf$Name,function(x) paste0(x,'_perm_',permNum))
  fileName <- paste0('perm_',permNum,'.gmt')
  data.table::fwrite(PermutedGMTDf %>% dplyr::select(Name,Size,Genes),file = paste0(outPath,fileName),
                     row.names = F,col.names = F,sep = '\t',quote = F)
}


trueGMT <- read.gmt('../../GWAS/PASCAL_New/resources/genesets/WGCNA_size3/AllWGCNAUnsigned_Soft4_Size3.gmt')
outPath <- '../../permuted_edge_list/AllWGCNAUnsigned_Soft4_Size3/'
numPerm <- 1000
for(i in 1:numPerm){
  print(i)
  PermutateNodeLabels(trueGMT = trueGMT,outPath = outPath,permNum=i)
}