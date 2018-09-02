library(data.table)
library(dplyr)
AppendOpenTarget <- function(JoinedDfMicroglia,csvPath,codingGenesInNetwork,allGenesInNetwork){
  library(parallel)
  diseases <- c('AD','ALS','MS','OS','PD','SLE')
  diseaseAssociatationGenes <- lapply(diseases,function(x) data.table::fread(paste0(csvPath,x,'_all.csv')) %>% {.$target.gene_info.symbol})
  drugTargetGenes <- lapply(diseases,function(x) data.table::fread(paste0(csvPath,x,'_drugs.csv')) %>% {.$target.gene_info.symbol})
  for(i in 1:length(diseases)){
    print(diseases[i])
    #Get current disease gene list from open target csv file
    curDisease <- diseases[i]

    #Find intersect with respective network for hypergeometric calculation
    numCodingAssociationGenes <- length(intersect(diseaseAssociatationGenes[[i]],codingGenesInNetwork))
    numCodingTargetGenes <- length(intersect(drugTargetGenes[[i]],codingGenesInNetwork))
    numAllAssociationGenes <- length(intersect(diseaseAssociatationGenes[[i]],allGenesInNetwork))
    numAllTargetGenes <- length(intersect(drugTargetGenes[[i]],allGenesInNetwork))
    
    #Find overlap of genes, and find length of overlap
    curAssociationGenes <- mclapply(JoinedDfMicroglia$GeneNames,function(x) intersect(x,diseaseAssociatationGenes[[i]]),mc.cores = 10)
    curTargetGenes <- mclapply(JoinedDfMicroglia$GeneNames,function(x) intersect(x,drugTargetGenes[[i]]),mc.cores = 10)
    curAssociationOverlap <- sapply(curAssociationGenes,length)
    curTargetOverlap <- sapply(curTargetGenes,length)
    
    #Calculate Open Target P-values. Use different base set for coding +/- non-coding genes
    Association_P <- rep(NA,nrow(JoinedDfMicroglia))
    Target_P <- rep(NA,nrow(JoinedDfMicroglia))
    for(j in 1:length(Association_P)){
      curSize <- JoinedDfMicroglia$Size[j]
      #Hypergeometirc test for enrichment (upper tail) of genes
      if(JoinedDfMicroglia$Biotype[j] == 'coding'){
        Association_P[j] <- phyper(q=curAssociationOverlap[j]-1,
                                 m=numCodingAssociationGenes,
                                 n=length(codingGenesInNetwork) - numCodingAssociationGenes,
                                 k=curSize,lower.tail = F)
        Target_P[j] <- phyper(q=curTargetOverlap[j]-1,
                          m=numCodingTargetGenes,
                          n=length(codingGenesInNetwork) - numCodingTargetGenes,
                          k=curSize,lower.tail = F)
        
      }else{
        Association_P[j] <- phyper(q=curAssociationOverlap[j]-1,
                                 m=numAllAssociationGenes,
                                 n=length(allGenesInNetwork) - numAllAssociationGenes,
                                 k=curSize,lower.tail = F)
        Target_P[j] <- phyper(q=curTargetOverlap[j]-1,
                          m=numAllTargetGenes,
                          n=length(allGenesInNetwork) - numAllTargetGenes,
                          k=curSize,lower.tail = F)
      }
    }
    curDf <- data_frame(curAssociationOverlap,
                        Association_P,
                        curAssociationGenes,
                        curTargetOverlap,
                        Target_P,
                        curTargetGenes)
    colnames(curDf) <- c(paste0(curDisease,'_target_overlap'),
                         paste0(curDisease,'_target_P'),
                         paste0(curDisease,'_target_Genes'),
                         paste0(curDisease,'_drug_overlap'),
                         paste0(curDisease,'_drug_P'),
                         paste0(curDisease,'_drug_Genes'))
    
    JoinedDfMicroglia <- cbind(JoinedDfMicroglia,curDf)
  }
  
  return(JoinedDfMicroglia)
}