library(data.table)
library(dplyr)
AppendOpenTargetEmpirical <- function(JoinedDfMicroglia,permPath,csvPath,codingGenesInNetwork,allGenesInNetwork){
  if(all(sapply(JoinedDfMicroglia$GeneNames,length)==1)){
    JoinedDfMicroglia$GeneNames <- lapply(JoinedDfMicroglia$GeneNames,function(x) unlist(strsplit(x=x,split = " ")))
  }
  diseases <- c('AD','ALS','MS','OS','PD','SLE')
  diseaseAssociatationGenes <- lapply(diseases,function(x) data.table::fread(paste0(csvPath,x,'_all.csv')) %>%
                                        dplyr::select(Genes=target.gene_info.symbol,Score=association_score.overall))
  names(diseaseAssociatationGenes) <- diseases
  drugTargetGenes <- lapply(diseases,function(x) data.table::fread(paste0(csvPath,x,'_drugs.csv')) %>% {.$target.gene_info.symbol})
  library(parallel)
  for(i in 1:length(diseases)){
    #Get current disease gene list from open target csv file
    curDisease <- diseases[i]
    print(curDisease)
    
    curOpenTargetHash <- hashmap::hashmap(keys = diseaseAssociatationGenes[[i]]$Genes,
                                          values =diseaseAssociatationGenes[[i]]$Score)
    
    #Find intersect with respective network for hypergeometric calculation
    numCodingTargetGenes <- length(intersect(drugTargetGenes[[i]],codingGenesInNetwork))
    numAllTargetGenes <- length(intersect(drugTargetGenes[[i]],allGenesInNetwork))
    
    #Find overlap of genes, and find length of overlap
    curAssociationGenes <- mclapply(JoinedDfMicroglia$GeneNames,function(x) intersect(x,diseaseAssociatationGenes[[i]]$Genes),mc.cores = 10)
    curTargetGenes <- mclapply(JoinedDfMicroglia$GeneNames,function(x) intersect(x,drugTargetGenes[[i]]),mc.cores = 10)
    curAssociationOverlap <- sapply(curAssociationGenes,length)
    curTargetOverlap <- sapply(curTargetGenes,length)
    
    Association_P <- rep(NA,nrow(JoinedDfMicroglia))
    Target_P <- rep(NA,nrow(JoinedDfMicroglia))
    
    #Load empirical score results.
    curEmpricalCodingList <- readRDS(paste0(permPath,curDisease,'_coding_perm.rds'))
    if(any(JoinedDfMicroglia$Biotype=='all')){
      curEmpricalAllList <- readRDS(paste0(permPath,curDisease,'_all_perm.rds'))
    }

    for(j in 1:length(Association_P)){
      curSize <- JoinedDfMicroglia$Size[j]
      curGeneScore <- sum(curOpenTargetHash[[JoinedDfMicroglia$GeneNames[[j]]]],na.rm = T)
      #Hypergeometirc test for enrichment (upper tail) of genes
      if(JoinedDfMicroglia$Biotype[j] == 'coding'){
        curEmpricalScore <- curEmpricalCodingList[[as.character(curSize)]]
        Association_P[j] <- sum(curEmpricalScore >= curGeneScore) / length(curEmpricalScore)
        Target_P[j] <- phyper(q=curTargetOverlap[j]-1,
                          m=numCodingTargetGenes,
                          n=length(codingGenesInNetwork) - numCodingTargetGenes,
                          k=curSize,lower.tail = F)
        
      }else{
        curEmpricalScore <- curEmpricalAllList[[as.character(curSize)]]
        Association_P[j] <- sum(curEmpricalScore >= curGeneScore) / length(curEmpricalScore)
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
# test <- AppendOpenTargetEmpirical(JoinedDfMicroglia = sigClustersCodingPearson,
#                           permPath = '../../Count_Data/OpenTarget/test/',
#                           csvPath='../../OpenTargets_scores/',
#                           codingGenesInNetwork = unique(unlist(JoinedDfMicrogliaPearson %>% dplyr::filter(Biotype=='coding') %>% {.$GeneNames})),
#                           allGenesInNetwork = unique(unlist(JoinedDfMicrogliaPearson %>% dplyr::filter(Biotype=='all') %>% {.$GeneNames})))

# test <- AppendOpenTargetEmpirical(JoinedDfMicroglia = sigClustersAllPearson,
#                           permPath = '../../Count_Data/OpenTarget/test/',
#                           csvPath='../../OpenTargets_scores/',
#                           codingGenesInNetwork = unique(unlist(JoinedDfMicrogliaPearson %>% dplyr::filter(Biotype=='coding') %>% {.$GeneNames})),
#                           allGenesInNetwork = unique(unlist(JoinedDfMicrogliaPearson %>% dplyr::filter(Biotype=='all') %>% {.$GeneNames})))