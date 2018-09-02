AppendGeneEnrichment <- function(JoinedDfMicroglia,codingGenesInNetwork,allGenesInNetwork){
  library(parallel)
  microgliaGenes <- data.table::fread('../../Count_Data/Galatro_Microglia_Core_Genes.txt',header = F)
  microgliaGenes <- unique(microgliaGenes$V2)
  
  kunkel_ad_tbl10 <- data.table::fread('../../Count_Data/AD_Kunkel_tbl10.txt',header = F)
  kunkel_ad_tbl15 <- data.table::fread('../../Count_Data/AD_Kunkel_tbl15.txt',header = F)
  AD_Genes <- unique(c(kunkel_ad_tbl10$V1,kunkel_ad_tbl15$V1))
  
  numMicrogliaCodingGenes <- length(intersect(microgliaGenes,codingGenesInNetwork))
  numMicrogliaAllGenes <- length(intersect(microgliaGenes,allGenesInNetwork))
  numADCodingGenes <- length(intersect(AD_Genes,codingGenesInNetwork))
  numADAllGenes <- length(intersect(AD_Genes,allGenesInNetwork))
  
  #Genes Names of those which are microglia core genes
  JoinedDfMicroglia$MicrogliaGenes <- lapply(JoinedDfMicroglia$GeneNames,function(x) intersect(x,microgliaGenes))
  #Number of genes which are microglia core genes
  JoinedDfMicroglia$MicrogliaOverlap <- sapply(JoinedDfMicroglia$MicrogliaGenes,length)
  #Genes Names of those which are AD genes
  JoinedDfMicroglia$ADGenes <- lapply(JoinedDfMicroglia$GeneNames,function(x) intersect(x,AD_Genes))
  #Number of genes which are AD  genes
  JoinedDfMicroglia$ADOverlap <- sapply(JoinedDfMicroglia$ADGenes,length)
  
  #P value of microglia core genes.
  Microglia_P <- rep(NA,nrow(JoinedDfMicroglia))
  AD_P <- rep(NA,nrow(JoinedDfMicroglia))
  
  for(i in 1:length(Microglia_P)){
    curMicrogliaOverlap <- JoinedDfMicroglia$MicrogliaOverlap[i]
    curADOverlap <- JoinedDfMicroglia$ADOverlap[i]
    curSize <- JoinedDfMicroglia$Size[i]
    #Hypergeometirc test for enrichment (upper tail) of genes
    if(JoinedDfMicroglia$Biotype[i] == 'coding'){
      Microglia_P[i] <- phyper(q=curMicrogliaOverlap-1,
                               m=numMicrogliaCodingGenes,
                               n=length(codingGenesInNetwork) - numMicrogliaCodingGenes,
                               k=curSize,lower.tail = F)
      AD_P[i] <- phyper(q=curADOverlap-1,
                        m=numADCodingGenes,
                        n=length(codingGenesInNetwork) - numADCodingGenes,
                        k=curSize,lower.tail = F)
      
    }else{
      Microglia_P[i] <- phyper(q=curMicrogliaOverlap-1,
                               m=numMicrogliaAllGenes,
                               n=length(allGenesInNetwork) - numMicrogliaAllGenes,
                               k=curSize,lower.tail = F)
      AD_P[i] <- phyper(q=curADOverlap-1,
                        m=numADAllGenes,
                        n=length(allGenesInNetwork) - numADAllGenes,
                        k=curSize,lower.tail = F)
    }
  }
  JoinedDfMicroglia$Microglia_P <- Microglia_P
  JoinedDfMicroglia$AD_P <- AD_P
  return(JoinedDfMicroglia)
}