set.seed(2)
library(dplyr)
library(data.table)
library(hashmap)
library(parallel)
#Calculate score of a random set of genes from gene list, using the provided open target hash.
GenerateSizePermutation <- function(openTargetHash,geneList,size){
  randGenes <- geneList[sample(1:length(geneList),size = size)]
  clusterScore <- sum(openTargetHash[[randGenes]],na.rm = T)
  return(clusterScore)
}

GeneratePermutations <- function(openTargetHash,geneList,sizes,disease,type){
  numPerm <- 500000
  result <- mclapply(1:length(sizes),function(i) sapply(1:numPerm,function(x) GenerateSizePermutation(openTargetHash,geneList,sizes[i])),mc.cores = 30)
  names(result) <- sizes
  saveRDS(result,file = paste0('../../Count_Data/OpenTarget/',disease,'_',type,'_perm.rds'))
}
library(hashmap)
diseases <- c('AD','ALS','MS','OS','PD','SLE')
csvPath = '../../OpenTargets_scores/'
diseaseAssociatationGenes <- lapply(diseases,function(x) data.table::fread(paste0(csvPath,x,'_all.csv')) %>%
                                                                             dplyr::select(Genes=target.gene_info.symbol,Score=association_score.overall))
load('../../Count_Data/PASCAL_Results/Microglia_Pearson_cor0p2_abs.rda')
load('../../Count_Data/PASCAL_Results/JoinedDfMicrogliaWGCNAUnsigned.rda')
codingGenesInNetwork <- unique(unlist(JoinedDfMicrogliaPearson %>% dplyr::filter(Biotype=='coding') %>% {.$GeneNames}))
allGenesInNetwork <- unique(unlist(JoinedDfMicrogliaPearson %>% dplyr::filter(Biotype=='all') %>% {.$GeneNames}))
allSizes <- unique(rbind(JoinedDfMicrogliaPearson,JoinedDfMicrogliaWGCNAUnsigned) %>% {.$Size})
names(diseaseAssociatationGenes) <- diseases
for(i in 1:length(diseaseAssociatationGenes)){
  curDisease <- names(diseaseAssociatationGenes)[i]
  #Create mapping between gene and open target score
  curOpenTargetHash <- hashmap::hashmap(keys = diseaseAssociatationGenes[[i]]$Genes,
                                        values =diseaseAssociatationGenes[[i]]$Score)
  GeneratePermutations(openTargetHash=curOpenTargetHash,
                       sizes = allSizes,
                       geneList = codingGenesInNetwork,
                       disease = curDisease,
                       type='coding')
  GeneratePermutations(openTargetHash=curOpenTargetHash,
                       sizes = allSizes,
                       geneList = allGenesInNetwork,
                       disease = curDisease,
                       type='all')
}
