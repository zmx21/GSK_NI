# load(file = '../Count_Data/geneGtfTableFull.rda')
# load(file='../Count_Data/Correlation_Matrices/Combat_Cor.rda')
TargetGenes <- c('CX3CR1','ITGAM','P2RY12','TYROBP')
TargetGenesEnsembl <- dplyr::left_join(data.frame(gene_name = TargetGenes,stringsAsFactors = F),geneGtfTableFull,by='gene_name') %>% 
  dplyr::select('gene_id') %>%
  {as.vector(t(.))}

TfDf <- read.csv(file='/local/data/public/zmx21/zmx21_private/GSK/Galatro_TF_Targets.csv') %>% 
  dplyr::select(-matches('X')) 
TfDf <- apply(TfDf,2,function(x) x[x%in%TargetGenes])

TFs <- names(TfDf)
TFsEnsembl <- dplyr::left_join(data.frame(gene_name = TFs,stringsAsFactors = F),geneGtfTableFull,by='gene_name') %>% 
  dplyr::select('gene_id') %>%
  {as.vector(t(.))}

allCorCombat <- data.frame(TF=character(),GeneTarget=character(),Cor=double())
for(i in 1:length(TFs))
{
  currentTargets <- TfDf[[i]]
  for(j in 1:length(currentTargets)){
    Cor <- tryCatch({SalmonCor_Gene_Combat[TFsEnsembl[i],
                          TargetGenesEnsembl[which(TargetGenes == currentTargets[j])]]},error = function(e){NA})
    allCorCombat <- rbind(allCorCombat,data.frame(TF = TFs[i],GeneTarget = currentTargets[j],Cor = Cor,stringsAsFactors = F))
    
  }
}

source('import_genefriends.R')
geneFriendsCorMatrix <- ImportGeneFriends(geneIDs = c(TargetGenesEnsembl,TFsEnsembl))
allCorGeneFriends <- data.frame(TF=character(),GeneTarget=character(),Cor=double())
for(i in 1:length(TFs))
{
  currentTargets <- TfDf[[i]]
  for(j in 1:length(currentTargets)){
    Cor <- tryCatch({geneFriendsCorMatrix[TFsEnsembl[i],
                                           TargetGenesEnsembl[which(TargetGenes == currentTargets[j])]]},error = function(e){NA})
    allCorGeneFriends <- rbind(allCorGeneFriends,data.frame(TF = TFs[i],GeneTarget = currentTargets[j],Cor = Cor,stringsAsFactors = F))
    
  }
}

library(ggplot2)
library(egg)
library(grid)

g1 <- ggplot(dplyr::filter(allCor,GeneTarget == TargetGenes[1]  ), aes(GeneTarget, Cor)) +   
  geom_bar(aes(fill = TF), position = "dodge", stat="identity") + ggtitle(paste0('Gene Target:',TargetGenes[1])) + scale_y_continuous(limits = c(-1,1))
g2 <- ggplot(dplyr::filter(allCor,GeneTarget == TargetGenes[2] ), aes(GeneTarget, Cor)) +   
  geom_bar(aes(fill = TF), position = "dodge", stat="identity") + ggtitle(paste0('Gene Target:',TargetGenes[2])) + scale_y_continuous(limits = c(-1,1))
g3 <- ggplot(dplyr::filter(allCor,GeneTarget == TargetGenes[3] ), aes(GeneTarget, Cor)) +   
  geom_bar(aes(fill = TF), position = "dodge", stat="identity") + ggtitle(paste0('Gene Target:',TargetGenes[3])) + scale_y_continuous(limits = c(-1,1))
g4 <- ggplot(dplyr::filter(allCor,GeneTarget == TargetGenes[4] ), aes(GeneTarget, Cor)) +   
  geom_bar(aes(fill = TF), position = "dodge", stat="identity") + ggtitle(paste0('Gene Target:',TargetGenes[4])) + scale_y_continuous(limits = c(-1,1))
ggarrange(g1,g2,g3,g4)

g5 <- ggplot(dplyr::filter(allCorGeneFriends,GeneTarget == TargetGenes[1] & !is.na(Cor)), aes(GeneTarget, Cor)) +   
  geom_bar(aes(fill = TF), position = "dodge", stat="identity") + ggtitle(paste0('Gene Target:',TargetGenes[1])) + scale_y_continuous(limits = c(-1,1))
g6 <- ggplot(dplyr::filter(allCorGeneFriends,GeneTarget == TargetGenes[2] & !is.na(Cor)), aes(GeneTarget, Cor)) +   
  geom_bar(aes(fill = TF), position = "dodge", stat="identity") + ggtitle(paste0('Gene Target:',TargetGenes[2])) + scale_y_continuous(limits = c(-1,1))
g7 <- ggplot(dplyr::filter(allCorGeneFriends,GeneTarget == TargetGenes[3] & !is.na(Cor)), aes(GeneTarget, Cor)) +   
  geom_bar(aes(fill = TF), position = "dodge", stat="identity") + ggtitle(paste0('Gene Target:',TargetGenes[3])) + scale_y_continuous(limits = c(-1,1))
g8 <- ggplot(dplyr::filter(allCorGeneFriends,GeneTarget == TargetGenes[4] & !is.na(Cor)), aes(GeneTarget, Cor)) +   
  geom_bar(aes(fill = TF), position = "dodge", stat="identity") + ggtitle(paste0('Gene Target:',TargetGenes[4])) + scale_y_continuous(limits = c(-1,1))
ggarrange(g5,g6,g7,g8)
