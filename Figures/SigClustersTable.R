##########################################################################################################
# Table 3.3.1.1 
##########################################################################################################
#Load significant modules
library(dplyr)
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingPearson.rda')
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersNonNeurologicalPearsonCoding.rda')
sigClustersCodingPearson <- rbind(sigClustersCodingPearson,sigClustersNonNeurologicalPearsonCoding)
sigClustersCodingPearson$Method <- rep('CPG',nrow(sigClustersCodingPearson))
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersAllPearson.rda')
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersNonNeurologicalPearsonAll.rda')
sigClustersAllPearson <- rbind(sigClustersAllPearson,sigClustersNonNeurologicalPearsonAll)
sigClustersAllPearson$Method <- rep('APG',nrow(sigClustersAllPearson))
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingWGCNAUnsigned.rda')
sigClustersCodingWGCNAUnsigned$Method <- rep('CWG',nrow(sigClustersCodingWGCNAUnsigned))
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersAllWGCNAUnsigned.rda')
sigClustersAllWGCNAUnsigned$Method <- rep('AWG',nrow(sigClustersAllWGCNAUnsigned))
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingPearsonTranscripts.rda')
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingNonNeurologicalPearsonTranscripts.rda')
sigClustersCodingPearsonTranscripts <- rbind(sigClustersCodingPearsonTranscripts,
                                             sigClustersCodingNonNeurologicalPearsonTranscripts)
sigClustersCodingPearsonTranscripts$Method <- rep('CPT',nrow(sigClustersCodingPearsonTranscripts))
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingWGCNATranscripts.rda')
sigClustersCodingWGCNATranscripts$Method <- rep('CWT',nrow(sigClustersCodingWGCNATranscripts))

AD_Endo <- c('AB42_Deming','ptau_Deming')
MS <- c('MS_IMSGC_2011','MS_IMSGC_2013','MS_IMSGC_2007')
VNR <- c('VNR_Davies_2018')
Memory <- c('Memory_Davies_2016')
ReactionTime <- c('ReactionTime_Davies_2018')
Education <- c('Education_Davies_2016')
Height <- c('Height_Wood')
BMD <- c('faBMD_Zheng')
StudyCategory <- list('AD Endophenotype'=AD_Endo,
                      MS=MS,VNR=VNR,Memory=Memory,
                      'Reaction Time'=ReactionTime,Education=Education,
                      'BMD'=BMD,'Height'=Height)
PhenotypeCategory <- list('Neurodegenerative'=c('AD Endophenotype','MS'),
                          'Neuronal'=c('VNR','Memory','Reaction Time','Education'),'Non-Neuronal'=c('BMD','Height'))


#Merge modules into single dataframe
clusterNames <-   c(paste0('CPG',1:nrow(sigClustersCodingPearson)),
                  paste0('APG',1:nrow(sigClustersAllPearson)),
                  paste0('CWG',1:nrow(sigClustersCodingWGCNAUnsigned)),
                  paste0('AWG',1:nrow(sigClustersAllWGCNAUnsigned)),
                  paste0('CPT',1:nrow(sigClustersCodingPearsonTranscripts)),
                  paste0('CWT',1:nrow(sigClustersCodingWGCNATranscripts)))
allMethodsDf <- rbind(sigClustersCodingPearson,
                      sigClustersAllPearson,
                      sigClustersCodingWGCNAUnsigned,
                      sigClustersAllWGCNAUnsigned,
                      sigClustersCodingPearsonTranscripts,
                      sigClustersCodingWGCNATranscripts) %>% dplyr::select(Name,Size,StudyName,Biotype,
                                                                           Method,adjPvalue,KEGG_Term=KEGG_2016.Term,
                                                                           KEGG_Overlap=KEGG_2016.Overlap,
                                                                           KEGG_P=KEGG_2016.Adjusted.P.value,
                                                                           contains('_target_P'),MicrogliaOverlap,Microglia_P,contains('_drug_overlap')) %>% dplyr::mutate(ClusterNameShort=clusterNames)
#Short version of KEGG term
allMethodsDf$KEGG_Term <- sapply(allMethodsDf$KEGG_Term,function(x) unlist(strsplit(x = x,split = '_'))[1])

#Intialize full table
collapsedDfFull <- data.frame(ClusterName=character(),ClusterNameShort=character(),P=numeric(),Size=numeric(),
                          Studies=character(),Method = character(),
                          KEGG=character(),MicrogliaEnriched=character(),
                          OpenTargetAssociation=character(),OpenTargetDrug=character())

collapsedDf <- data.frame(ClusterName=character(),ClusterNameShort=character(),P=numeric(),Size=numeric(),
                              Studies=character(),numStudies=numeric(),Method = character(),
                              KEGG=character(),MicrogliaEnriched=character(),
                              OpenTargetAssociation=character(),OpenTargetDrug=character())
uniqueClusterNames <- unique(allMethodsDf$Name)
target_P <- colnames(sigClustersCodingPearson)[grepl(colnames(sigClustersCodingPearson),pattern = '_target_P')]
drug_Overlap <- colnames(sigClustersCodingPearson)[grepl(colnames(sigClustersCodingPearson),pattern = '_drug_overlap')]

for(i in 1:length(uniqueClusterNames)){
  currentClusterName <- uniqueClusterNames[i]
  currentDf <- dplyr::filter(allMethodsDf,Name==currentClusterName)
  currentNumStudies <- length(currentDf$StudyName)
  currentStudies <- paste(currentDf$StudyName,collapse = '\n')
  currentDf <- currentDf[1,]
  currentKEGG <- paste0(currentDf$KEGG_Term,'(P=',signif(currentDf$KEGG_P,2),')')
  currentTargetAssociationFull <- paste(sapply(target_P,function(x) paste0(unlist(strsplit(x,split = '_'))[1],'(P=',signif(currentDf[,x],2),')')),collapse = '\n')
  signifTargets <- target_P[which(currentDf[,target_P]<0.05)]
  currentTargetAssociation <- paste(sapply(signifTargets,function(x) paste0(unlist(strsplit(x,split = '_'))[1],'(P=',signif(currentDf[,x],2),')')),collapse = '\n')
  
  currentMicrogliaOverlap <- paste0(currentDf$MicrogliaOverlap,'(P=',signif(currentDf$Microglia_P,2),')')
  currentDrugOverlapFull <- paste(sapply(drug_Overlap,function(x) paste0(unlist(strsplit(x,split = '_'))[1],'(P=',signif(currentDf[,x],2))),collapse = '\n')
  signifOverlap <- drug_Overlap[which(currentDf[,drug_Overlap]>0)]
  currentDrugOverlap <- paste(sapply(signifOverlap,function(x) paste0(unlist(strsplit(x,split = '_'))[1],'(',currentDf[,x],')')),collapse = '\n')
  
  collapsedDfFull <- rbind(collapsedDfFull,data.frame(ClusterName=currentClusterName,
                                                      ClusterNameShort=currentDf$ClusterNameShort,
                                                  P=signif(currentDf$adjPvalue,2),
                                                  Size=currentDf$Size,
                                                  Studies=currentStudies,
                                                  Method=currentDf$Method,
                                                  KEGG=currentKEGG,
                                                  MicrogliaEnriched=currentMicrogliaOverlap,
                                                  OpenTargetAssociation=currentTargetAssociationFull,
                                                  OpenTargetDrug=currentDrugOverlapFull))
  collapsedDf <- rbind(collapsedDf,data.frame(ClusterName=currentClusterName,
                                              ClusterNameShort=currentDf$ClusterNameShort,
                                              P=signif(currentDf$adjPvalue,2),
                                              Size=currentDf$Size,
                                              Studies=currentStudies,
                                              numStudies=currentNumStudies,
                                              Method=currentDf$Method,
                                              KEGG=currentKEGG,
                                              MicrogliaEnriched=currentMicrogliaOverlap,
                                              OpenTargetAssociation=currentTargetAssociation,
                                              OpenTargetDrug=currentDrugOverlap))
  
}
# write.csv(collapsedDfFull,file='../../FinalFigures/Supplementary/FullSigClusterTable.csv')
collapsedDf <- collapsedDf %>% dplyr::mutate(Disease=sapply(as.character(Studies),function(x) names(StudyCategory)[which(sapply(StudyCategory,function(y) unlist(strsplit(x=x,split = '\n'))[1]%in% y))]))
collapsedDf <- collapsedDf %>% dplyr::mutate('PhenotypeCategory'=sapply(Disease,function(x) names(PhenotypeCategory)[which(sapply(PhenotypeCategory,function(y) x%in% y))]))
collapsedDf <- collapsedDf %>% dplyr::arrange(PhenotypeCategory,Disease,Method,P)
collapsedDf <- collapsedDf %>% dplyr::select(Category=PhenotypeCategory,
                                             Disease=Disease,
                                             'Network Type'=Method,
                                             'Cluster Name'=ClusterNameShort,
                                             FDR=P,
                                             Size=Size,
                                             Studies=Studies,
                                             'KEGG \n Pathway'=KEGG,
                                             'Microglia \n Enriched Genes'=MicrogliaEnriched,
                                             'Open Target \n Disease Association'=OpenTargetAssociation,
                                             'Open Target \n Drug Targets'=OpenTargetDrug)
collapsedDf$Studies <- sapply(collapsedDf$Studies,function(x) gsub(pattern = '_',replacement = '.',x = x))

write.csv(collapsedDf,file='../../FinalFigures/SigClusterTable.csv',row.names = F)
library(xtable)
print.xtable(xtable(collapsedDf),include.rownames = F)

detailedDf <- rbind(sigClustersCodingPearson,
      sigClustersAllPearson,
      sigClustersCodingWGCNAUnsigned,
      sigClustersAllWGCNAUnsigned,
      sigClustersCodingPearsonTranscripts,
      sigClustersCodingWGCNATranscripts) %>% dplyr::mutate(ClusterNameShort=clusterNames)
detailedDf <- detailedDf %>% dplyr::select(SLE_drug_Genes,MS_drug_Genes,MS_target_Genes,MicrogliaGenes,Name=Name,ShortName="ClusterNameShort",KEGG_Term=KEGG_2016.Term,"GeneNames",KEGG_Genes=KEGG_2016.Genes)

microgliaSigAndKeggSig <- dplyr::filter(allMethodsDf,KEGG_2016.Adjusted.P.value < 0.05 & Microglia_P < 0.05 & StudyName %in% c(StudyCategory$`AD Endophenotype`,StudyCategory$MS))
detailedDf <- detailedDf %>% dplyr::filter(Name %in% microgliaSigAndKeggSig$Name | ShortName=='CPG12')

write.csv(detailedDf %>% dplyr::select(-Name) %>% dplyr::rename(Name=ShortName),file='../../FinalFigures/allGenes.csv',row.names = F)
