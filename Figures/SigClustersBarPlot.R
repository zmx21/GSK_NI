##########################################################################################################
# Fig 3.3.1.1 
##########################################################################################################
#Load significant modules
library(dplyr)
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingPearson.rda')
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersNonNeurologicalPearsonCoding.rda')
sigClustersCodingPearson <- rbind(sigClustersCodingPearson,sigClustersNonNeurologicalPearsonCoding)
sigClustersCodingPearson$Method <- rep('Pearson Coding Genes',nrow(sigClustersCodingPearson))
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersAllPearson.rda')
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersNonNeurologicalPearsonAll.rda')
sigClustersAllPearson <- rbind(sigClustersAllPearson,sigClustersNonNeurologicalPearsonAll)
sigClustersAllPearson$Method <- rep('Pearson Coding & lncRNA Genes',nrow(sigClustersAllPearson))
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingWGCNAUnsigned.rda')
sigClustersCodingWGCNAUnsigned$Method <- rep('WGCNA Coding Genes',nrow(sigClustersCodingWGCNAUnsigned))
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersAllWGCNAUnsigned.rda')
sigClustersAllWGCNAUnsigned$Method <- rep('WGCNA Coding & lncRNA Genes',nrow(sigClustersAllWGCNAUnsigned))
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingPearsonTranscripts.rda')
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingNonNeurologicalPearsonTranscripts.rda')
sigClustersCodingPearsonTranscripts <- rbind(sigClustersCodingPearsonTranscripts,
                                             sigClustersCodingNonNeurologicalPearsonTranscripts)
sigClustersCodingPearsonTranscripts$Method <- rep('Pearson Coding Transcripts',nrow(sigClustersCodingPearsonTranscripts))
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingWGCNATranscripts.rda')
sigClustersCodingWGCNATranscripts$Method <- rep('WGCNA Coding Transcripts',nrow(sigClustersCodingWGCNATranscripts))

#Get overlapping cluster information
allMethodsDf <- rbind(sigClustersCodingPearson,
                      sigClustersAllPearson,
                      sigClustersCodingWGCNAUnsigned,
                      sigClustersAllWGCNAUnsigned,
                      sigClustersCodingPearsonTranscripts,
                      sigClustersCodingWGCNATranscripts) %>% dplyr::select(Name,StudyName,Method,Genes,adjPvalue,Microglia_P)
clusterNames <- c(paste0('CPG',1:nrow(sigClustersCodingPearson)),
                  paste0('APG',1:nrow(sigClustersAllPearson)),
                  paste0('CWG',1:nrow(sigClustersCodingWGCNAUnsigned)),
                  paste0('AWG',1:nrow(sigClustersAllWGCNAUnsigned)),
                  paste0('CPT',1:nrow(sigClustersCodingPearsonTranscripts)),
                  paste0('CWT',1:nrow(sigClustersCodingWGCNATranscripts)))
rownames(allMethodsDf) <- clusterNames


allSigClustersGenes <- allMethodsDf$Genes
numSigClusters <- nrow(allMethodsDf)
similiarityMatrix <- matrix(NA,nrow = numSigClusters,ncol = numSigClusters)
rownames(similiarityMatrix) <- rownames(allMethodsDf);colnames(similiarityMatrix) <- rownames(allMethodsDf)
for(i in 1:nrow(similiarityMatrix)){
  for(j in 1:ncol(similiarityMatrix)){
    iGenes <- allSigClustersGenes[[i]]
    jGenes <- allSigClustersGenes[[j]]
    #Similiarity is based on length of intersect / length of smallest cluster within the pair.
    similiarityMatrix[i,j] <- length(intersect(iGenes,jGenes)) / min(c(length(iGenes),length(jGenes)))
  }
}
#Remove overlapping clusters (one which is a complete subcluster of another) if they map to same study
hc <- hclust(as.dist(1-similiarityMatrix))
overlappingClusters <- cutree(tree = hc,h=0)
adjPValues <- allMethodsDf$adjPvalue; names(adjPValues) <- rownames(allMethodsDf)

uniqueClusters <- unique(overlappingClusters)
clustersToKeep <- c()
for(i in 1:length(uniqueClusters)){
  currentGroup <- names(which(overlappingClusters==uniqueClusters[i]))
  #Remove overlaping clusters from same study
  currentGroupStudies <- allMethodsDf[currentGroup,]$StudyName
  uniqueStudies <- unique(currentGroupStudies)
  for(j in 1:length(uniqueStudies)){
    curClusterToKeep <- allMethodsDf[allMethodsDf$StudyName==uniqueStudies[j] & rownames(allMethodsDf) %in% currentGroup,]
    curClusterToKeep <- curClusterToKeep[order(curClusterToKeep$adjPvalue,decreasing = F),]
    curClusterToKeep <- rownames(curClusterToKeep)[1]
    clustersToKeep <- c(clustersToKeep,curClusterToKeep)
  }
}

#Append disease type
allMethodsDf <- allMethodsDf[clustersToKeep,] %>% dplyr::mutate(Microglia_Sig=ifelse(Microglia_P<0.04,T,F)) %>% dplyr::select(StudyName,Method,Microglia_Sig,adjPvalue)
AD <- c('Alzheimer_IGAP_Stage1','Alzheimer_Maternal_Marioni',
        'Alzheimer_Meta_Marioni','Alzheimer_Parental_Marioni',
        'Alzheimer_Paternal_Marioni','Alzheimer_PCA_Schott')
AD_Endo <- c('AB42_Deming','ptau_Deming','tau_Deming')
ALS <- 'ALS_vanRheenen'
MS <- c('MS_IMSGC_2011','MS_IMSGC_2013','MS_IMSGC_2007')
PD <- c('Parkinsons_Pankratz','Parkinson_Chang')
VNR <- c('VNR_Davies_2018','VNR_Davies_2016')
Memory <- c('Memory_Davies_2016')
ReactionTime <- c('ReactionTime_Davies_2018','ReactionTime_Davies_2016')
Education <- c('Education_Davies_2016')
Height <- c('Height_Wood')
BMD <- c('faBMD_Zheng')
StudyCategory <- list('AD Endophenotype'=AD_Endo,ALS=ALS,PD=PD,AD=AD,
                      MS=MS,VNR=VNR,Memory=Memory,
                      'Reaction Time'=ReactionTime,Education=Education,
                      'BMD'=BMD,'Height'=Height)

allMethodsDf <- allMethodsDf %>% dplyr::mutate(Disease=sapply(StudyName,function(x) names(StudyCategory)[which(sapply(StudyCategory,function(y) x%in% y))]))

countDf <- data.frame(Disease=character(),Method=character(),Count=numeric(),Type=character())
uniqueMethods <- unique(allMethodsDf$Method)
uniqueDisease <- unique(allMethodsDf$Disease)
for(i in 1:length(uniqueMethods)){
  for(j in 1:length(uniqueDisease)){
    currentCount <- dplyr::filter(allMethodsDf,Disease==uniqueDisease[j] & Method==uniqueMethods[i]) %>% nrow()
    currentCountMicroglia <- dplyr::filter(allMethodsDf,Disease==uniqueDisease[j] & Method==uniqueMethods[i] & Microglia_Sig == T) %>% nrow()
    countDf <- rbind(countDf,data.frame(Disease=uniqueDisease[j],Method=uniqueMethods[i],Count=currentCount,Type='All'))
    countDf <- rbind(countDf,data.frame(Disease=uniqueDisease[j],Method=uniqueMethods[i],Count=currentCountMicroglia,Type='Microglia Enriched'))
  }
}

#Append method information
countDf$Method <- factor(countDf$Method, levels = rev(c("Pearson Coding Genes",
                                                      'Pearson Coding & lncRNA Genes',
                                                      'Pearson Coding Transcripts',
                                                      'WGCNA Coding Genes',
                                                      'WGCNA Coding & lncRNA Genes',
                                                      'WGCNA Coding Transcripts')))
#Append phenotype information
PhenotypeCategory <- list('Neurodegenerative'=c('AD Endophenotype','MS','AD','PD','ALS'),
                          'Neuronal'=c('VNR','Memory','Reaction Time','Education'),'Non-Neuronal'=c('BMD','Height'))
countDf <- countDf %>% dplyr::mutate(PhenotypeCategory=sapply(Disease,function(x) names(PhenotypeCategory)[which(sapply(PhenotypeCategory,function(y) x%in% y))]))
allMethodsDf <- allMethodsDf %>% dplyr::mutate(PhenotypeCategory=sapply(Disease,function(x) names(PhenotypeCategory)[which(sapply(PhenotypeCategory,function(y) x%in% y))]))

library(ggplot2)
tiff(file =  '../../FinalFigures/SigClusterBarPlot.tiff',width = 1200,height = 600)

#Generate barplot
p1 <- ggplot(data = countDf, aes(x=Disease,fill=Method,y=Count)) + ylab('Number of Significant \n Trait Associated Modules') + 
  scale_fill_manual('Network Type',values=c("Pearson Coding Genes"="navy",
                                            'Pearson Coding & lncRNA Genes'="royalblue",
                                            'Pearson Coding Transcripts'="skyblue3",
                                            'WGCNA Coding Genes'='darkgreen',
                                            'WGCNA Coding & lncRNA Genes'='green4',
                                            'WGCNA Coding Transcripts'='green3'),
                    labels=c("Pearson \nCoding Genes\n",
                             'Pearson \nCoding & lncRNA Genes\n',
                             'Pearson \nCoding Transcripts\n',
                             'WGCNA \nCoding Genes\n',
                             'WGCNA \nCoding & lncRNA Genes\n',
                             'WGCNA \nCoding Transcripts\n'),
                    breaks=c("Pearson Coding Genes",
                             'Pearson Coding & lncRNA Genes',
                             'Pearson Coding Transcripts',
                             'WGCNA Coding Genes',
                             'WGCNA Coding & lncRNA Genes',
                             'WGCNA Coding Transcripts')) + 
  scale_color_manual('Network Type',values=c("Pearson Coding Genes"="navy",
                                             'Pearson Coding & lncRNA Genes'="royalblue",
                                             'Pearson Coding Transcripts'="skyblue3",
                                             'WGCNA Coding Genes'='darkgreen',
                                             'WGCNA Coding & lncRNA Genes'='green4',
                                             'WGCNA Coding Transcripts'='green3'),
                     labels=c("Pearson \nCoding Genes\n",
                              'Pearson \nCoding & lncRNA Genes\n',
                              'Pearson \nCoding Transcripts\n',
                              'WGCNA \nCoding Genes\n',
                              'WGCNA \nCoding & lncRNA Genes\n',
                              'WGCNA \nCoding Transcripts\n'),
                     breaks=c("Pearson Coding Genes",
                              'Pearson Coding & lncRNA Genes',
                              'Pearson Coding Transcripts',
                              'WGCNA Coding Genes',
                              'WGCNA Coding & lncRNA Genes',
                              'WGCNA Coding Transcripts')) + 
  xlab('Phenotype') + 
  geom_bar(stat = "identity", colour = "black", position = position_dodge(width = 0.85), width = 0.6) + coord_flip() + 
  geom_text(aes(label = ifelse(Count==1e-10,'NA',Count),colour=Method), position = position_dodge(width = 0.85), hjust = -1)
facet_names <- c(
  'Neurodegenerative' = "Neuro\ndegenerative",
  'Neuronal' = "Neuronal",
  "Non-Neuronal" = "Non-\nNeuronal",
  'All' = "All Modules",
  'Microglia Enriched' = "Microglia Enriched Modules"
)

p1 <- p1 + facet_grid(PhenotypeCategory~Type,scales = "free_y", space = "free",switch = "y",labeller = as_labeller(facet_names)) + 
  theme(strip.text.y = element_text(angle = 180,size = 14),
        strip.text.x = element_text(size = 14),
        strip.placement = 'outside',
        axis.text.x =  element_text(size=12),
        axis.text.y =  element_text(size=14),
        text = element_text(size=13)) 


library(egg)
ggsave(plot = p1,file =  '../../FinalFigures/SigClusterBarPlot.pdf',device = 'pdf',width = 14,height = 12,units = 'in')
tiff(file =  '../../FinalFigures/SigClusterBarPlot.tiff',width = 1400,height = 700)
ggarrange(p1)
dev.off()