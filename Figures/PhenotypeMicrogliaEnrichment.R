##########################################################################################################
# Table 3.3.3.1 
##########################################################################################################
library(dplyr)
#Load all results
load('../../Count_Data/Final_Results/allPASCALResults.rda')
#Studies within each phenotype
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
allPASCALResults <- allPASCALResults %>%
  dplyr::mutate(Disease=sapply(StudyName,function(x) names(StudyCategory)[which(sapply(StudyCategory,function(y) x%in% y))]))
PhenotypeCategory <- list('Neurodegenerative'=c('AD Endophenotype','MS','AD','PD','ALS'),
                          'Neuronal'=c('VNR','Memory','Reaction Time','Education'),'Non-Neuronal'=c('BMD','Height'))
allPASCALResults <- allPASCALResults %>%
  dplyr::mutate(PhenotypeCategory=sapply(Disease,function(x) names(PhenotypeCategory)[which(sapply(PhenotypeCategory,function(y) x%in% y))]))

#Modules of each phenotype
NeurodegenerativeClusters <- dplyr::filter(allPASCALResults,PhenotypeCategory=='Neurodegenerative')
NeuronalClusters <- dplyr::filter(allPASCALResults,PhenotypeCategory=='Neuronal')
NonNeuronalClusters <- dplyr::filter(allPASCALResults,PhenotypeCategory=='Non-Neuronal')

#Enrichment of Microglia enriched genes within each phenotype
GetEnrichment <- function(NeurodegenerativeClusters,NeuronalClusters,NonNeuronalClusters,adjPValueCutOff){
  NeurodegenerativeMicrogliaEnrichment <- phyper(nrow(dplyr::filter(NeurodegenerativeClusters,adjPvalue<adjPValueCutOff & Microglia_P<0.05))-1,
                                                 nrow(dplyr::filter(NeurodegenerativeClusters,Microglia_P<=0.05)),
                                                 nrow(dplyr::filter(NeurodegenerativeClusters,Microglia_P>0.05)),
                                                 nrow(dplyr::filter(NeurodegenerativeClusters,adjPvalue<adjPValueCutOff)),lower.tail = F)
  NeuronalMicrogliaEnrichment <- phyper(nrow(dplyr::filter(NeuronalClusters,adjPvalue<adjPValueCutOff & Microglia_P<0.05))-1,
                                        nrow(dplyr::filter(NeuronalClusters,Microglia_P<=0.05)),
                                        nrow(dplyr::filter(NeuronalClusters,Microglia_P>0.05)),
                                        nrow(dplyr::filter(NeuronalClusters,adjPvalue<adjPValueCutOff)),lower.tail = F)
  NonNeuronalMicrogliaEnrichment <- phyper(nrow(dplyr::filter(NonNeuronalClusters,adjPvalue<adjPValueCutOff & Microglia_P<0.05))-1,
                                           nrow(dplyr::filter(NonNeuronalClusters,Microglia_P<=0.05)),
                                           nrow(dplyr::filter(NonNeuronalClusters,Microglia_P>0.05)),
                                           nrow(dplyr::filter(NonNeuronalClusters,adjPvalue<adjPValueCutOff)),lower.tail = F)
  return(c(NeurodegenerativeMicrogliaEnrichment,NeuronalMicrogliaEnrichment,NonNeuronalMicrogliaEnrichment))
}

#Microglia enrichment table
enrichmentResults <- data.frame(do.call(rbind,lapply(c(0.01,0.05,0.1),function(x) GetEnrichment(NeurodegenerativeClusters,NeuronalClusters,NonNeuronalClusters,x))))
enrichmentResults <- t(enrichmentResults)
colnames(enrichmentResults) <- c('FDR < 0.01','FDR < 0.05','FDR < 0.1')

enrichmentResults <- cbind(data.frame(Category=c('Neurodegenrative','Neuronal','Non-Neuronal')),enrichmentResults)
write.csv(enrichmentResults,file='../../FinalFigures/MicrogliaEnrichmentSigClusters.csv')
