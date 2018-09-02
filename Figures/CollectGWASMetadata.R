##############################################################################################
#Table 3.3.1.2
##############################################################################################
#Load signafication modules
source('../PASCAL/write_GWAS_metadata.R')
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingPearson.rda')
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersAllPearson.rda')
PearsonGWASMetda <- GetGWASMetdata(rbind(sigClustersCodingPearson,sigClustersAllPearson)) %>% {.[-nrow(.),]}

GWASSummary <- GetGWASMetdata(rbind(sigClustersCodingPearson,sigClustersAllPearson),clusterLevel = F)
GWASSummary <- GWASSummary %>% dplyr::select(-'Num Subjects')

source('../PASCAL/write_GWAS_metadata_non_neurological.R')
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersNonNeurologicalPearsonCoding.rda')
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersNonNeurologicalPearsonAll.rda')
PearsonGWASMetdaNonNeuro <- GetGWASMetdataNonNeuro(rbind(sigClustersNonNeurologicalPearsonCoding,
                                               sigClustersNonNeurologicalPearsonAll),clusterLevel = T) %>% {.[1:2,]}
PearsonGWASMetda <- rbind(PearsonGWASMetda,PearsonGWASMetdaNonNeuro)

source('../PASCAL/write_GWAS_metadata.R')
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingWGCNAUnsigned.rda')
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersAllWGCNAUnsigned.rda')
WGCNAGWASMetda <- GetGWASMetdata(rbind(sigClustersCodingPearson,sigClustersAllPearson))  %>% dplyr::select(-contains('In Network'))

overallMetadata <- dplyr::left_join(GWASSummary,PearsonGWASMetda,by=c('Study Name'='Study Name'))
overallMetadata <- overallMetadata[,c(1,2,4,5,3,6,7)]

write.csv(overallMetadata,file='../../FinalFigures/GWASSummary.csv')
