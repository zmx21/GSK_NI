source('analyze_PASCAL_output.R')
library(parallel)
# library(enrichR)
# library(biomartr)
clustersInfo <-  unlist(lapply('../../Louvain_results/microglia_gene_clusters.gmt',function(x) qusage::read.gmt(file=x)),recursive = F)
clustersDf <- data_frame(ClusterName = names(clustersInfo),Genes = clustersInfo) %>%
  dplyr::mutate(size=sapply(Genes,length)) %>% 
  dplyr::mutate(Level=sapply(ClusterName,function(x){
    splitString <- unlist(strsplit(x,'_'));
    return(as.numeric(splitString[which(splitString=='level')+1]))})) #Gene or transcript


# LevelFourClusters <- clustersDf %>% dplyr::filter(Level == 4)
# allEnrichRLib <- enrichR::listEnrichrDbs()$libraryName
# GOAnnoList <- vector(mode='list',length = nrow(LevelFourClusters))
# for(i in 1:length(GOAnnoList)){
#   print(i)
#   GOAnnoList[[i]] <- enrichR::enrichr(genes = LevelFourClusters$Genes[[i]],databases = allEnrichRLib)
# }
# save(GOAnnoList,file='../../Louvain_results/EnrichR_Microglia.rda')

library("biomaRt")
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",host = "http://grch37.ensembl.org")
biomartResults <- getBM(attributes = c("ensembl_gene_id", "go_id","name_1006","namespace_1003"),
                 filters = "ensembl_gene_id",values = unique(unlist(clustersDf$Genes,recursive = T)),
                 mart = mart) %>% filter(go_id!="")
save(biomartResults,file='../../Louvain_results/BioMart_Microglia.rda')


ExtractGOAnno <- function(genes,GOTbl){
  return(dplyr::filter(GOTbl,ensembl_gene_id %in% genes) %>% dplyr::select(GOId=go_id))
}
load('../../Louvain_results/BioMart_Microglia.rda')
load('../../Louvain_results/JoinedDfMicroglia.rda')
JoinedDfMicroglia <- JoinedDfMicroglia %>% filter(Level > 4 & Level < 8)
studies <- unique(JoinedDfMicroglia$StudyName)
topClusters <- lapply(studies,function(x) JoinedDfMicroglia %>% dplyr::filter(StudyName==x) %>% dplyr::filter(empPvalue == min(empPvalue)))
topClustersGenes <- lapply(topClusters,function(x)x  %>% dplyr::select(Genes) %>% unlist(use.names = F))
GOAnno <- lapply(topClustersGenes,function(x) ExtractGOAnno(genes = x,GOTbl = biomartResults))
for(i in 1:length(GOAnno)){
  write.table(GOAnno[[i]],
              file = paste0('../../Louvain_results/',topClusters[[i]]$StudyName,'_TopClusterGenes.txt'),row.names = F,col.names = F,quote = F)
}
