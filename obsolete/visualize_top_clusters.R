WriteNetworkFile <- function(cluster,expMatrix){
  allGenes <- unlist(cluster$Genes)
  load(file='../../Count_Data/geneGtfTableFull.rda')
  geneIdToName <- hashmap::hashmap(keys = geneGtfTableFull$gene_id,values = geneGtfTableFull$gene_name)
  allGeneNames <- geneIdToName[[allGenes]]
  
  studyName <- cluster$StudyName
  geneScoreTbl <- data.table::fread(file=paste0('../../GWAS/PASCAL_results/KEGG_Ensembl_All/',
                                                gsub('_','.',studyName),'.sum.genescores.txt'))
  if(cluster$Biotype=='coding'){
   expMatrix <- expMatrix$coding 
  }else{
   expMatrix <- rbind(expMatrix$noncoding,expMatrix$coding)
  }
  nodeTbl <- left_join(data.frame(gene_id=allGenes),geneScoreTbl,by=c('gene_id'='gene_id')) %>% dplyr::select(gene_id,pvalue) %>% 
    dplyr::left_join(allGeneNames,by=c('gene_id'='ensembl_gene_id'))
  nodeTbl$pvalue[is.na(nodeTbl$pvalue)] <- 1
  nodeTbl$pvalue <- -1*log10(nodeTbl$pvalue)
  
  edgeInteractions <- rcorr(t(expMatrix[allGenes,]))
  upperTri <- which(upper.tri(edgeInteractions$r,diag = F),arr.ind = T)
  edgeTbl <- data.frame(Node1=character(),Node2=character(),W=numeric())
  
  for(i in 1:nrow(upperTri)){
    x <- upperTri[i,]
    edgeTbl <- rbind(edgeTbl,data.frame(Node1=allGenes[x[1]],
                                        Node2=allGenes[x[2]],
                                        W=abs(edgeInteractions$r[x[1],x[2]])))
  }
  # edgeTbl <- edgeTbl %>% filter(W>0.5)
  write.table(nodeTbl,file=paste0('../../Cytoscape_Networks/nodeTbl_',studyName,cluster$Name,'.txt'),sep = '\t',quote = F,row.names = F)
  write.table(edgeTbl,file=paste0('../../Cytoscape_Networks/edgeTbl_',studyName,cluster$Name,'.txt'),sep= '\t',quote = F,row.names = F)
  
}
load('../../Count_Data/CV_Filtered/MicrogliaGeneCVFiltered.rda')
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingPearson.rda')
WriteNetworkFile(sigClustersCodingPearson[1,],MicrogliaGeneCVFiltered)

# allStudiesMicrogliaPath <- GetPathToPASCALResults('../../GWAS/PASCAL_results/microglia_gene/')
# JoinedDfMicroglia <- ParsePASCALFile(allStudiesMicrogliaPath,'../../Louvain_results/microglia_gene_clusters.gmt')
# JoinedDfMicroglia <- AppendAdjustedPValue(JoinedDfMicroglia)
# #keep only non repeating clusters. 
# JoinedDfMicroglia <- JoinedDfMicroglia %>% dplyr::distinct(StudyName,chi2Pvalue,Biotype,.keep_all = T)
# SignificantClusters <- JoinedDfMicroglia %>% dplyr::filter(adjPval < 0.1 & Level > 3)
# load('../../Count_Data/CV_Filtered/MicrogliaGeneCVFiltered.rda')
# library(Hmisc)
# WriteNetworkFile(SignificantClusters[1,],MicrogliaGeneCVFiltered)
# WriteNetworkFile(SignificantClusters[2,],MicrogliaGeneCVFiltered)
# WriteNetworkFile(SignificantClusters[3,],MicrogliaGeneCVFiltered)
# WriteNetworkFile(SignificantClusters[4,],MicrogliaGeneCVFiltered)
# WriteNetworkFile(SignificantClusters[5,],MicrogliaGeneCVFiltered)
# WriteNetworkFile(SignificantClusters[6,],MicrogliaGeneCVFiltered)
# WriteNetworkFile(SignificantClusters[7,],MicrogliaGeneCVFiltered)
# 
# WriteNetworkFile(SignificantClusters[10,],MicrogliaGeneCVFiltered)
