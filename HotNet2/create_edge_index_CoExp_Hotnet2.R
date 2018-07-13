##########################################################################################################
#Create edgelist, genelist, and heatlist for HotNet2 
##########################################################################################################
library(data.table)
library(dplyr)
library(hashmap)
WriteHotnetFiles <- function(edgeListPath,GWASDb,path,weightCutOff=0){
  #Only take significant edges, above weight cutoff
  edgeList <- data.table::fread(edgeListPath) %>% dplyr::filter(V3>weightCutOff) %>% select(GeneA=V1,GeneB=V2)
  #Genes must be in GWAS study
  edgeList <- edgeList %>% dplyr::filter(GeneA %in% GWASDb$gene_id & GeneB %in% GWASDb$gene_id)
  system(command = paste0('mkdir ',path)) #Make directory to write filtes
  
  #Get all unique genes
  allgenes <- union(unique(edgeList$GeneA),unique(edgeList$GeneB))
  geneIdToIndex <- hashmap::hashmap(keys = allgenes,values = seq(1,length(allgenes),1))
  
  #Get index of each gene of the ENSEMBL edge list, based on the order of allgenes.
  GeneAIndex <- geneIdToIndex[[edgeList$GeneA]]
  GeneBIndex <- geneIdToIndex[[edgeList$GeneB]]
  #Edge list based on index
  indexEdgeList <- data.frame(GeneA=GeneAIndex,GeneB=GeneBIndex)
  #Mapping between index and gene name (ENSEMBL ID)
  geneList <- data.frame(Index = 1:length(allgenes),Gene = allgenes)
  #Mapping between each gene, and it's heat score (-log10(pval))
  heatList <- dplyr::left_join(geneList,GWASDb,by=c('Gene'='gene_id')) %>% dplyr::mutate(Heat = -1*log10(pvalue)) %>% dplyr::select(Gene,Heat)
  
  #Write files to the directory.
  data.table::fwrite(indexEdgeList,quote = F,row.names = F,col.names = F,file = paste0(path,'edgelist.txt'),sep = " ")
  data.table::fwrite(geneList,quote = F,row.names = F,col.names = F,file = paste0(path,'genelist.txt'),sep = " ")
  data.table::fwrite(heatList,quote = F,row.names = F,col.names = F,file = paste0(path,'heatlist.txt'),sep = " ")
}
StudyNames <- c('AB42.Deming','ALS.vanRheenen','Alzheimer.IGAP.Stage1','MS.IMSGC.2007','Parkinsons.Pankratz','ptau.Deming','tau.Deming')
for(i in 1:length(StudyNames)){
  GWASPath <- '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_results/microglia_coding_genes/'
  GWASDb <- data.table::fread(file=paste0(GWASPath,StudyNames[i],'.sum.genescores.txt'))
  curStudyString <- gsub('[.]','_',StudyNames[i])
  WriteHotnetFiles('/local/data/public/zmx21/zmx21_private/GSK/Louvain_Edge_List/Jaccard/CodingGenesEdgeListMicroglia_Jaccard_pval0p05_cor0p25_abs.txt',
                   GWASDb,weightCutOff = 0,
                   paste0('/local/data/public/zmx21/zmx21_private/GSK/HotNet_CoExp/',curStudyString,'_coding_genes_Jaccard0/'))
  edgeList <- data.table::fread('/local/data/public/zmx21/zmx21_private/GSK/Louvain_Edge_List/Jaccard/CodingGenesEdgeListMicroglia_Jaccard_pval0p05_cor0p25_abs.txt')
  weigthCutOff <- quantile(edgeList$V3,0.25)
  WriteHotnetFiles('/local/data/public/zmx21/zmx21_private/GSK/Louvain_Edge_List/Jaccard/CodingGenesEdgeListMicroglia_Jaccard_pval0p05_cor0p25_abs.txt',
                   GWASDb,weightCutOff = weigthCutOff,
                   paste0('/local/data/public/zmx21/zmx21_private/GSK/HotNet_CoExp/',curStudyString,'_coding_genes_Jaccard25/'))
  
  # GWASPath <- '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_results/microglia_all_genes/'
  # GWASDb <- data.table::fread(file=paste0(GWASPath,StudyNames[i],'.sum.genescores.txt'))
  # WriteHotnetFiles('/local/data/public/zmx21/zmx21_private/GSK/Louvain_Edge_List/AllGenesEdgeListMicroglia.txt',
  #                  GWASDb,paste0('/local/data/public/zmx21/zmx21_private/GSK/HotNet_CoExp/',curStudyString,'_all_genes/'))
  # GWASPath <- '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_results/microglia_coding_transcripts/'
  # GWASDb <- data.table::fread(file=paste0(GWASPath,StudyNames[i],'.sum.genescores.txt'))
  # WriteHotnetFiles('/local/data/public/zmx21/zmx21_private/GSK/Louvain_Edge_List/CodingTranscriptsEdgeListMicroglia.txt',
  #                  GWASDb,paste0('/local/data/public/zmx21/zmx21_private/GSK/HotNet_CoExp/',curStudyString,'_coding_transcripts/'))
  # WriteHotnetFiles('/local/data/public/zmx21/zmx21_private/GSK/Louvain_Edge_List/AllTranscriptsEdgeListMicroglia.txt',
  #                  GWASDb,paste0('/local/data/public/zmx21/zmx21_private/GSK/HotNet_CoExp/',curStudyString,'_all_transcripts/'))
}



