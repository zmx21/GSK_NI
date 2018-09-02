#Write files of module genes, with input dataframe of significant modules
WriteResultClusters <- function(sigClusters,name){
  outpath <- paste0('../../results/',name)
  system(paste0('mkdir -p ',outpath))
  system(paste0('mkdir -p ',outpath,'/ENSEMBL'))
  system(paste0('mkdir -p ',outpath,'/HGNC'))
  
  for(i in 1:nrow(sigClusters)){
    currentOutFileEnsembl <- file(paste0(outpath,'/ENSEMBL/row',i,'.txt'))
    currentOutFileHGNC <- file(paste0(outpath,'/HGNC/row',i,'.txt'))
    
    write(sigClusters$Genes[[i]],sep = '\n',file = currentOutFileEnsembl,append = F)
    write(unlist(strsplit(sigClusters$GeneNames[[i]],split = ' ')),sep = '\n',file = currentOutFileHGNC,append = F)
    
    close.connection(currentOutFileHGNC)
    close.connection(currentOutFileEnsembl)
  }
}
WriteResultClusters(sigClustersCodingPearson %>% dplyr::arrange(KEGG_2016.Term),name='Pearson_Coding')
WriteResultClusters(sigClustersAllPearson %>% dplyr::arrange(KEGG_2016.Term),name='Pearson_Coding_lncRNA')
WriteResultClusters(sigClustersCodingWGCNAUnsigned %>% dplyr::arrange(KEGG_2016.Term),name='WGCNA_Coding')
WriteResultClusters(sigClustersAllWGCNAUnsigned %>% dplyr::arrange(KEGG_2016.Term),name='WGCNA_Coding_lncRNA')
