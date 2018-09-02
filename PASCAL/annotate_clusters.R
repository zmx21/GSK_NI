#Use EnrichR to find KEGG annotation for each module.
library(gridExtra)
#Called by other functions to call enrichR for a vector of genes. 
GetSingleClusterAnnotation <- function(genes,librariesToRun){
  #Call enrichR
  invisible(capture.output(df <- enrichR::enrichr(genes = genes,databases = librariesToRun)))
  #Collapse all reuslts into data frame
  df <- do.call(rbind,df)
  #Shorten library name
  df$Lib <- sapply(rownames(df),function(x){
    splitStr <- strsplit(x,'[.]')
    unlist(splitStr)[1]
  })
  return(df)
}
#Get top result of enrichR results, for each library.
GetTopTermsForEachLib <- function(enrichrResult,librariesToRun){
  allDf <- vector(mode='list',length = length(librariesToRun))
  #If enrichR did not return results, return single row data frame. 
  if(nrow(enrichrResult) == 0){
    for(i in 1:length(librariesToRun)){
      lib <- librariesToRun[i]
      curDf <- data.frame(Term=NA,Adjusted.P.value=1,Overlap=0,Genes=NA) 
      colnames(curDf) <- sapply(colnames(curDf),function(x)paste0(lib,'.',x))
      allDf[[i]] <- curDf
    }
    return(do.call(cbind,allDf))
  #If enrichR returned results, then resturn actual results.
  }else{
    for(i in 1:length(librariesToRun)){
      lib <- librariesToRun[i]
      if(!'Adjusted.P.value' %in% colnames(enrichrResult %>% dplyr::filter(Lib==lib))){
        return(do.call(cbind,allDf))
      }
      #Sort by pvalue, select appropriate columns
      curDf <- enrichrResult %>% dplyr::filter(Lib==lib) %>% 
        dplyr::arrange(Adjusted.P.value) %>% dplyr::select(Term,Adjusted.P.value,Overlap,Genes) %>% 
        {.[1,]}
      colnames(curDf) <- sapply(colnames(curDf),function(x)paste0(lib,'.',x))
      allDf[[i]] <- curDf
    }
    return(do.call(cbind,allDf))
  }
}
#Get annotation of multiple cluster (From a dataframe)
#Do so for the modules specified in topPercentile
GetClustersAnnotation <- function(JoinedDfMicroglia,codingGenes,geneIdToName,librariesToRun,topPercentile=0.1,save = F,outprefix=''){
  uniqueStudies <- unique(JoinedDfMicroglia$StudyName)
  clusterResults <- vector(mode = 'list',length = length(uniqueStudies))
  names(clusterResults) <- uniqueStudies
  for(i in 1:length(uniqueStudies)){
    print(uniqueStudies[i])
    #Get all clusters, coding and all, as level
    currentStudyDf <- JoinedDfMicroglia %>% dplyr::filter(StudyName==uniqueStudies[i])
    #Arrange by increasing p value
    currentStudyDf <- currentStudyDf %>% dplyr::arrange(adjPvalue)
    #If raw number rather than percentile specified.
    if(topPercentile > 1){
      topClusters <- currentStudyDf[1:topPercentile,]
    }else{
      topClusters <- currentStudyDf[1:floor(topPercentile*nrow(currentStudyDf)),]
    }
    allGenes <- lapply(topClusters$Genes,function(x) geneIdToName[[x]])
    topClusters$GeneNames <- sapply(allGenes,function(x) paste(x,collapse = ' '))
    #Run enrichR for each cluster
    enrichrResults <- vector(mode='list',length = length(allGenes))
    pb <- txtProgressBar(min=0,max=length(enrichrResults))
    for(j in 1:length(enrichrResults)){
      setTxtProgressBar(pb, j)
      currentGenes <- allGenes[[j]]
      currentGenes <- currentGenes[which(codingGenes[[currentGenes]])]
      enrichrResults[[j]] <- tryCatch({GetSingleClusterAnnotation(genes = allGenes[[j]],librariesToRun = librariesToRun)},
                                      warning = function(w) {
                                      }, error = function(e) {
                                        enrichrResults[[j-1]][0,]
                                      }, finally = {
                                      })
    }
    close(pb)
    names(enrichrResults) <- topClusters$Name
    #Keep the top ranking term, p value, overlap and genes for the cluster.
    topForEachCluster <- do.call(rbind.fill,lapply(enrichrResults,function(x) GetTopTermsForEachLib(x,librariesToRun)))
    currentStudyResult <- list(df=cbind(topClusters,topForEachCluster),enrichrResult=enrichrResults)
    if(save){
      save(currentStudyResult,file=paste0('../../Count_Data/PASCAL_Results/Annotations/ClusterAnnot_',outprefix,'_',uniqueStudies[i],'.rda'))
    }
    clusterResults[[i]] <- currentStudyResult
  }
  return(clusterResults)
}

#Get annotation of multiple cluster (From a dataframe)
#Do so for all modules
GetAllClustersAnnotation <- function(JoinedDfMicroglia,codingGenes,geneIdToName,librariesToRun,save = F,outprefix=''){
    allGenes <- lapply(JoinedDfMicroglia$Genes,function(x) geneIdToName[[x]])
    JoinedDfMicroglia$GeneNames <- sapply(allGenes,function(x) paste(x,collapse = ' '))
    #Run enrichR for each cluster
    enrichrResults <- vector(mode='list',length = length(allGenes))
    print(length(enrichrResults))
    pb <- txtProgressBar(min=0,max=length(enrichrResults))
    for(j in 1:length(enrichrResults)){
      setTxtProgressBar(pb, j)
      currentGenes <- allGenes[[j]]
      currentGenes <- currentGenes[which(codingGenes[[currentGenes]])]
      enrichrResults[[j]] <- tryCatch({GetSingleClusterAnnotation(genes = allGenes[[j]],librariesToRun = librariesToRun)},
                                      warning = function(w) {
                                      }, error = function(e) {
                                        print('error')
                                        enrichrResults[[j-1]]
                                      }, finally = {
                                      })
    }
    close(pb)
    names(enrichrResults) <- JoinedDfMicroglia$Name
    #Keep the top ranking term, p value, overlap and genes for the cluster.
    topForEachCluster <- do.call(rbind.fill,lapply(enrichrResults,function(x) GetTopTermsForEachLib(x,librariesToRun)))
    topForEachCluster$Name <- names(enrichrResults)
  return(topForEachCluster)
}

library(enrichR)
library(hashmap)
library(plyr)
library(dplyr)
source('../BatchCorrection_and_QC/load_GTF.R')
load(file='../../Count_Data/geneGtfTableFull.rda')
#Obtain coding genes
geneIdToName <- hashmap::hashmap(keys = geneGtfTableFull$gene_id,values = geneGtfTableFull$gene_name)
codingGenes <- hashmap::hashmap(keys = geneGtfTableFull$gene_name,values = ifelse(geneGtfTableFull$gene_biotype=='protein_coding',T,F))
# librariesToRun <- c('KEGG_2016','TRANSFAC_and_JASPAR_PWMs','Jensen_DISEASES')
librariesToRun <- c('KEGG_2016')

#Main Function, For each study, obtain annotations of modules with top disease associations.
GetTopPercentileOfEachStudy <- function(){
  method <- c('JoinedDfMicrogliaWGCNATranscripts','JoinedDfMicrogliaPearsonTranscripts')
  for(i in 1:length(method)){
    JoinedDfMicroglia <- readRDS(paste0('../../Count_Data/PASCAL_Results/Empirical_Df/',method[i],'.rds'))
    JoinedDfMicroglia <- JoinedDfMicroglia %>% dplyr::arrange(adjPvalue)
    ClusterAnnotResults <- GetClustersAnnotation(JoinedDfMicroglia%>% dplyr::filter(Biotype=='coding'),codingGenes,geneIdToName,librariesToRun,topPercentile = 0.1,save = T,outprefix = paste0(method[i],'_coding'))
    saveRDS(ClusterAnnotResults,file=paste0('../../Count_Data/PASCAL_Results/Annotations_New/Annotation_',method[i],'_coding.rds'))
    if(any(JoinedDfMicroglia$Biotype=='all')){
      ClusterAnnotResults <-  GetClustersAnnotation(JoinedDfMicroglia%>% dplyr::filter(Biotype=='all'),codingGenes,geneIdToName,librariesToRun,topPercentile = 0.1,save = T,outprefix = paste0(method[i],'_all'))
      saveRDS(ClusterAnnotResults,file=paste0('../../Count_Data/PASCAL_Results/Annotations_New/Annotation_',method[i],'_all.rds'))
    }
  }
}

#Main Function, For each method, obtain annotations of modules with top disease associations.
GetTopPercentileOfEachMethod <- function(topPercentile,methodPath=NULL,method=NULL){
  #Path to all dataframes created by each methods
  if(is.null(methodPath)){
    method <- c('JoinedDfMicrogliaWGCNAUnsigned',
                'JoinedDfMicrogliaWGCNATranscripts',
                'Microglia_Pearson_cor0p2_abs',
                'JoinedDfMicrogliaPearsonTranscripts')
    methodPath <- paste0('../../Count_Data/PASCAL_Results/Empirical_Df/',method,'.rds')
  }
  #Get annotation for each module
  for(i in 1:length(method)){
    print(method[i])
    JoinedDfMicroglia <- readRDS(methodPath[i])
    JoinedDfMicrogliaCoding <- JoinedDfMicroglia %>% dplyr::filter(Biotype=='coding') %>% dplyr::arrange(adjPvalue)
    ClusterAnnotResults <- GetAllClustersAnnotation(JoinedDfMicrogliaCoding %>% dplyr::arrange(adjPvalue) %>%
                                                    dplyr::filter(adjPvalue <= quantile(unique(adjPvalue),topPercentile)),
                                                 codingGenes,geneIdToName,librariesToRun,save = F)
    saveRDS(ClusterAnnotResults,file=paste0('../../Count_Data/PASCAL_Results/Annotations/Annotation_',method[i],'_Top50_coding.rds'))
    
    if(any(JoinedDfMicroglia$Biotype=='all')){
      JoinedDfMicrogliaAll <- JoinedDfMicroglia %>% dplyr::filter(Biotype=='all') %>% dplyr::arrange(adjPvalue)
      ClusterAnnotResults <- GetAllClustersAnnotation(JoinedDfMicrogliaAll %>% dplyr::arrange(adjPvalue) %>%
                                                      dplyr::filter(adjPvalue <= quantile(unique(adjPvalue),topPercentile)),
                                                      codingGenes,geneIdToName,librariesToRun,save = F)
      saveRDS(ClusterAnnotResults,file=paste0('../../Count_Data/PASCAL_Results/Annotations/Annotation_',method[i],'_Top50_all.rds'))
    }
  }
}

#Main Function, obtain annotations of modules with significant disease associations.
GetAnnotationForSignificantClusters <- function(JoinedDfMicroglia,isPath=T,pCutOff=0.05){
  if(isPath){
    load(JoinedDfMicroglia)
  }
  ClusterAnnotSigResults <- GetClustersAnnotation(JoinedDfMicroglia %>% filter(adjPvalue < pCutOff),codingGenes,
                                                  geneIdToName,librariesToRun,topPercentile = 1,save = F)
  library(gridExtra)
  dfFull <- do.call(rbind,lapply(ClusterAnnotSigResults,function(x) x$df))
  rownames(dfFull) <- as.character(seq(1,nrow(dfFull)))
  #Organize annotation result headers
  dfFull <- dfFull %>% dplyr::arrange(KEGG_2016.Term) %>% dplyr::mutate('KEGG Term'=sapply(KEGG_2016.Term,function(x) unlist(strsplit(x=x,split='_'))[1]))
  df <- dfFull %>% dplyr::select(adjP=adjPvalue,Study=StudyName,Size,Biotype,'KEGG Term',
                  'KEGG Overlap'=KEGG_2016.Overlap,'KEGG P'=KEGG_2016.Adjusted.P.value,
                  MicrogliaOverlap,Microglia_P,ADOverlap,AD_P)
  #Append Microglia gene enrichment information
  
  
  grid.arrange(tableGrob(df,theme = ttheme_default(base_size = 9)))
  
  return(dfFull)
}
#Main Function, obtain annotations of modules with significant disease associations.
GetAnnotationTranscriptClusters <- function(){
  load('../../Count_Data/PASCAL_Results/JoinedDfMicrogliaPathwayTranscript.rda')
  ClusterAnnotSigResults <- GetClustersAnnotation(JoinedDfMicrogliaTranscript %>% filter(adjPvalue < 0.2),codingGenes,
                                                  geneIdToName,librariesToRun,topPercentile = 1,save = F)
  return(ClusterAnnotSigResults)
}
