GetSingleClusterAnnotation <- function(genes,librariesToRun){
  invisible(capture.output(df <- enrichR::enrichr(genes = genes,databases = librariesToRun)))
  df <- do.call(rbind,df)
  df$Lib <- sapply(rownames(df),function(x){
    splitStr <- strsplit(x,'[.]')
    unlist(splitStr)[1]
  })
  return(df)
}
GetTopTermsForEachLib <- function(enrichrResult,librariesToRun){
  allDf <- vector(mode='list',length = length(librariesToRun))
  if(nrow(enrichrResult) == 0){
    for(i in 1:length(librariesToRun)){
      lib <- librariesToRun[i]
      curDf <- data.frame(Term=NA,Adjusted.P.value=NA,Overlap=NA,Genes=NA) 
      colnames(curDf) <- sapply(colnames(curDf),function(x)paste0(lib,'.',x))
      allDf[[i]] <- curDf
    }
    return(do.call(cbind,allDf))
  }else{
    for(i in 1:length(librariesToRun)){
      lib <- librariesToRun[i]
      curDf <- enrichrResult %>% dplyr::filter(Lib==lib) %>% 
        dplyr::arrange(Adjusted.P.value) %>% dplyr::select(Term,Adjusted.P.value,Overlap,Genes) %>% 
        {.[1,]}
      colnames(curDf) <- sapply(colnames(curDf),function(x)paste0(lib,'.',x))
      allDf[[i]] <- curDf
    }
    return(do.call(cbind,allDf))
  }
}

GetClustersAnnotation <- function(JoinedDfMicroglia,geneIdToName,librariesToRun,topPercentile=0.1,save = F){
  uniqueStudies <- unique(JoinedDfMicroglia$StudyName)
  clusterResults <- vector(mode = 'list',length = length(uniqueStudies))
  names(clusterResults) <- uniqueStudies
  for(i in 1:length(uniqueStudies)){
    print(uniqueStudies[i])
    #Get all clusters, coding and all, as level
    currentStudyDf <- JoinedDfMicroglia %>% dplyr::filter(StudyName==uniqueStudies[i])
    #Arrange by increasing p value
    currentStudyDf <- currentStudyDf %>% dplyr::arrange(adjPvalue)
    #Take top 10 ranked clusters. 
    topClusters <- currentStudyDf[1:floor(topPercentile*nrow(currentStudyDf)),]
    allGenes <- lapply(topClusters$Genes,function(x) geneIdToName[[x]])
    topClusters$GeneNames <- sapply(allGenes,function(x) paste(x,collapse = ' '))
    #Run enrichR for each cluster
    enrichrResults <- vector(mode='list',length = length(allGenes))
    for(j in 1:length(enrichrResults)){
      print(j/length(enrichrResults))
      enrichrResults[[j]] <- GetSingleClusterAnnotation(genes = allGenes[[j]],
                                                        librariesToRun = librariesToRun)
    }
    names(enrichrResults) <- topClusters$Name
    #Keep the top ranking term, p value, overlap and genes for the cluster.
    topForEachCluster <- do.call(rbind.fill,lapply(enrichrResults,function(x) GetTopTermsForEachLib(x,librariesToRun)))
    currentStudyResult <- list(df=cbind(topClusters,topForEachCluster),enrichrResult=enrichrResults)
    if(save){
      save(currentStudyResult,file=paste0('../../GWAS/PASCAL_results/microglia_gene/ClusterAnnot_',uniqueStudies[i],'.rda'))
    }
    clusterResults[[i]] <- currentStudyResult
  }
  return(clusterResults)
}
library(enrichR)
library(hashmap)
library(plyr)
library(dplyr)
source('../BatchCorrection_and_QC/load_GTF.R')
geneGtfTableFull <- LoadGTF(full=T) %>% dplyr::distinct(gene_id,.keep_all=T)
geneGtfTableFull <- geneGtfTableFull %>% dplyr::select(gene_name,gene_id)
geneIdToName <- hashmap::hashmap(keys = geneGtfTableFull$gene_id,values = geneGtfTableFull$gene_name)
librariesToRun <- c('KEGG_2016','TRANSFAC_and_JASPAR_PWMs','Jensen_DISEASES')

GetTopPercentileOfEachStudy <- function(){
  load('../../Count_Data/JoinedDfMicroglia.rda')
  JoinedDfMicroglia <- JoinedDfMicroglia %>% dplyr::filter(Biotype=='coding' & Level < 10) %>%
    dplyr::arrange(adjPvalue)
  ClusterAnnotResults <- GetClustersAnnotation(JoinedDfMicroglia,geneIdToName,librariesToRun,topPercentile = 0.1,save = T)
  save(ClusterAnnotResults,file='../../GWAS/PASCAL_results/microglia_gene/ClusterAnnot.rda')
}
GetAnnotationForSignificantClusters <- function(){
  load('../../Count_Data/PASCAL_Results/JoinedDfMicrogliaPathwayGenes.rda')
  #Keep only clusters greater than 10 and less than 200
  JoinDfMicroglia <- JoinDfMicroglia %>% dplyr::filter(Size > 10 & Size < 200 & Level > 3)
  ClusterAnnotSigResults <- GetClustersAnnotation(JoinedDfMicroglia %>% filter(adjPvalue < 0.1),
                                                  geneIdToName,librariesToRun,topPercentile = 1,save = F)
  return(ClusterAnnotSigResults)
}
GetAnnotationTranscriptClusters <- function(){
  load('../../Count_Data/PASCAL_Results/JoinedDfMicrogliaPathwayTranscript.rda')
  ClusterAnnotSigResults <- GetClustersAnnotation(JoinedDfMicrogliaTranscript %>% filter(adjPvalue < 0.2),
                                                  geneIdToName,librariesToRun,topPercentile = 1,save = F)
  return(ClusterAnnotSigResults)
}


# allEnrichRLib <- enrichR::listEnrichrDbs()$libraryName 
# #Remove duplicate libraries of different versions.
# allEnrichRLib <- allEnrichRLib[allEnrichRLib!='GeneSigDB' & allEnrichRLib!='Genes_Associated_with_NIH_Grants']
# allVersions <- sapply(allEnrichRLib,function(x){
#  splitString <- unlist(strsplit(x,'_'))
#  return(substr(splitString[length(splitString)],1,2) == '20')
# })
# allEnrichRLibNoVersion <- allEnrichRLib[!allVersions]
# allEnrichRLibWithVersion <- allEnrichRLib[allVersions]
# 
# suffixToKeep <- c()
# for(str in allEnrichRLibWithVersion){
#   splitString <- unlist(strsplit(str,'_'))
#   preFix <- paste0(splitString[1:(length(splitString)-1)],collapse = '')
#   allMatchingPrefixes <- grepl(preFix,allEnrichRLibWithVersion)
#   allMatchingSuffix <- sapply(allEnrichRLibWithVersion[allMatchingPrefixes],function(x){
#     splitString <- unlist(strsplit(x,'_'))
#     return(splitString[length(splitString)])
#   })
#   suffixToKeep <- c(suffixToKeep,
#                     which(names(allMatchingSuffix)[which.max(allMatchingSuffix)] == allEnrichRLibWithVersion))
# }
# suffixToKeep <- unique(suffixToKeep)
# allEnrichRLib <- c(allEnrichRLibNoVersion,allEnrichRLibWithVersion[suffixToKeep],'GO_Biological_Process_2017',
#                    'GO_Biological_Process_2017b','GO_Cellular_Component_2017','GO_Cellular_Component_2017b',
#                    'GO_Molecular_Function_2017','GO_Molecular_Function_2017b')

