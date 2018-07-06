GetSingleClusterAnnotation <- function(genes,librariesToRun){
  invisible(capture.output(df <- enrichR::enrichr(genes = genes,databases = librariesToRun)))
  df <- do.call(rbind,df)
  df$Lib <- sapply(rownames(df),function(x){
    splitStr <- strsplit(x,'[.]')
    unlist(splitStr)[1]
  })
  return(df)
}
GetTopTermsForEachLib <- function(enrichrResult){
  libraries <- unique(enrichrResult$Lib)
  allDf <- vector(mode='list',length = length(libraries))
  for(i in 1:length(libraries)){
    lib <- libraries[i]
    curDf <- enrichrResult %>% dplyr::filter(Lib==lib) %>% 
      dplyr::arrange(Adjusted.P.value) %>% dplyr::select(Term,Adjusted.P.value,Overlap,Genes) %>% 
    {.[1,]}
    colnames(curDf) <- sapply(colnames(curDf),function(x)paste0(lib,'.',x))
    allDf[[i]] <- curDf
  }
  return(do.call(cbind,allDf))
}

GetClustersAnnotation <- function(JoinedDfMicroglia,geneIdToName,librariesToRun){
  uniqueStudies <- unique(JoinedDfMicroglia$StudyName)
  clusterResults <- vector(mode = 'list',length = length(uniqueStudies))
  names(clusterResults) <- uniqueStudies
  for(i in 1:length(uniqueStudies)){
    #Get all clusters, coding and all, as level
    currentStudyDf <- JoinedDfMicroglia %>% dplyr::filter(StudyName==uniqueStudies[i])
    #Remove duplicated clusters
    currentStudyDf <- currentStudyDf %>% distinct(chi2Pvalue,Size,Biotype,.keep_all = T)
    #Arrange by increasing p value
    currentStudyDf <- currentStudyDf %>% dplyr::arrange(adjPvalue)
    #Take top 20 ranked clusters. 
    topClusters <- currentStudyDf[1:floor(0.10*nrow(currentStudyDf)),]
    allGenes <- lapply(topClusters$Genes,function(x) geneIdToName[[x]])
    #Run enrichR for each cluster
    enrichrResults <- mclapply(allGenes,function(x) GetSingleClusterAnnotation(genes = x,
                                                                            librariesToRun = librariesToRun),mc.cores = 20)
    names(enrichrResults) <- topClusters$Name
    #Keep the top ranking term, p value, overlap and genes for the cluster.
    topForEachCluster <- do.call(rbind,lapply(enrichrResults,function(x) GetTopTermsForEachLib(x)))
    clusterResults[[i]] <- list(df=cbind(topClusters,topForEachCluster),enrichrResult=enrichrResults)
  }  
}
library(enrichR)
library(hashmap)
load('../../Count_Data/geneGtfTableFull.rda')
geneGtfTableFull <- geneGtfTableFull %>% dplyr::select(gene_name,gene_id)
geneIdToName <- hashmap::hashmap(keys = geneGtfTableFull$gene_id,values = geneGtfTableFull$gene_name)

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

librariesToRun <- c('KEGG_2016','TRANSFAC_and_JASPAR_PWMs','Jensen_DISEASES')
ClusterAnnotResults <- GetClustersAnnotation(JoinedDfMicroglia,geneIdToName,librariesToRun)