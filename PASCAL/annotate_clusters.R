library(gridExtra)
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

GetClustersAnnotation <- function(JoinedDfMicroglia,codingGenes,geneIdToName,librariesToRun,topPercentile=0.1,save = F){
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
      currentGenes <- allGenes[[j]]
      currentGenes <- currentGenes[which(codingGenes[[currentGenes]])]
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
# geneGtfTableFull <- LoadGTF(full=T) %>% dplyr::distinct(gene_id,.keep_all=T)
# geneGtfTableFull <- geneGtfTableFull %>% dplyr::select(gene_name,gene_id,gene_biotype)
load(file='../../Count_Data/geneGtfTableFull.rda')
geneIdToName <- hashmap::hashmap(keys = geneGtfTableFull$gene_id,values = geneGtfTableFull$gene_name)
codingGenes <- hashmap::hashmap(keys = geneGtfTableFull$gene_name,values = ifelse(geneGtfTableFull$gene_biotype=='protein_coding',T,F))
librariesToRun <- c('KEGG_2016','TRANSFAC_and_JASPAR_PWMs','Jensen_DISEASES')

GetTopPercentileOfEachStudy <- function(){
  # method <- c('CodingGenesMicroglia_Pearson_cor0p2_abs',
  #             'CodingGenesMicroglia_Jaccard_cor0p2_abs')
  method <- 'CodingWGCNAUnsigned_Soft4_Size3'
  for(i in 1:length(method)){
    load(paste0('../../Count_Data/PASCAL_Results/',method[i],'.rda'))
    JoinedDfMicroglia <- JoinedDfMicroglia %>% dplyr::arrange(adjPvalue)
    ClusterAnnotResults <- GetClustersAnnotation(JoinedDfMicroglia,codingGenes,geneIdToName,librariesToRun,topPercentile = 0.1,save = F)
    save(ClusterAnnotResults,file=paste0('../../GWAS/PASCAL_results/Annotation_',method[i],'.rda'))
  }
}
GetAnnotationForSignificantClusters <- function(JoinedDfMicroglia,isPath=T){
  if(isPath){
    load(JoinedDfMicroglia)
  }
  ClusterAnnotSigResults <- GetClustersAnnotation(JoinedDfMicroglia %>% filter(adjPvalue < 0.1),codingGenes,
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
GetAnnotationTranscriptClusters <- function(){
  load('../../Count_Data/PASCAL_Results/JoinedDfMicrogliaPathwayTranscript.rda')
  ClusterAnnotSigResults <- GetClustersAnnotation(JoinedDfMicrogliaTranscript %>% filter(adjPvalue < 0.2),codingGenes,
                                                  geneIdToName,librariesToRun,topPercentile = 1,save = F)
  return(ClusterAnnotSigResults)
}
# p0p05 <- '../../Count_Data/PASCAL_Results/JoinedDfMicrogliaPathwayGenes_0p05.rda'
# GetAnnotationForSignificantClusters(p0p05)
# 
# p0p01 <- '../../Count_Data/PASCAL_Results/JoinedDfMicrogliaPathwayGenes_0p01.rda'
# sigClusters0p01 <- do.call(rbind,lapply(GetAnnotationForSignificantClusters(p0p01),function(x) x$df))
# 
# signed0p05 <- '../../Count_Data/PASCAL_Results/JoinedDfMicrogliaPathwayGenesSigned_0p05.rda'
# sigClustersSigned0p05 <- do.call(rbind,lapply(GetAnnotationForSignificantClusters(signed0p05),function(x) x$df))



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

