source('pval_correction.R')

ParsePASCALFilePermuted <- function(resultPaths,clusterPaths,appendInfo=T,empirical){
  #Load all PASCAL results
  studyNames <- sapply(resultPaths,function(x){
    splitStr <- unlist(strsplit(x,'/'));
    fileName <- splitStr[length(splitStr)];
    fileNameSplit <- unlist(strsplit(fileName,'[.]'));
    endIndex <- grep('PathwaySet',fileNameSplit);
    studyName <- paste(fileNameSplit[1:(endIndex-1)],collapse = '_');
    return(studyName)
  })
  resultDf <- do.call(rbind,lapply(1:length(resultPaths),function(i) {
    df <- data.table::fread(resultPaths[[i]],header = T);
    df$StudyName <- rep(studyNames[i],nrow(df));
    return(df)}))
  resultDf <- resultDf %>%
    dplyr::mutate(Perm = unlist(lapply(Name,function(x) as.numeric(gsub(x=x,pattern = '.*_perm_',replacement = '')))),
                  RealName = unlist(lapply(Name,function(x) gsub(x=x,pattern = '_perm_.*',replacement = ''))))
  #For each cluster, load it's respective gene ID
  clustersInfo <-  unlist(lapply(clusterPaths,function(x) qusage::read.gmt(file=x)),recursive = F)
  clustersDf <- data_frame(ClusterName = names(clustersInfo),Genes = clustersInfo)
  
  #Join result and gene info for each cluster. Filter for those with P-values. Then, extract metainfo from cluster Name (level,type)
  joinedDf <- dplyr::left_join(resultDf,clustersDf,by=c('Name'='ClusterName')) %>% filter(!is.na(chi2Pvalue)) %>%
    dplyr::mutate(Size=sapply(Genes,length)) %>%  #Size info
    dplyr::mutate(Biotype=factor(sapply(Name,function(x) unlist(strsplit(x,'_'))[1])))  %>% #Biotype
    dplyr::mutate(CellType=factor(sapply(Name,function(x) unlist(strsplit(x,'_'))[2])))  %>% #Cell Type
    dplyr::mutate(AggregationType=factor(sapply(Name,function(x) unlist(strsplit(x,'_'))[3]))) %>% #Gene or transcript
    dplyr::mutate(Level=sapply(Name,function(x){
      splitString <- unlist(strsplit(x,'_'));
      return(as.numeric(splitString[which(splitString=='level')+1]))})) #Gene or transcript
  #Have a gene Name column for convenience
  load(file='../../Count_Data/geneGtfTableFull.rda')
  geneIdToName <- hashmap::hashmap(keys = geneGtfTableFull$gene_id,values = geneGtfTableFull$gene_name)
  library(parallel)
  joinedDf$GeneNames <- mclapply(joinedDf$Genes,function(x) geneIdToName[[x]],mc.cores = 1)
  codingGenesInNetwork <- unique(unlist(joinedDf %>% dplyr::filter(Biotype=='coding') %>% {.$GeneNames}))
  allGenesInNetwork <- unique(unlist(joinedDf %>% dplyr::filter(Biotype=='all') %>% {.$GeneNames}))
  
  if(appendInfo){
    source('append_gene_enrichment.R')
    joinedDf <- AppendGeneEnrichment(joinedDf,codingGenesInNetwork=codingGenesInNetwork,allGenesInNetwork=allGenesInNetwork)
    if(empirical){
      source('append_open_target_empirical.R')
      joinedDf <- AppendOpenTargetEmpirical(joinedDf,
                                            permPath = '../../Count_Data/OpenTarget/',
                                            csvPath = '/local/data/public/zmx21/zmx21_private/GSK/OpenTargets_scores/',
                                            codingGenesInNetwork=codingGenesInNetwork,
                                            allGenesInNetwork=allGenesInNetwork)
      
    }else{
      source('append_open_target.R')
      joinedDf <- AppendOpenTarget(joinedDf,
                                   csvPath = '/local/data/public/zmx21/zmx21_private/GSK/OpenTargets_scores/',
                                   codingGenesInNetwork=codingGenesInNetwork,
                                   allGenesInNetwork=allGenesInNetwork)
      
    }
  }
  return(joinedDf)
}
#Get path to all studies in a directory
GetPathToPASCALResults <- function(PASCALResultPath){
  allStudies <- dir(PASCALResultPath)
  allStudies <- allStudies[grep('PathwaySet',allStudies)]
  allStudiesPath <- sapply(allStudies,function(x) paste0(PASCALResultPath,x))
  return(allStudiesPath)
}

LoadPASCALResultsPermuted <- function(filter=T,resultPaths,clusterPaths,sizeLower=10,sizeUpper=300,appendInfo=T,empirical=F){
  source('pval_correction.R')
  microgliaAllGenesPath <- GetPathToPASCALResults(resultPaths[1])
  JoinedDfMicroglia <- ParsePASCALFilePermuted(microgliaAllGenesPath,clusterPaths[1],appendInfo,empirical)
  if(length(resultPaths)== 2 & length(clusterPaths) == 2){
    microgliaCodingGenesPath <- GetPathToPASCALResults(resultPaths[2])
    JoinedDfMicroglia <- rbind(JoinedDfMicroglia,ParsePASCALFilePermuted(microgliaCodingGenesPath,clusterPaths[2],appendInfo,empirical))
  }
  JoinedDfMicroglia <- JoinedDfMicroglia %>% dplyr::filter(Size > sizeLower & Size < sizeUpper)
  JoinedDfMicroglia <- AppendCorrectedPVal(JoinedDfMicroglia)
  if(filter){
    #Remove repeating rows
    JoinedDfMicroglia <- JoinedDfMicroglia %>% dplyr::distinct(StudyName,chi2Pvalue,Size,Biotype,.keep_all=T)
  }
  return(JoinedDfMicroglia)
}
JoinedDfMicrogliaPearsonPermuted<- LoadPASCALResultsPermuted(filter=T,c('../../GWAS/PASCAL_results/Permuted_Pearson_100/permuted_Coding_Pearson/'),
                                                         c('../../genesets/Permuted_Pearson_100/permuted_Coding_Pearson.gmt'),
                                                     sizeLower = 3,sizeUpper = 300,appendInfo = F,empirical = F)
