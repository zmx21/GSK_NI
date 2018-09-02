library(qusage)
library(dplyr)
library(ggplot2)
library(parallel)
source('annotate_clusters.R')
source('pval_correction.R')
source('write_GWAS_metadata.R')
#Merged result from PASCAL with cluster info
ParsePASCALFile <- function(resultPaths,clusterPaths,appendInfo=T,empirical){
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
  joinedDf$GeneNames <- lapply(joinedDf$Genes,function(x) geneIdToName[[x]])
  codingGenesInNetwork <- unique(unlist(joinedDf %>% dplyr::filter(Biotype=='coding') %>% {.$GeneNames}))
  allGenesInNetwork <- unique(unlist(joinedDf %>% dplyr::filter(Biotype=='all') %>% {.$GeneNames}))
  
  #Append additional information
  if(appendInfo){
    #Append microglia gene enrichment
    source('append_gene_enrichment.R')
    joinedDf <- AppendGeneEnrichment(joinedDf,codingGenesInNetwork=codingGenesInNetwork,allGenesInNetwork=allGenesInNetwork)
    #Appen Open Target enrichment, based on empirical or hypergeometric.
    if(empirical){
      source('append_open_target_empirical.R')
      joinedDf <- AppendOpenTargetEmpirical(joinedDf,
                                            permPath = '../../Count_Data/OpenTarget/',
                                            csvPath = '../../OpenTargets_scores/',
                                            codingGenesInNetwork=codingGenesInNetwork,
                                            allGenesInNetwork=allGenesInNetwork)
      
    }else{
      source('append_open_target.R')
      joinedDf <- AppendOpenTarget(joinedDf,
                                   csvPath = '../../OpenTargets_scores/',
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

#Plots of p-value distribution at different levels
ConstructPlots <- function(JoinedDfMicroglia){
  library(egg)
  
  #plot number of clusters at different levels.
  levelVSNumber <- JoinedDfMicroglia %>% dplyr::filter(StudyName=='AB42_Deming')%>% dplyr::select(Level,Biotype)
  tiff('../../Figures/PASCAL/NumOfClustersVSLevel.tiff',width = 1200,height = 500)
  p1 <- ggplot(levelVSNumber) + aes(fill=factor(Biotype),x=factor(Level)) + geom_bar(position = 'dodge') + 
    scale_y_continuous(breaks = round(seq(0, 1700, by = 100),1),limits = c(0,1650)) + labs(y='Number of Clusters',x='Louvain Recursion Level',fill='Biotype')
  ggarrange(p1)
  dev.off()
  
  #plot level of clustering vs size of clusters.
  levelvsSize <- JoinedDfMicroglia %>% dplyr::select(Level,Size,Biotype)
  levelvsSize <- levelvsSize[!duplicated(levelvsSize),] 
  p1 <- ggplot(levelvsSize) + aes(x=factor(Level),y=Size,fill=factor(Biotype)) + geom_boxplot(width=0.3,position='dodge') +
    scale_y_continuous(breaks = round(seq(0, 600, by = 20),1),limits = c(0,600)) + ggtitle('Distribution of Cluster Sizes at Each Level') + 
    xlab('Louvain Recursion Level') + ylab('Cluster Size') + labs(fill='Biotype')
  tiff('../../Figures/PASCAL/SizeVSLevel.tiff',width = 1200,height = 500)
  ggarrange(p1)
  dev.off()
  
  #Plot density of P-Value for different Levels
  PlotDensityByLevel <- function(joinedDf,title,pValType='emp'){
    if(pValType=='emp'){
      return(ggplot(joinedDf,aes(x=-1*log10(empPvalue))) + geom_density(trim=T,fill="#4271AE") + 
               facet_grid(Level ~ Biotype,labeller = labeller(Biotype = as_labeller(c(coding='Coding Genes',all='All Genes')),
                                                              Level = as_labeller(c(`4`='Level 4',`5`='Level 5',`6`='Level 6',`7`='Level 7',`8`='Level 8',`9`='Level 9',`10`='Level 10')))) + 
               theme(strip.text = element_text(size=15))+ xlim(1,5) + labs(x='-log10(Empirical P-Value)',y='Density',fill='Louvain\nRecursion\nLevel') + ggtitle(title))
    }else if (pValType=='fdr'){
      return(ggplot(joinedDf,aes(x=adjPvalue)) + geom_density(trim=T,fill="#4271AE") + geom_vline(xintercept=0.05,colour='red') +
               facet_grid(Level ~ Biotype,labeller = labeller(Biotype = as_labeller(c(coding='Coding Genes',all='All Genes')),
                                                              Level = as_labeller(c(`4`='Level 4',`5`='Level 5',`6`='Level 6',`7`='Level 7',`8`='Level 8',`9`='Level 9',`10`='Level 10')))) + 
               theme(strip.text = element_text(size=15))+ xlim(0,1) + labs(x='BH Adjusted P-Value',y='Density',fill='Louvain\nRecursion\nLevel') + ggtitle(title))
      
    }
  } 
  
  #Look at Microglia distribution across levels, for all studies
  p1 <- PlotDensityByLevel(JoinedDfMicroglia %>% dplyr::filter(Level >= 4 & Level < 11),
                           title = c('Microglia Coding Genes-\nP-val Density Across Levels for all Studies'))
  tiff(filename = '../../Figures/PASCAL/LevelVSPVal.tiff',width = 1200,height = 500)
  ggarrange(p1,nrow=1)
  dev.off()
  
  
  #Plot density for inidividual studies
  library(scales)
  allStudies <- lapply(1:length(unique(JoinedDfMicroglia$StudyName)),function(i) PlotDensityByLevel(JoinedDfMicroglia %>% 
                                                                                                            dplyr::filter((Level %in% c(4,5,6,7)) & StudyName == unique(JoinedDfMicroglia$StudyName)[i]),
                                                                                                          title = paste0('Microglia Genes -\n',unique(JoinedDfMicroglia$StudyName)[i])) + scale_y_continuous(breaks = pretty_breaks(n = 3)))
  tiff(filename = '../../Figures/PASCAL/PValIndivStudies.tiff',width = 1600,height = 600)
  ggarrange(plots = allStudies,ncol = 4,nrow = 2)
  dev.off()
  
  tiff(filename = '../../Figures/PASCAL/PValAdjIndivStudies.tiff',width = 1600,height = 600)
  allStudiesFDR <- lapply(1:length(unique(JoinedDfMicroglia$StudyName)),function(i) PlotDensityByLevel(JoinedDfMicroglia %>% 
                                                                                                      dplyr::filter((Level %in% c(4,5,6,7)) & StudyName == unique(JoinedDfMicroglia$StudyName)[i]),
                                                                                                    title = paste0('Microglia Genes AdjPVal-\n',unique(JoinedDfMicroglia$StudyName)[i]),pValType = 'fdr') + scale_y_continuous(breaks = pretty_breaks(n = 3)))
  ggarrange(plots=allStudiesFDR,nrow=2,ncol = 4)
  dev.off()
  
  #Plot density vs random cluster
  allStudiesRandomMicrogliaPath <- GetPathToPASCALResults('../../GWAS/PASCAL_results/randommicroglia_gene/')
  JoinedDfMicrogliaRandom <- ParsePASCALFile(allStudiesRandomMicrogliaPath,'../../Louvain_results/randommicroglia_clusters.gmt')
  JoinedDfMicrogliaRandom <- AppendCorrectedPVal(JoinedDfMicrogliaRandom)
  source('geom_split_violin.R')
  tiff(filename = '../../Figures/PASCAL/RandomVSTrueNetwork.tiff',width = 1200,height = 600)
  ggplot(rbind(JoinedDfMicrogliaRandom,JoinedDfMicroglia) %>% 
           dplyr::filter(Biotype!='randomall' & Biotype!='all' & Level > 3 & Level <= 7) %>% 
           dplyr::mutate(Type=ifelse(Biotype=='randomcoding','Random','True')), 
         aes(x=adjPvalue)) + geom_density(trim=T,fill="#4271AE") + facet_grid(Level ~ Type,labeller = labeller(Type = as_labeller(c(Random='Random Network',True='True Network')),
                                                                                             Level = as_labeller(c(`4`='Level 4',`5`='Level 5',`6`='Level 6',`7`='Level 7',`8`='Level 8',`9`='Level 9',`10`='Level 10')))) + 
     xlab('BH Adjusted P-Value') + theme(strip.text = element_text(size=15))+
    labs(fill='Random \n or \nTrue') + ggtitle('Randomly permuted network VS True Network')
  dev.off()
  
  #Return significant clusters, remove those with exact same genes (chi2PValue) and in the same study with the same biotype
  SignificantClustersFDR5 <- JoinedDfMicroglia %>% dplyr::filter(adjPvalue<=0.05 & Level > 3) %>% 
    dplyr::distinct(chi2Pvalue,Size,StudyName,.keep_all=T)
  SignificantClustersFDR10 <- JoinedDfMicroglia %>% dplyr::filter(adjPvalue<=0.1 & Level > 3) %>% 
    dplyr::distinct(chi2Pvalue,Size,StudyName,.keep_all=T)
  
  #P values against study name
  
}
#compare distribution of true p-value and permuted network p-value of GWAS 
CompareTrueWithRandom <- function(global=F){
  source('pval_correction.R')
  allStudiesMicrogliaPath <- GetPathToPASCALResults('../../GWAS/PASCAL_results/microglia_gene/')
  JoinedDfMicroglia <- ParsePASCALFile(allStudiesMicrogliaPath,'../../Louvain_results/microglia_gene_clusters.gmt')
  JoinedDfMicroglia <- AppendCorrectedPVal(JoinedDfMicroglia)
  allStudiesRandomMicrogliaPath <- GetPathToPASCALResults('../../GWAS/PASCAL_results/randommicroglia_gene/')
  JoinedDfMicrogliaRandom <- ParsePASCALFile(allStudiesRandomMicrogliaPath,'../../Louvain_results/randommicroglia_clusters.gmt')
  JoinedDfMicrogliaRandom <- AppendCorrectedPVal(JoinedDfMicrogliaRandom)
  
  #Look at different in quantile between random and true network Pval Distribution
  library(WRS)
  if(!global){
    uniqueLevels <- unique(JoinedDfMicrogliaRandom$Level)
    print('Testing Random VS True')
    quantileTests <- mclapply(1:length(uniqueLevels),function(i) {
      curLevel <- JoinedDfMicroglia %>% filter(Level== uniqueLevels[i] & Biotype=='coding');
      curLevelRandom <- JoinedDfMicrogliaRandom %>% filter(Level==uniqueLevels[i] & Biotype=='randomcoding');
      qcomhd(curLevel$adjPvalue,curLevelRandom$adjPvalue, q=seq(0.001,0.005,0.001),plotit = F);
    },mc.cores=length(uniqueLevels))
    save(quantileTests,file='../../Louvain_results/RandomVsTrue.rda')
  }else{
    JoinedDfMicroglia <- JoinedDfMicroglia %>% filter(Biotype=='coding')
    JoinedDfMicrogliaRandom <- JoinedDfMicrogliaRandom %>% filter(Biotype=='randomcoding')
    quantileTest <- qcomhd(JoinedDfMicroglia$adjPvalue,JoinedDfMicrogliaRandom$adjPvalue,q=seq(0.001,0.005,0.001))
    save(quantileTest,file='../../Louvain_results/RandomVsTrueGlobal.rda')
  }
}

#Main function, which loads file output of PASCAL, and construct a dataframe. 
#option allows append other info (KEGG, Open Target), and Open Target can be generated empirically.
LoadPASCALResults <- function(filter=T,resultPaths,clusterPaths,sizeLower=10,sizeUpper=200,appendInfo=T,empirical=F){
  source('pval_correction.R')
  microgliaAllGenesPath <- GetPathToPASCALResults(resultPaths[1])
  JoinedDfMicroglia <- ParsePASCALFile(microgliaAllGenesPath,clusterPaths[1],appendInfo,empirical)
  if(length(resultPaths)== 2 & length(clusterPaths) == 2){
    microgliaCodingGenesPath <- GetPathToPASCALResults(resultPaths[2])
    JoinedDfMicroglia <- rbind(JoinedDfMicroglia,ParsePASCALFile(microgliaCodingGenesPath,clusterPaths[2],appendInfo,empirical))
  }
  # JoinedDfMicroglia <- ParsePASCALFile(GetPathToPASCALResults('../../GWAS/PASCAL_results/Old/microglia_gene/'),'../../Louvain_results/Old/microglia_gene_clusters.gmt')
  JoinedDfMicroglia <- JoinedDfMicroglia %>% dplyr::filter(Size > sizeLower & Size < sizeUpper)
  JoinedDfMicroglia <- AppendCorrectedPVal(JoinedDfMicroglia)
  if(filter){
    #Remove repeating rows
    JoinedDfMicroglia <- JoinedDfMicroglia %>% dplyr::distinct(StudyName,chi2Pvalue,Size,Biotype,.keep_all=T)
  }
  return(JoinedDfMicroglia)
}
#Same as LoadPASCALResults, but for transcripts
LoadPASCALTranscriptResults <- function(filter=T){
  source('pval_correction.R')
  microgliaAllGenesPath <- GetPathToPASCALResults('../../GWAS/PASCAL_results/microglia_coding_transcripts/')
  JoinedDfMicroglia <- ParsePASCALFile(microgliaAllGenesPath,'../../Louvain_results/CodingMicrogliaTranscripts/CodingMicrogliaTranscriptsAsGenes.gmt')

  JoinedDfMicroglia <- AppendCorrectedPVal(JoinedDfMicroglia)
  if(filter){
    #Remove repeating rows
    JoinedDfMicroglia <- JoinedDfMicroglia %>% dplyr::distinct(StudyName,chi2Pvalue,Size,Biotype,.keep_all=T)
  }
  return(JoinedDfMicroglia)
}
#Append OpenTarget results to an existing dataframe
ViewOpenTargetResults <- function(JoinedDfMicroglia,type,empirical=F){
  round_df <- function(x, digits) {
    # round all numeric variables
    # x: data frame 
    # digits: number of digits to round
    numeric_columns <- sapply(x, mode) == 'numeric'
    x[numeric_columns] <-  round(x[numeric_columns], digits)
    return(x)
  }
  #Append empricialy results if specified.
  if(empirical){
    source('append_open_target_empirical.R')
    JoinedDfMicroglia <- JoinedDfMicroglia %>% dplyr::select(-contains('_drug_P'),
                                                             -contains('_drug_overlap'),
                                                             -contains('_target_P'),
                                                             -contains('_target_overlap'),
                                                             -contains('_target_Genes'),
                                                             -contains('_drug_Genes'))
    load('../../Count_Data/PASCAL_Results/Microglia_Pearson_cor0p2_abs.rda')
    codingGenesInNetwork <- unique(unlist(JoinedDfMicrogliaPearson %>% dplyr::filter(Biotype=='coding') %>% {.$GeneNames}))
    allGenesInNetwork <- unique(unlist(JoinedDfMicrogliaPearson %>% dplyr::filter(Biotype=='all') %>% {.$GeneNames}))
    
    JoinedDfMicroglia <- AppendOpenTargetEmpirical(JoinedDfMicroglia,
                                                   permPath = '../../Count_Data/OpenTarget/',
                                                   csvPath = '../../OpenTargets_scores/',
                                                   codingGenesInNetwork=codingGenesInNetwork,
                                                   allGenesInNetwork=allGenesInNetwork)
  }
  #Construct table to show results, along with KEGG Term.
  if(type=='target'){
    df <- JoinedDfMicroglia %>% dplyr::select(StudyName,Size,adjPvalue,'KEGG_Term'=KEGG_2016.Term,'KEGG_Overlap'=KEGG_2016.Overlap,'KEGG_P'=KEGG_2016.Adjusted.P.value,
                                              Microglia_Overlap=MicrogliaOverlap,Microglia_P,
                                              contains('target'),-contains('target_Genes')) %>% dplyr::mutate('KEGG_Term'=sapply(KEGG_Term,function(x) unlist(strsplit(x=x,split='_'))[1]))
    df <- dplyr::arrange(df,KEGG_Term)
    colnames(df) <- sapply(colnames(df),function(x) gsub(pattern = '_',replacement = '\n',x = x))
    df <- round_df(df,3)
    rownames(df) <- c()
    #set all cells as black
    cols <- matrix("black", nrow(df), ncol(df))
    
    #Set significant cells as red.
    pColumns <- sapply(colnames(df),function(x) grepl(pattern = '\nP',x = x))
    pColumns[3] <- F
    sigCells <- which(df<0.05,arr.ind = T)
    sigCells <- sigCells[sigCells[,2]%in%which(pColumns),]   #Only keep interested columns
    if(is.null(nrow(sigCells))){
      sigCells <- t(as.matrix(sigCells))
    }
    for(i in 1:nrow(sigCells)){
      cols[sigCells[i,1],sigCells[i,2]] <- c("red")
    }

  }else if(type=='drug'){
    df <- JoinedDfMicroglia %>% dplyr::select(StudyName,Size,adjPvalue,'KEGG_Term'=KEGG_2016.Term,'KEGG_Overlap'=KEGG_2016.Overlap,'KEGG_P'=KEGG_2016.Adjusted.P.value,
                                              Microglia_Overlap=MicrogliaOverlap,Microglia_P,
                                              contains('drug'),-contains('drug_Genes'))%>% dplyr::mutate('KEGG_Term'=sapply(KEGG_Term,function(x) unlist(strsplit(x=x,split='_'))[1]))
    df <- dplyr::arrange(df,KEGG_Term)
    colnames(df) <- sapply(colnames(df),function(x) gsub(pattern = '_',replacement = '\n',x = x))
    df <- round_df(df,3)
    rownames(df) <- c()
    #set all cells as black
    cols <- matrix("black", nrow(df), ncol(df))
    #Set significant cells as red.
    pColumns <- sapply(colnames(df),function(x) grepl(pattern = '\nP',x = x))
    pColumns[3] <- F
    sigCells <- which(df<0.05,arr.ind = T)
    sigCells <- sigCells[sigCells[,2]%in%which(pColumns),]   #Only keep interested columns
    if(is.null(nrow(sigCells))){
      sigCells <- t(as.matrix(sigCells))
    }
    for(i in 1:nrow(sigCells)){
      cols[sigCells[i,1],sigCells[i,2]] <- c("red")
    }
    
  }
  grid.arrange(tableGrob(df,theme = ttheme_default(base_size = 9,core=list(fg_params = list(col = cols),
                                                                           bg_params = list(col=NA)))))
}
###################PEARSON Jaccard################################

# JoinedDfMicrogliaJaccard <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results/Jaccard_Cor0p2/CodingGenesMicroglia_Jaccard_cor0p2_abs/',
#                                                   '../../GWAS/PASCAL_results/Jaccard_Cor0p2/AllGenesMicroglia_Jaccard_cor0p2_abs/'),
#                                        c('../../Louvain_results/Jaccard_Cor0p2/CodingGenesMicroglia_Jaccard_cor0p2_abs.gmt',
#                                          '../../Louvain_results/Jaccard_Cor0p2/AllGenesMicroglia_Jaccard_cor0p2_abs.gmt'),sizeLower = 3,sizeUpper = 300)
# JoinedDfMicroglia <- JoinedDfMicrogliaJaccard %>% dplyr::filter(Biotype=='coding')
# save(JoinedDfMicroglia,file='../../Count_Data/PASCAL_Results/CodingGenesMicroglia_Jaccard_cor0p2_abs.rda')
# sigClustersCodingJaccard <- do.call(rbind,lapply(GetAnnotationForSignificantClusters(JoinedDfMicrogliaJaccard %>% dplyr::filter(Biotype=='coding')
#                                                                               ,isPath = F),function(x) x$df))
# sigClustersAllJaccard <- do.call(rbind,lapply(GetAnnotationForSignificantClusters(JoinedDfMicrogliaJaccard %>% dplyr::filter(Biotype=='all')
#                                                                            ,isPath = F),function(x) x$df))
# GetGWASMetdata(rbind(sigClustersCodingJaccard,sigClustersAllJaccard))

###################PEARSON GENES################################

# JoinedDfMicrogliaPearson <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results/Pearson_Cor0p2/CodingGenesMicroglia_Pearson_cor0p2_abs/',
#                                                   '../../GWAS/PASCAL_results/Pearson_Cor0p2/AllGenesMicroglia_Pearson_cor0p2_abs/'),
#                                        c('../../Louvain_results/Pearson_Cor0p2/CodingGenesMicroglia_Pearson_cor0p2_abs.gmt',
#                                          '../../Louvain_results/Pearson_Cor0p2/AllGenesMicroglia_Pearson_cor0p2_abs.gmt'),sizeLower = 3,sizeUpper = 300)

# save(JoinedDfMicrogliaPearson,file='../../Count_Data/PASCAL_Results/Microglia_Pearson_cor0p2_abs.rda')
# load('../../Count_Data/PASCAL_Results/Microglia_Pearson_cor0p2_abs.rda')
# sigClustersCodingPearson <- GetAnnotationForSignificantClusters(JoinedDfMicrogliaPearson %>% dplyr::filter(Biotype=='coding'),isPath = F)
# save(sigClustersCodingPearson,file='../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingPearson.rda')
# ViewOpenTargetResults(sigClustersCodingPearson,type='target',empirical = T)
# ViewOpenTargetResults(sigClustersCodingPearson,type='drug')
# sigClustersAllPearson <- GetAnnotationForSignificantClusters(JoinedDfMicrogliaPearson %>% dplyr::filter(Biotype=='all'),isPath = F)
# save(sigClustersAllPearson,file='../../Count_Data/PASCAL_Results/sigClusters/sigClustersAllPearson.rda')
# ViewOpenTargetResults(sigClustersAllPearson,type='target',empirical = T)
# ViewOpenTargetResults(sigClustersAllPearson,type='drug')
# GetGWASMetdata(rbind(sigClustersCodingPearson,sigClustersAllPearson))

###################WGCNA GENES################################

# JoinedDfMicrogliaWGCNAUnsigned <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results/WGCNA_size3/CodingWGCNAUnsigned_Soft4_Size3//',
#                                                              '../../GWAS/PASCAL_results/WGCNA_size3/AllWGCNAUnsigned_Soft4_Size3//'),
#                                                   c('../../GWAS/PASCAL_New/resources/genesets/WGCNA_size3/CodingWGCNAUnsigned_Soft4_Size3.gmt',
#                                                     '../../GWAS/PASCAL_New/resources/genesets/WGCNA_size3/AllWGCNAUnsigned_Soft4_Size3.gmt'),sizeLower = 3,sizeUpper = 100)
# save(JoinedDfMicrogliaWGCNAUnsigned,file='../../Count_Data/PASCAL_Results/JoinedDfMicrogliaWGCNAUnsigned.rda')
# 
# sigClustersCodingWGCNAUnsigned <- GetAnnotationForSignificantClusters(JoinedDfMicrogliaWGCNAUnsigned %>% dplyr::filter(Biotype=='coding'),isPath = F)
# save(sigClustersCodingWGCNAUnsigned,file='../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingWGCNAUnsigned.rda')
# 
# ViewOpenTargetResults(sigClustersCodingWGCNAUnsigned,type='target',empirical = T)
# ViewOpenTargetResults(sigClustersCodingWGCNAUnsigned,type='drug')
# sigClustersAllWGCNAUnsigned <- GetAnnotationForSignificantClusters(JoinedDfMicrogliaWGCNAUnsigned %>% dplyr::filter(Biotype=='all'),isPath = F)
# save(sigClustersAllWGCNAUnsigned,file='../../Count_Data/PASCAL_Results/sigClusters/sigClustersAllWGCNAUnsigned.rda')
# 
# ViewOpenTargetResults(sigClustersAllWGCNAUnsigned,type='target',empirical = T)
# ViewOpenTargetResults(sigClustersAllWGCNAUnsigned,type='drug')
# GetGWASMetdata(rbind(sigClustersCodingWGCNAUnsigned,sigClustersCodingWGCNAUnsigned))

###################Non Neurological Genes################################

# JoinedDfNonNeurologicalPearson <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results_non_neurological/Pearson_Cor0p2_Pathway/CodingGenesMicroglia_Pearson_cor0p2_abs/',
#                                                              '../../GWAS/PASCAL_results_non_neurological/Pearson_Cor0p2_Pathway/AllGenesMicroglia_Pearson_cor0p2_abs/'),
#                                                   c('../../genesets/Pearson_Cor0p2/CodingGenesMicroglia_Pearson_cor0p2_abs.gmt',
#                                                     '../../genesets/Pearson_Cor0p2/AllGenesMicroglia_Pearson_cor0p2_abs.gmt'),sizeLower = 3,sizeUpper = 300,appendInfo = T,empirical = T)
# 
# save(JoinedDfNonNeurologicalPearson,file='../../Count_Data/PASCAL_Results/JoinedDfNonNeurologicalPearsonEmp.rda')
# 
# sigClustersNonNeurologicalPearsonCoding <- GetAnnotationForSignificantClusters(JoinedDfNonNeurologicalPearson %>% dplyr::filter(Biotype=='coding'),isPath = F)
# save(sigClustersNonNeurologicalPearsonCoding,file='../../Count_Data/PASCAL_Results/sigClusters/sigClustersNonNeurologicalPearsonCoding.rda')
# 
# sigClustersNonNeurologicalPearsonAll<- GetAnnotationForSignificantClusters(JoinedDfNonNeurologicalPearson %>% dplyr::filter(Biotype=='all'),isPath = F)
# save(sigClustersNonNeurologicalPearsonAll,file='../../Count_Data/PASCAL_Results/sigClusters/sigClustersNonNeurologicalPearsonAll.rda')

# JoinedDfNonNeurologicalWGCNA <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results_non_neurological/WGCNA_size3/CodingWGCNAUnsigned_Soft4_Size3/',
#                                                         '../../GWAS/PASCAL_results_non_neurological/WGCNA_size3/AllWGCNAUnsigned_Soft4_Size3/'),
#                                        c('../../genesets/WGCNA_size3//CodingWGCNAUnsigned_Soft4_Size3.gmt',
#                                          '../../genesets/WGCNA_size3/AllWGCNAUnsigned_Soft4_Size3.gmt'),sizeLower = 3,sizeUpper = 300,appendInfo = T,empirical = T)
# save(JoinedDfNonNeurologicalWGCNA,file='../../Count_Data/PASCAL_Results/JoinedDfNonNeurologicalWGCNAEmp.rda')
# sigClustersNonNeurologicalWGCNACoding <- GetAnnotationForSignificantClusters(JoinedDfNonNeurologicalWGCNA %>% dplyr::filter(Biotype=='coding'),isPath = F)
# save(sigClustersNonNeurologicalWGCNACoding,file='../../Count_Data/PASCAL_Results/sigClusters/sigClustersNonNeurologicalWGCNACoding.rda')
# sigClustersNonNeurologicalWGCNAAll <- GetAnnotationForSignificantClusters(JoinedDfNonNeurologicalWGCNA %>% dplyr::filter(Biotype=='all'),isPath = F)
# save(sigClustersNonNeurologicalWGCNAAll,file='../../Count_Data/PASCAL_Results/sigClusters/sigClustersNonNeurologicalWGCNAAll.rda')


# GetGWASMetdata(sigClustersNonNeurologicalWGCNA)
# 

###################PEARSON TRANSCRIPT################################
# JoinedDfMicrogliaPearsonTranscripts <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results/Pearson_Cor0p2_Transcripts/CodingTranscriptsMicroglia_Pearson_cor0p2_abs/',
#                                                                     '../../GWAS/PASCAL_results/Pearson_Cor0p2_Transcripts/AllTranscriptsMicroglia_Pearson_cor0p2_abs/'),
#                                               c('../../genesets/Pearson_Cor0p2_Transcripts/CodingTranscriptsMicroglia_Pearson_cor0p2_abs.gmt',
#                                                 '../../genesets/Pearson_Cor0p2_Transcripts/AllTranscriptsMicroglia_Pearson_cor0p2_abs.gmt'),sizeLower = 3,sizeUpper = 300,appendInfo = F)
# save(JoinedDfMicrogliaPearsonTranscripts,file='../../Count_Data/PASCAL_Results/JoinedDfMicrogliaPearsonTranscripts.rda')
# codingGenesInNetworkTranscripts <- unique(unlist(JoinedDfMicrogliaPearsonTranscripts %>% dplyr::filter(Biotype=='coding') %>% {.$GeneNames}))
# allGenesInNetworkTranscripts <- unique(unlist(JoinedDfMicrogliaPearsonTranscripts %>% dplyr::filter(Biotype=='all') %>% {.$GeneNames}))
# sigClustersCodingPearsonTranscripts <- JoinedDfMicrogliaPearsonTranscripts %>% dplyr::filter(adjPvalue<0.05 & Biotype=='coding') %>%
# {AppendOpenTarget(JoinedDfMicroglia = .,
#                   csvPath = '/local/data/public/zmx21/zmx21_private/GSK/OpenTargets_scores/',
#                   codingGenesInNetwork=codingGenesInNetworkTranscripts,
#                   allGenesInNetwork=allGenesInNetworkTranscripts)} %>%
#                   {AppendGeneEnrichment(JoinedDfMicroglia = .,
#                                         codingGenesInNetwork=codingGenesInNetworkTranscripts,
#                                         allGenesInNetwork=allGenesInNetworkTranscripts)}
# 
# sigClustersCodingPearsonTranscripts <- GetAnnotationForSignificantClusters(sigClustersCodingPearsonTranscripts,isPath = F)
# save(sigClustersCodingPearsonTranscripts,file='../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingPearsonTranscripts.rda')
# 
# sigClustersAllPearsonTranscripts <- JoinedDfMicrogliaPearsonTranscripts %>% dplyr::filter(adjPvalue<0.05 & Biotype=='all') %>%
# {AppendOpenTarget(JoinedDfMicroglia = .,
#                   csvPath = '/local/data/public/zmx21/zmx21_private/GSK/OpenTargets_scores/',
#                   codingGenesInNetwork=codingGenesInNetworkTranscripts,
#                   allGenesInNetwork=allGenesInNetworkTranscripts)} %>%
#                   {AppendGeneEnrichment(JoinedDfMicroglia = .,
#                                         codingGenesInNetwork=codingGenesInNetworkTranscripts,
#                                         allGenesInNetwork=allGenesInNetworkTranscripts)}
# 
# sigClustersAllPearsonTranscripts <- GetAnnotationForSignificantClusters(sigClustersAllPearsonTranscripts,isPath = F)
# save(sigClustersAllPearsonTranscripts,file='../../Count_Data/PASCAL_Results/sigClusters/sigClustersAllPearsonTranscripts.rda')

# ViewOpenTargetResults(sigClustersCodingPearsonTranscripts,type = 'target')
# ViewOpenTargetResults(sigClustersCodingPearsonTranscripts,type = 'drug')

###################WGNCA TRANSCRIPT################################
# JoinedDfMicrogliaWGCNATranscripts <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results/WGCNA_size3_transcripts/CodingWGCNAUnsigned_Soft4_Size3_DeepSplit2/'),
#                                                          c('../../GWAS/PASCAL_New/resources/genesets/WGCNA_size3_transcripts/CodingWGCNAUnsigned_Soft4_Size3_DeepSplit2.gmt'),sizeLower = 3,sizeUpper = 300,appendInfo = F)
# save(JoinedDfMicrogliaWGCNATranscripts,file='../../Count_Data/PASCAL_Results/JoinedDfMicrogliaWGCNATranscripts.rda')
# codingGenesInNetworkTranscripts <- unique(unlist(JoinedDfMicrogliaWGCNATranscripts %>% dplyr::filter(Biotype=='coding') %>% {.$GeneNames}))
# allGenesInNetworkTranscripts <- unique(unlist(JoinedDfMicrogliaWGCNATranscripts %>% dplyr::filter(Biotype=='all') %>% {.$GeneNames}))
# 
# sigClustersCodingWGCNATranscripts <- JoinedDfMicrogliaWGCNATranscripts %>% dplyr::filter(adjPvalue<0.05) %>%
# {AppendOpenTarget(JoinedDfMicroglia = .,
#                   csvPath = '/local/data/public/zmx21/zmx21_private/GSK/OpenTargets_scores/',
#                   codingGenesInNetwork=codingGenesInNetworkTranscripts,
#                   allGenesInNetwork=allGenesInNetworkTranscripts)} %>% 
#                   {AppendGeneEnrichment(JoinedDfMicroglia = .,
#                                         codingGenesInNetwork=codingGenesInNetworkTranscripts,
#                                         allGenesInNetwork=allGenesInNetworkTranscripts)}
# 
# sigClustersCodingWGCNATranscripts <- GetAnnotationForSignificantClusters(sigClustersCodingWGCNATranscripts,isPath = F)
# save(sigClustersCodingWGCNATranscripts,file='../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingWGCNATranscripts.rda')

# ViewOpenTargetResults(sigClustersCodingWGCNATranscripts,type = 'target')
# ViewOpenTargetResults(sigClustersCodingWGCNATranscripts,type = 'drug')


#############################################OPEN TARGET PEARSON#########################################
# load('../../Count_Data/PASCAL_Results/Microglia_Pearson_cor0p2_abs.rda')


# JoinedDfMicrogliaPearson <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results/Pearson_Cor0p2/AllGenesMicroglia_Pearson_cor0p2_abs/'),
#                                         c('../../Louvain_results/Pearson_Cor0p2/AllGenesMicroglia_Pearson_cor0p2_abs.gmt'),sizeLower = 3,sizeUpper = 300,empirical = T)
# save(JoinedDfMicrogliaPearson,file='../../Count_Data/PASCAL_Results/Microglia_Pearson_cor0p2_absAllEmp.rda')


# PearsonOpenTargetSig <- JoinedDfMicrogliaPearson %>%
#   dplyr::filter(adjPvalue > 0.05) %>%
#   dplyr::filter(AD_target_P < 0.05,ALS_target_P < 0.05 ,MS_target_P < 0.05 ,PD_target_P < 0.05) %>%
#   dplyr::arrange(AD_target_P,ALS_target_P,MS_target_P,PD_target_P) %>% dplyr::distinct(Name,.keep_all=T) %>%
#   {GetAnnotationForSignificantClusters(.,isPath = F,pCutOff = 1)}
# for(i in 1:nrow(PearsonOpenTargetSig)){
#   PearsonOpenTargetSig$StudyName[i] <- dplyr::filter(JoinedDfMicrogliaPearson,Name==PearsonOpenTargetSig$Name[i]) %>%
#     dplyr::arrange(adjPvalue) %>% {.$StudyName[1]}
#   PearsonOpenTargetSig$adjPvalue[i] <- dplyr::filter(JoinedDfMicrogliaPearson,Name==PearsonOpenTargetSig$Name[i]) %>%
#     dplyr::arrange(adjPvalue) %>% {.$adjPvalue[1]}
# }
# PearsonOpenTargetSig <- PearsonOpenTargetSig %>% dplyr::filter(adjPvalue > 0.05)
# ViewOpenTargetResults(PearsonOpenTargetSig %>% dplyr::filter(Biotype=='coding'),type = 'target',empirical = T)
# ViewOpenTargetResults(PearsonOpenTargetSig %>% dplyr::filter(Biotype=='all'),type = 'target',empirical = T)
#############################################OPEN TARGET WGCNA#########################################

# load('../../Count_Data/PASCAL_Results/JoinedDfMicrogliaWGCNAUnsigned.rda')

# load('../../Count_Data/PASCAL_Results/JoinedDfMicrogliaWGCNAUnsigned_AllEmp.rda')
# WGCNAOpenTargetSig <- JoinedDfMicrogliaWGCNAUnsigned %>%
#   dplyr::filter(adjPvalue > 0.05) %>%
#   dplyr::filter(AD_target_P < 0.05,ALS_target_P < 0.05 ,MS_target_P < 0.05 ,PD_target_P < 0.05) %>%
#   dplyr::arrange(AD_target_P,ALS_target_P,MS_target_P,PD_target_P) %>% dplyr::distinct(Name,.keep_all=T) %>%
#   {GetAnnotationForSignificantClusters(.,isPath = F,pCutOff = 1)}
# for(i in 1:nrow(WGCNAOpenTargetSig)){
#   WGCNAOpenTargetSig$StudyName[i] <- dplyr::filter(JoinedDfMicrogliaWGCNAUnsigned,Name==WGCNAOpenTargetSig$Name[i]) %>%
#     dplyr::arrange(adjPvalue) %>% {.$StudyName[1]}
#   WGCNAOpenTargetSig$adjPvalue[i] <- dplyr::filter(JoinedDfMicrogliaWGCNAUnsigned,Name==WGCNAOpenTargetSig$Name[i]) %>%
#     dplyr::arrange(adjPvalue) %>% {.$adjPvalue[1]}
# }
# WGCNAOpenTargetSig <- WGCNAOpenTargetSig %>% dplyr::filter(adjPvalue > 0.05)
# ViewOpenTargetResults(WGCNAOpenTargetSig %>% dplyr::filter(Biotype=='coding'),type = 'target')
# ViewOpenTargetResults(WGCNAOpenTargetSig %>% dplyr::filter(Biotype=='all'),type = 'target',empirical = T)


# JoinedDfMicrogliaPearson <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results/Pearson_Cor0p2/CodingGenesMicroglia_Pearson_cor0p2_abs/',
#                                                   '../../GWAS/PASCAL_results/Pearson_Cor0p2/AllGenesMicroglia_Pearson_cor0p2_abs/'),
#                                        c('../../genesets/Pearson_Cor0p2/CodingGenesMicroglia_Pearson_cor0p2_abs.gmt',
#                                          '../../genesets/Pearson_Cor0p2/Pearson_Cor0p2/AllGenesMicroglia_Pearson_cor0p2_abs.gmt'),sizeLower = 3,sizeUpper = 300,appendInfo = T,empirical = T)
# 
# saveRDS(JoinedDfMicrogliaPearson,file='../../Count_Data/PASCAL_Results/Empricial_Df/Microglia_Pearson_cor0p2_abs.rds')
# 
# JoinedDfMicrogliaWGCNAUnsigned <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results/WGCNA_size3/CodingWGCNAUnsigned_Soft4_Size3//',
#                                                              '../../GWAS/PASCAL_results/WGCNA_size3/AllWGCNAUnsigned_Soft4_Size3//'),
#                                                   c('../../genesets/WGCNA_size3/CodingWGCNAUnsigned_Soft4_Size3.gmt',
#                                                     '../../genesets/WGCNA_size3/AllWGCNAUnsigned_Soft4_Size3.gmt'),sizeLower = 3,sizeUpper = 100,appendInfo = T,empirical = T)
# saveRDS(JoinedDfMicrogliaWGCNAUnsigned,file='../../Count_Data/PASCAL_Results/Empricial_Df/JoinedDfMicrogliaWGCNAUnsigned.rds')
# 
# JoinedDfNonNeurologicalPearson <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results_non_neurological/Pearson_Cor0p2_Pathway/CodingGenesMicroglia_Pearson_cor0p2_abs/',
#                                                              '../../GWAS/PASCAL_results_non_neurological/Pearson_Cor0p2_Pathway/AllGenesMicroglia_Pearson_cor0p2_abs/'),
#                                                   c('../../genesets/Pearson_Cor0p2/CodingGenesMicroglia_Pearson_cor0p2_abs.gmt',
#                                                     '../../genesets/Pearson_Cor0p2/AllGenesMicroglia_Pearson_cor0p2_abs.gmt'),sizeLower = 3,sizeUpper = 300,appendInfo = T,empirical = T)
# 
# saveRDS(JoinedDfNonNeurologicalPearson,file='../../Count_Data/PASCAL_Results/Empricial_Df/JoinedDfNonNeurologicalPearson.rds')
# 
# JoinedDfNonNeurologicalWGCNA <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results_non_neurological/WGCNA_size3/CodingWGCNAUnsigned_Soft4_Size3/',
#                                                         '../../GWAS/PASCAL_results_non_neurological/WGCNA_size3/AllWGCNAUnsigned_Soft4_Size3/'),
#                                        c('../../genesets/WGCNA_size3//CodingWGCNAUnsigned_Soft4_Size3.gmt',
#                                          '../../genesets/WGCNA_size3/AllWGCNAUnsigned_Soft4_Size3.gmt'),sizeLower = 3,sizeUpper = 300,appendInfo = T,empirical = T)
# saveRDS(JoinedDfNonNeurologicalWGCNA,file='../../Count_Data/PASCAL_Results/Empricial_Df/JoinedDfNonNeurologicalWGCNA.rds')
# 
# JoinedDfMicrogliaPearsonTranscripts <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results/Pearson_Cor0p2_Transcripts/CodingTranscriptsMicroglia_Pearson_cor0p2_abs/',
#                                                                     '../../GWAS/PASCAL_results/Pearson_Cor0p2_Transcripts/AllTranscriptsMicroglia_Pearson_cor0p2_abs/'),
#                                               c('../../genesets/Pearson_Cor0p2_Transcripts/CodingTranscriptsMicroglia_Pearson_cor0p2_abs.gmt',
#                                                 '../../genesets/Pearson_Cor0p2_Transcripts/AllTranscriptsMicroglia_Pearson_cor0p2_abs.gmt'),sizeLower = 3,sizeUpper = 300,appendInfo = T,empirical = T)
# saveRDS(JoinedDfMicrogliaPearsonTranscripts,file='../../Count_Data/PASCAL_Results/Empricial_Df/JoinedDfMicrogliaPearsonTranscripts.rds')
# 
# JoinedDfMicrogliaWGCNATranscripts <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results/WGCNA_size3_transcripts/CodingWGCNAUnsigned_Soft4_Size3_DeepSplit2/'),
#                                                          c('../../genesets/WGCNA_size3_transcripts/CodingWGCNAUnsigned_Soft4_Size3_DeepSplit2.gmt'),sizeLower = 3,sizeUpper = 300,appendInfo = T,empirical = T)
# saveRDS(JoinedDfMicrogliaWGCNATranscripts,file='../../Count_Data/PASCAL_Results/Empricial_Df/JoinedDfMicrogliaWGCNATranscripts.rds')
