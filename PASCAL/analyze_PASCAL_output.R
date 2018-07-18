library(qusage)
library(dplyr)
library(ggplot2)
library(parallel)
source('annotate_clusters.R')
source('pval_correction.R')
source('write_GWAS_metadata.R')
#Merged result from PASCAL with cluster info
ParsePASCALFile <- function(resultPaths,clusterPaths){
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
  joinedDf$GeneNames <- mclapply(joinedDf$Genes,function(x) geneIdToName[[x]],mc.cores = 20)
  source('append_gene_enrichment.R')
  joinedDf <- AppendGeneEnrichment(joinedDf)
  source('append_open_target.R')
  joinedDf <- AppendOpenTarget(joinedDf,csvPath = '/local/data/public/zmx21/zmx21_private/GSK/OpenTargets_scores/')
  # 
  # famMonoDiseases <- c('familial_ALS','familial_MS','monogenic_AD','monogenic_recessive_PD')
  # famMonoGenes <- lapply(famMonoDiseases,function(x) data.table::fread(paste0('/local/data/public/zmx21/zmx21_private/GSK/Monogenic_disease/',x,'.txt')) %>% {.$V1})
  # 
  return(joinedDf)
}

#Get path to all studies in a directory
GetPathToPASCALResults <- function(PASCALResultPath){
  allStudies <- dir(PASCALResultPath)
  allStudies <- allStudies[grep('PathwaySet',allStudies)]
  allStudiesPath <- sapply(allStudies,function(x) paste0(PASCALResultPath,x))
  return(allStudiesPath)
}


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

LoadPASCALResults <- function(filter=T,resultPaths,clusterPaths,sizeLower=10,sizeUpper=200){
  source('pval_correction.R')
  microgliaAllGenesPath <- GetPathToPASCALResults(resultPaths[1])
  JoinedDfMicroglia <- ParsePASCALFile(microgliaAllGenesPath,clusterPaths[1])
  if(length(resultPaths)== 2 & length(clusterPaths) == 2){
    microgliaCodingGenesPath <- GetPathToPASCALResults(resultPaths[2])
    JoinedDfMicroglia <- rbind(JoinedDfMicroglia,ParsePASCALFile(microgliaCodingGenesPath,clusterPaths[2]))
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

ViewOpenTargetResults <- function(JoinedDfMicroglia,type){
  round_df <- function(x, digits) {
    # round all numeric variables
    # x: data frame 
    # digits: number of digits to round
    numeric_columns <- sapply(x, mode) == 'numeric'
    x[numeric_columns] <-  round(x[numeric_columns], digits)
    return(x)
  }
  
  if(type=='target'){
    df <- JoinedDfMicroglia %>% dplyr::select(StudyName,Size,adjPvalue,KEGG=KEGG_2016.Term,'KEGG_Overlap'=KEGG_2016.Overlap,'KEGG_P'=KEGG_2016.Adjusted.P.value,
                                              contains('target'),-contains('target_Genes'))
    colnames(df) <- sapply(colnames(df),function(x) gsub(pattern = '_',replacement = '\n',x = x))
    df <- round_df(df,2)
    df <- dplyr::arrange(df,KEGG)
    rownames(df) <- c()
  }else if(type=='drug'){
    df <- JoinedDfMicroglia %>% dplyr::select(StudyName,Size,adjPvalue,KEGG=KEGG_2016.Term,'KEGG_Overlap'=KEGG_2016.Overlap,'KEGG_P'=KEGG_2016.Adjusted.P.value,
                                              contains('drug'),-contains('drug_Genes'))
    colnames(df) <- sapply(colnames(df),function(x) gsub(pattern = '_',replacement = '\n',x = x))
    df <- round_df(df,2)
    df <- dplyr::arrange(df,KEGG)
    rownames(df) <- c()
  }
  grid.arrange(tableGrob(df,theme = ttheme_default(base_size = 9)))
}

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


JoinedDfMicrogliaPearson <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results/Pearson_Cor0p2/CodingGenesMicroglia_Pearson_cor0p2_abs/',
                                                  '../../GWAS/PASCAL_results/Pearson_Cor0p2/AllGenesMicroglia_Pearson_cor0p2_abs/'),
                                       c('../../Louvain_results/Pearson_Cor0p2/CodingGenesMicroglia_Pearson_cor0p2_abs.gmt',
                                         '../../Louvain_results/Pearson_Cor0p2/AllGenesMicroglia_Pearson_cor0p2_abs.gmt'),sizeLower = 3,sizeUpper = 300)
save(JoinedDfMicrogliaPearson,file='../../Count_Data/PASCAL_Results/Microglia_Pearson_cor0p2_abs.rda')
load('../../Count_Data/PASCAL_Results/Microglia_Pearson_cor0p2_abs.rda')
sigClustersCodingPearson <- GetAnnotationForSignificantClusters(JoinedDfMicrogliaPearson %>% dplyr::filter(Biotype=='coding'),isPath = F)
save(sigClustersCodingPearson,file='../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingPearson.rda')
ViewOpenTargetResults(sigClustersCodingPearson,type='target')
ViewOpenTargetResults(sigClustersCodingPearson,type='drug')
sigClustersAllPearson <- GetAnnotationForSignificantClusters(JoinedDfMicrogliaPearson %>% dplyr::filter(Biotype=='all'),isPath = F)
save(sigClustersAllPearson,file='../../Count_Data/PASCAL_Results/sigClusters/sigClustersAllPearson.rda')
ViewOpenTargetResults(sigClustersAllPearson,type='target')
ViewOpenTargetResults(sigClustersAllPearson,type='drug')
GetGWASMetdata(rbind(sigClustersCodingPearson,sigClustersAllPearson))



# JoinedDfMicrogliaWGCNASigned <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results/WGCNA_size3/CodingWGCNASigned_Soft4_Size3//',
#                                                          '../../GWAS/PASCAL_results/WGCNA_size3/AllWGCNASigned_Soft4_Size3//'),
#                                               c('../../GWAS/PASCAL_New/resources/genesets/WGCNA_size3/CodingWGCNASigned_Soft4_Size3.gmt',
#                                                 '../../GWAS/PASCAL_New/resources/genesets/WGCNA_size3/AllWGCNASigned_Soft4_Size3.gmt'),sizeLower = 3,sizeUpper = 100)
# sigClustersCodingWGCNASigned <- do.call(rbind,lapply(GetAnnotationForSignificantClusters(JoinedDfMicrogliaWGCNASigned %>% dplyr::filter(Biotype=='coding')
#                                                                                      ,isPath = F),function(x) x$df))
# sigClustersAllWGCNASigned <- do.call(rbind,lapply(GetAnnotationForSignificantClusters(JoinedDfMicrogliaWGCNASigned %>% dplyr::filter(Biotype=='all')
#                                                                                   ,isPath = F),function(x) x$df))
# GetGWASMetdata(rbind(sigClustersCodingWGCNASigned,sigClustersCodingWGCNASigned))


JoinedDfMicrogliaWGCNAUnsigned <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results/WGCNA_size3/CodingWGCNAUnsigned_Soft4_Size3//',
                                                             '../../GWAS/PASCAL_results/WGCNA_size3/AllWGCNAUnsigned_Soft4_Size3//'),
                                                  c('../../GWAS/PASCAL_New/resources/genesets/WGCNA_size3/CodingWGCNAUnsigned_Soft4_Size3.gmt',
                                                    '../../GWAS/PASCAL_New/resources/genesets/WGCNA_size3/AllWGCNAUnsigned_Soft4_Size3.gmt'),sizeLower = 3,sizeUpper = 100)
save(JoinedDfMicrogliaWGCNAUnsigned,file='../../Count_Data/PASCAL_Results/JoinedDfMicrogliaWGCNAUnsigned.rda')

sigClustersCodingWGCNAUnsigned <- GetAnnotationForSignificantClusters(JoinedDfMicrogliaWGCNAUnsigned %>% dplyr::filter(Biotype=='coding'),isPath = F)
save(sigClustersCodingWGCNAUnsigned,file='../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingWGCNAUnsigned.rda')

ViewOpenTargetResults(sigClustersCodingWGCNAUnsigned,type='target')
ViewOpenTargetResults(sigClustersCodingWGCNAUnsigned,type='drug')
sigClustersAllWGCNAUnsigned <- GetAnnotationForSignificantClusters(JoinedDfMicrogliaWGCNAUnsigned %>% dplyr::filter(Biotype=='all'),isPath = F)
save(sigClustersAllWGCNAUnsigned,file='../../Count_Data/PASCAL_Results/sigClusters/sigClustersAllWGCNAUnsigned.rda')

ViewOpenTargetResults(sigClustersAllWGCNAUnsigned,type='target')
ViewOpenTargetResults(sigClustersAllWGCNAUnsigned,type='drug')
GetGWASMetdata(rbind(sigClustersCodingWGCNAUnsigned,sigClustersCodingWGCNAUnsigned))

# JoinedDfMicroglia <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results/microglia_all_genes/','../../GWAS/PASCAL_results/microglia_coding_genes/'),
#                                        c('../../Louvain_results/AllMicrogliaGenes/microgliaAllGenesClusters.gmt','../../Louvain_results/CodingMicrogliaGenes/microgliaCodingGenesClusters.gmt'))
# 
# JoinedDfMicroglia <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results/microglia_all_genes_pval0p01/','../../GWAS/PASCAL_results/microglia_coding_genes_pval0p01/'),
#                                        c('../../Louvain_results/AllMicrogliaGenes_pval0p01/allMicrogliaGenesClusters_pval0p01.gmt','../../Louvain_results/CodingMicrogliaGenes_pval0p01/codingMicrogliaGenesClusters_pval0p01.gmt'))

# JoinedDfMicroglia <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results/microglia_coding_genes_signed/'),
#                                        c('../../Louvain_results/CodingMicrogliaGenesSigned/microgliaCodingGenesSigned.gmt'))



# JoinedDfMicroglia <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results2/CodingGenesMicroglia_Jaccard_pval0p05_cor0p25_abs/'),
#                                        c('../../Louvain_results/Jaccard/CodingGenesMicroglia_Jaccard_pval0p05_cor0p25_abs.gmt'),sizeLower = 5)
# save(JoinedDfMicroglia,file='../../Count_Data/PASCAL_Results/CodingGenesMicroglia_Jaccard_pval0p05_cor0p25_abs.rda')
# sigClusters <- do.call(rbind,lapply(GetAnnotationForSignificantClusters(JoinedDfMicroglia,isPath = F),function(x) x$df))
# GetGWASMetdata(sigClusters)
# 
# JoinedDfMicroglia <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results2/CodingGenesMicroglia_Jaccard_pval0p05_cor0p25_noabs/'),
#                                        c('../../Louvain_results/Jaccard/CodingGenesMicroglia_Jaccard_pval0p05_cor0p25_noabs.gmt'))
# save(JoinedDfMicroglia,file='../../Count_Data/PASCAL_Results/CodingGenesMicroglia_Jaccard_pval0p05_cor0p25_noabs.rda')
# sigClusters <- do.call(rbind,lapply(GetAnnotationForSignificantClusters('../../Count_Data/PASCAL_Results/CodingGenesMicroglia_Jaccard_pval0p05_cor0p25_noabs.rda'),function(x) x$df))
# GetGWASMetdata(sigClusters)
# 
# 
# JoinedDfMicroglia <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results2/CodingGenesMicroglia_Jaccard_top1milpos/'),
#                                        c('../../Louvain_results/Jaccard/CodingGenesMicroglia_Jaccard_top1milpos.gmt'))
# save(JoinedDfMicroglia,file='../../Count_Data/PASCAL_Results/CodingGenesMicroglia_Jaccard_top1milpos.rda')
# sigClusters <- do.call(rbind,lapply(GetAnnotationForSignificantClusters('../../Count_Data/PASCAL_Results/CodingGenesMicroglia_Jaccard_top1milpos.rda'),function(x) x$df))
# GetGWASMetdata(sigClusters)


# JoinedDfMicroglia <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results2/CodingGenesMicroglia_Pearson_pval0p05_cor0p25_abs/'),
#                                        c('../../Louvain_results/Pearson/CodingGenesMicroglia_Pearson_pval0p05_cor0p25_abs.gmt'))
# save(JoinedDfMicroglia,file='../../Count_Data/PASCAL_Results/CodingGenesMicroglia_Pearson_pval0p05_cor0p25_abs.rda')
# sigClusters <- do.call(rbind,lapply(GetAnnotationForSignificantClusters('../../Count_Data/PASCAL_Results/CodingGenesMicroglia_Pearson_pval0p05_cor0p25_abs.rda'),function(x) x$df))
# GetGWASMetdata(sigClusters)
# 
# 
# JoinedDfMicroglia <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results2/CodingGenesMicroglia_Pearson_pval0p05_cor0p25_noabs/'),
#                                        c('../../Louvain_results/Pearson/CodingGenesMicroglia_Pearson_pval0p05_cor0p25_noabs.gmt'))
# save(JoinedDfMicroglia,file='../../Count_Data/PASCAL_Results/CodingGenesMicroglia_Pearson_pval0p05_cor0p25_noabs.rda')
# sigClusters <- do.call(rbind,lapply(GetAnnotationForSignificantClusters('../../Count_Data/PASCAL_Results/CodingGenesMicroglia_Pearson_pval0p05_cor0p25_noabs.rda'),function(x) x$df))
# GetGWASMetdata(sigClusters)
# 
# 
# JoinedDfMicroglia <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results2/CodingGenesMicroglia_Pearson_top1milpos/'),
#                                        c('../../Louvain_results/Pearson/CodingGenesMicroglia_Pearson_top1milpos.gmt'))
# save(JoinedDfMicroglia,file='../../Count_Data/PASCAL_Results/CodingGenesMicroglia_Pearson_top1milpos.rda')
# sigClusters <- do.call(rbind,lapply(GetAnnotationForSignificantClusters('../../Count_Data/PASCAL_Results/CodingGenesMicroglia_Pearson_top1milpos.rda'),function(x) x$df))
# GetGWASMetdata(sigClusters)



# JoinedDfMicroglia <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results_WGCNA/CodingWGCNASigned/'),
#                                        c('../../GWAS/PASCAL_New/resources/genesets/WGCNA_clusters/CodingWGCNASigned.gmt'))
# save(JoinedDfMicroglia,file='../../Count_Data/PASCAL_Results/WGCNASigned.rda')
# sigClusters <- do.call(rbind,lapply(GetAnnotationForSignificantClusters('../../Count_Data/PASCAL_Results/WGCNASigned.rda'),function(x) x$df))
# GetGWASMetdata(sigClusters)
# 
# 
# JoinedDfMicroglia <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results_WGCNA/CodingWGCNAUnsigned/'),
#                                        c('../../GWAS/PASCAL_New/resources/genesets/WGCNA_clusters/CodingWGCNAUnsigned.gmt'))
# save(JoinedDfMicroglia,file='../../Count_Data/PASCAL_Results/WGCNAUnsigned.rda')
# sigClusters <- do.call(rbind,lapply(GetAnnotationForSignificantClusters('../../Count_Data/PASCAL_Results/WGCNAUnsigned.rda'),function(x) x$df))
# GetGWASMetdata(sigClusters)
# 

# JoinedDfMicroglia <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results_WGCNA_soft4/CodingWGCNASigned_Soft4/'),
#                                        c('../../GWAS/PASCAL_New/resources/genesets/WGCNA_clusters_soft4/CodingWGCNASigned_Soft4.gmt'))
# save(JoinedDfMicroglia,file='../../Count_Data/PASCAL_Results/WGCNASigned_Soft4.rda')
# sigClusters <- do.call(rbind,lapply(GetAnnotationForSignificantClusters(JoinedDfMicroglia,isPath = F),function(x) x$df))
# GetGWASMetdata(sigClusters)


JoinedDfNonNeurologicalPearson <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results_non_neurological/Pearson_Cor0p2_Pathway/CodingGenesMicroglia_Pearson_cor0p2_abs/',
                                                             '../../GWAS/PASCAL_results_non_neurological/Pearson_Cor0p2_Pathway/AllGenesMicroglia_Pearson_cor0p2_abs/'),
                                                  c('../../GWAS/PASCAL_New/resources/genesets/Pearson_Cor0p2/CodingGenesMicroglia_Pearson_cor0p2_abs.gmt',
                                                    '../../GWAS/PASCAL_New/resources/genesets/Pearson_Cor0p2/AllGenesMicroglia_Pearson_cor0p2_abs.gmt'))
sigClustersNonNeurologicalPearson<- do.call(rbind,lapply(GetAnnotationForSignificantClusters(JoinedDfNonNeurologicalPearson,isPath = F),function(x) x$df))

JoinedDfNonNeurologicalWGCNA <- LoadPASCALResults(filter=T,c('../../GWAS/PASCAL_results_non_neurological/WGCNA_size3/CodingWGCNAUnsigned_Soft4_Size3/',
                                                        '../../GWAS/PASCAL_results_non_neurological/WGCNA_size3/AllWGCNAUnsigned_Soft4_Size3/'),
                                       c('../../GWAS/PASCAL_New/resources/genesets/WGCNA_size3/CodingWGCNAUnsigned_Soft4_Size3.gmt',
                                         '../../GWAS/PASCAL_New/resources/genesets/WGCNA_size3/AllWGCNAUnsigned_Soft4_Size3.gmt'))
sigClustersNonNeurologicalWGCNA <- do.call(rbind,lapply(GetAnnotationForSignificantClusters(JoinedDfNonNeurologicalWGCNA,isPath = F),function(x) x$df))
GetGWASMetdata(sigClustersNonNeurologicalWGCNA)

