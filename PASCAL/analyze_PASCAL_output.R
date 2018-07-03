library(qusage)
library(dplyr)

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
  
  #For each cluster, load it's respective gene names.
  clustersInfo <-  unlist(lapply(clusterPaths,function(x) qusage::read.gmt(file=x)),recursive = F)
  clustersDf <- data_frame(ClusterName = names(clustersInfo),Genes = clustersInfo)
  
  #Join result and gene info for each cluster. Filter for those with P-values. Then, extract metainfo from cluster Name (level,type)
  joinedDf <- dplyr::left_join(resultDf,clustersDf,by=c('Name'='ClusterName')) %>% dplyr::filter(!is.na(chi2Pvalue)) %>% 
    dplyr::mutate(Size=sapply(Genes,length)) %>%  #Size info
    dplyr::mutate(Biotype=factor(sapply(Name,function(x) unlist(strsplit(x,'_'))[1])))  %>% #Biotype
    dplyr::mutate(CellType=factor(sapply(Name,function(x) unlist(strsplit(x,'_'))[2])))  %>% #Cell Type
    dplyr::mutate(AggregationType=factor(sapply(Name,function(x) unlist(strsplit(x,'_'))[3]))) %>% #Gene or transcript
    dplyr::mutate(Level=sapply(Name,function(x){
      splitString <- unlist(strsplit(x,'_'));
      return(as.numeric(splitString[which(splitString=='level')+1]))})) #Gene or transcript
  return(joinedDf)
}

#Get path to all studies in a directory
GetPathToPASCALResults <- function(PASCALResultPath){
  allStudies <- dir(PASCALResultPath)
  allStudies <- allStudies[grep('PathwaySet',allStudies)]
  allStudiesPath <- sapply(allStudies,function(x) paste0(PASCALResultPath,x))
  return(allStudiesPath)
}

allStudiesMicrogliaPath <- GetPathToPASCALResults('../../GWAS/PASCAL_results/microglia_gene/')
JoinedDfMicroglia <- ParsePASCALFile(allStudiesMicrogliaPath,'../../Louvain_results/microglia_gene_clusters.gmt')
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
PlotDensityByLevel <- function(joinedDf,title){
  return(ggplot(joinedDf,aes(x=-1*log10(empPvalue),fill=factor(Level))) + geom_density(trim=F) + 
    coord_flip() + 
    facet_grid(. ~ Level) + 
    xlim(1,5) + labs(x='-log10(Empirical P-Value)',y='Density',fill='Louvain\nRecursion\nLevel') + ggtitle(title))
} 

#Look at coding only Microglia distribution across levels, for all studies
p1 <- PlotDensityByLevel(JoinedDfMicroglia %>% dplyr::filter(Biotype=='coding' & Level >= 3 & Level < 11),title = c('Microglia Coding Genes-\nP-val Density Across Levels for all Studies'))
#Look at coding and non coding Microglia distribution across levels, for all studies
p2 <- PlotDensityByLevel(JoinedDfMicroglia %>% dplyr::filter(Biotype=='all' & Level >= 3 & Level < 11),title = c('Microglia All Genes -\nP-val Density Across Levels for all Studies'))
tiff(filename = '../../Figures/PASCAL/LevelVSPVal.tiff',width = 1200,height = 600)
ggarrange(p1,p2,nrow=2)
dev.off()


#Plot density for inidividual studies
library(scales)
allStudiescoding <- lapply(1:length(unique(JoinedDfMicroglia$StudyName)),function(i) PlotDensityByLevel(JoinedDfMicroglia %>% 
                           dplyr::filter(Biotype=='coding' & (Level %in% c(4,5,6,7)) & StudyName == unique(JoinedDfMicroglia$StudyName)[i]),
                         title = paste0('Microglia Coding Genes -\n',unique(JoinedDfMicroglia$StudyName)[i])) + scale_y_continuous(breaks = pretty_breaks(n = 3)))
tiff(filename = '../../Figures/PASCAL/CodingPValIndivStudies.tiff',width = 1600,height = 600)
ggarrange(plots = allStudiescoding,ncol = 4,nrow = 2)
dev.off()

allStudiesallgenes <- lapply(1:length(unique(JoinedDfMicroglia$StudyName)),function(i) PlotDensityByLevel(JoinedDfMicroglia %>% 
                                                                                                          dplyr::filter(Biotype=='all' & (Level %in% c(4,5,6,7)) & StudyName == unique(JoinedDfMicroglia$StudyName)[i]),
                                                                                                        title = paste0('Microglia All Genes -\n',unique(JoinedDfMicroglia$StudyName)[i])) + scale_y_continuous(breaks = pretty_breaks(n = 3)))
tiff(filename = '../../Figures/PASCAL/AllPValIndivStudies.tiff',width = 1600,height = 600)
ggarrange(plots = allStudiesallgenes,ncol = 4,nrow = 2)
dev.off()


#Plot density vs random cluster
allStudiesRandomMicrogliaPath <- GetPathToPASCALResults('../../GWAS/PASCAL_results/randommicroglia_gene/')
JoinedDfMicrogliaRandom <- ParsePASCALFile(allStudiesRandomMicrogliaPath,'../../Louvain_results/randommicroglia_clusters.gmt')
tiff(filename = '../../Figures/PASCAL/RandomVSTrueNetwork.tiff',width = 1200,height = 400)
ggplot(rbind(JoinedDfMicrogliaRandom,JoinedDfMicroglia) %>% 
         dplyr::filter(Biotype!='randomall' & Biotype!='all' & Level > 3) %>% 
         dplyr::mutate(Type=ifelse(Biotype=='randomcoding','Random','True')), 
       aes(y=-1*log10(empPvalue),fill=factor(Type),x=factor(Level))) + geom_split_violin() + 
  ylim(2.3,5) + xlab('Louvain Recursion Level') + ylab('-log10(Empirical P Value)') +
  labs(fill='Random \n or \nTrue') + ggtitle('Randomly permuted network VS True Network')
dev.off()

