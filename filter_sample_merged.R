FilterCountMatrix <- function(countMatrixInput,meanAbundances,cv,cvCutOffAbsolute,meanCutOffAbsolute,mappingTable=NULL){
  #Filter for genes according to cutoff
  belowMeanCutOff <- which(meanAbundances <= meanCutOffAbsolute)
  belowCVCutoff <- which(cv <= cvCutOffAbsolute | is.nan(cv))
  rowsToFilter <- union(belowMeanCutOff,belowCVCutoff)
  
  filteredCountMatrix <- countMatrixInput[setdiff(1:nrow(countMatrixInput),rowsToFilter),]
  
  #Decide whether to return categorial counts.
  if(is.null(mappingTable)){
    return(list(countMatrix = filteredCountMatrix,numFiltered = length(rowsToFilter)))
  }else{
    type <- ifelse(grepl('ENST',rownames(countMatrixInput)[1]),'transcript_id','gene_id')
    allFilteredNames <- data.frame(ids = rownames(filteredCountMatrix),stringsAsFactors = F)
    colnames(allFilteredNames) <- type
    categoricalCounts <- as.data.frame(dplyr::left_join(allFilteredNames,mappingTable,by=c(type,type)))
    return(list(countMatrix = filteredCountMatrix,numFiltered = length(rowsToFilter),categoricalCounts = categoricalCounts))
  }
}

load(file = '../Count_Data/gtfTables.rda')

source('load_GTF.R')
# gtfTables <- LoadBiotypeMapping()


load('../Count_Data/SalmonTPM_Gene_Merged_ReadFilt.rda')
load('../Count_Data/SalmonTPM_Transcript_Merged_ReadFilt.rda')
library(dplyr)
source('batch_adjustment.R')
source('pca_analysis.R')

# #Don't consider Olah samples
# SalmonTPM_Gene_Merged <- SalmonTPM_Gene_Merged[,!sapply(colnames(SalmonTPM_Gene_Merged),function(x) grepl('H5KN',x))]
# SalmonTPM_Transcript_Merged <- SalmonTPM_Transcript_Merged[,!sapply(colnames(SalmonTPM_Transcript_Merged),function(x) grepl('H5KN',x))]

#Remove samples from Galatro which have bad mapping. 
BadMappingGalatro <- read.table(file='/local/data/public/zmx21/zmx21_private/GSK/Galatro/STAR_aligned_merged/multiqc_data/multiqc_rseqc_read_distribution.txt',
                      header = T,stringsAsFactors = F) %>% dplyr::select('Sample','cds_exons_tag_count','total_tags','other_intergenic_tag_count') %>% 
  dplyr::mutate(Sample_Name = sapply(.$Sample,function(x) unlist(strsplit(x,'[.]'))[1]),
                PercExonic = round(as.numeric(cds_exons_tag_count/total_tags),2),
                ExonicTags = round(as.numeric(cds_exons_tag_count),2),
                PercIntergenic = round(as.numeric((other_intergenic_tag_count)/total_tags),2)) %>% 
  dplyr::filter(ExonicTags < 4e6 | PercIntergenic > 0.25) %>% dplyr::select(Sample_Name) %>% t() %>% as.vector() #Filter by mapping rate 
#Remove samples from Olah which have bad mapping.
BadMappingOlah <- read.table(file='/local/data/public/zmx21/zmx21_private/GSK/Olah/QC_Reports/star_multiqc_data/multiqc_rseqc_read_distribution.txt',
                         header = T,stringsAsFactors = F) %>% dplyr::select('Sample_Name'='Sample','cds_exons_tag_count','total_tags') %>% 
  dplyr::mutate(PercExonic = round(as.numeric(cds_exons_tag_count/total_tags),2),
                ExonicTags = round(as.numeric(cds_exons_tag_count),2)) %>% 
  dplyr::filter(ExonicTags < 4e6) %>% dplyr::select(Sample_Name) %>% t() %>% as.vector() #Filter by mapping rate 


#Remove bad mapping samples
SalmonTPM_Gene_Merged <- SalmonTPM_Gene_Merged[,!colnames(SalmonTPM_Gene_Merged)%in% c(BadMappingGalatro,BadMappingOlah)]
SalmonTPM_Transcript_Merged <- SalmonTPM_Transcript_Merged[,!colnames(SalmonTPM_Transcript_Merged)%in% c(BadMappingGalatro,BadMappingOlah)]

#Filter according to expression of genes
#The removal criteria is that more than 50% of the samples are unexpressed
#(expression level below the 25th percentile of the respective sample's distribution).
sampleQuantilesGenes <- apply(SalmonTPM_Gene_Merged,2,function(x) quantile(x,0.25))
sampleQuantilesTranscripts <- apply(SalmonTPM_Transcript_Merged,2,function(x) quantile(x,0.25))
genesToKeep <- sapply(1:nrow(SalmonTPM_Gene_Merged),function(i) (sum(SalmonTPM_Gene_Merged[i,] <= sampleQuantilesGenes)/ncol(SalmonTPM_Gene_Merged)) < 0.5)
transcriptsToKeep <- sapply(1:nrow(SalmonTPM_Transcript_Merged),function(i) (sum(SalmonTPM_Transcript_Merged[i,] <= sampleQuantilesTranscripts)/ncol(SalmonTPM_Transcript_Merged)) < 0.5)

SalmonTPM_Gene_Merged <- SalmonTPM_Gene_Merged[rownames(SalmonTPM_Gene_Merged)[genesToKeep],]
SalmonTPM_Transcript_Merged <- SalmonTPM_Transcript_Merged[rownames(SalmonTPM_Transcript_Merged)[transcriptsToKeep],]

#First, correct each batch independently
GalatroSamples_Genes <- SalmonTPM_Gene_Merged[,sapply(colnames(SalmonTPM_Gene_Merged),function(x) substr(x,1,3) == 'GSM')]
GosselinSamples_Genes <- SalmonTPM_Gene_Merged[,sapply(colnames(SalmonTPM_Gene_Merged),function(x) substr(x,1,3) == 'SRR')]
OlahSamples_Genes <- SalmonTPM_Gene_Merged[,sapply(colnames(SalmonTPM_Gene_Merged),function(x) substr(x,1,3) == 'H5K')]
GalatroSamples_Transcripts <- SalmonTPM_Transcript_Merged[,sapply(colnames(SalmonTPM_Transcript_Merged),function(x) substr(x,1,3) == 'GSM')]
GosselinSamples_Transcripts <- SalmonTPM_Transcript_Merged[,sapply(colnames(SalmonTPM_Transcript_Merged),function(x) substr(x,1,3) == 'SRR')]
OlahSamples_Transcripts <- SalmonTPM_Transcript_Merged[,sapply(colnames(SalmonTPM_Transcript_Merged),function(x) substr(x,1,3) == 'H5K')]
RunBatchCorrection <- function(GeneMatrix,TranscriptMatrix,Sample,expType=F,full=F){
  #Filter any rows with all zero counts
  GeneMatrix <- GeneMatrix[apply(GeneMatrix,1,function(x) !all(x==0)),]
  TranscriptMatrix <- TranscriptMatrix[apply(TranscriptMatrix,1,function(x) !all(x==0)),]
  
  MatrixCorrected <- RunCombat(GeneMatrix,TranscriptMatrix,Samples = Sample,expType = expType,full= full)
  PCA_PostCorrection_Gene <- CalcPCA(MatrixCorrected$SalmonTPM_Gene_Combat_Merged)
  PCA_PostCorrection_Gene$Df$batch <- MatrixCorrected$Batches$batch
  
  PCA_PostCorrection_Transcript <- CalcPCA(MatrixCorrected$SalmonTPM_Transcript_Combat_Merged)
  PCA_PostCorrection_Transcript$Df$batch <- MatrixCorrected$Batches$batch
  
  
  PCA_PostCorrectionPlot <- list(autoplot(PCA_PostCorrection_Gene$PCA,data=PCA_PostCorrection_Gene$Df,colour='batch',size=4,shape=F) + ggtitle('Post Batch Correction - Gene'),
                                 autoplot(PCA_PostCorrection_Transcript$PCA,data=PCA_PostCorrection_Transcript$Df,colour='batch',size=4,shape=F) + ggtitle('Post Batch Correction - Transcript'))
                                 
  
  PCA_PreCorrection_Gene <- CalcPCA(GeneMatrix)
  PCA_PreCorrection_Gene$Df$batch <- MatrixCorrected$Batches$batch
  
  PCA_PreCorrection_Transcript <- CalcPCA(TranscriptMatrix)
  PCA_PreCorrection_Transcript$Df$batch <- MatrixCorrected$Batches$batch
  
  
  PCA_PreCorrectionPlot <- list(autoplot(PCA_PreCorrection_Gene$PCA,data=PCA_PreCorrection_Gene$Df,colour='batch',size=4,shape=F) + ggtitle('Pre Batch Correction - Gene'),
                                autoplot(PCA_PreCorrection_Transcript$PCA,data=PCA_PreCorrection_Transcript$Df,colour='batch',size=4,shape=F) + ggtitle('Pre Batch Correction - Transcript'))
  library(egg)
  ggarrange(plots = c(PCA_PreCorrectionPlot,PCA_PostCorrectionPlot),ncol=2)

  return(list(PCA_Plots = list(PCA_PreCorrectionPlot,PCA_PostCorrectionPlot),MatrixCorrected=MatrixCorrected))
}

GalatroCorrected <- RunBatchCorrection(log2(GalatroSamples_Genes+1),log2(GalatroSamples_Transcripts+1),Sample='Galatro',expType = F,full=F)
Gosselin_ExpType_Corrected <- RunBatchCorrection(log2(GosselinSamples_Genes+1),log2(GosselinSamples_Transcripts+1),Sample='Gosselin',expType = T,full=F)
Gosselin_Instrument_Corrected <- RunBatchCorrection(Gosselin_ExpType_Corrected$MatrixCorrected$SalmonTPM_Gene_Combat_Merged,
                                                             Gosselin_ExpType_Corrected$MatrixCorrected$SalmonTPM_Transcript_Combat_Merged,
                                                             Sample='Gosselin',expType = F,full=F)
MergedIndivCorrected <- list(Transcript = 
                               lapply(list(GalatroCorrected,Gosselin_Instrument_Corrected),
                                      function(x) x$MatrixCorrected$SalmonTPM_Transcript_Combat_Merged) %>% {c(.,list(log2(OlahSamples_Transcripts + 1)))},
                             Gene = 
                               lapply(list(GalatroCorrected,Gosselin_Instrument_Corrected),
                                      function(x) x$MatrixCorrected$SalmonTPM_Gene_Combat_Merged) %>% {c(.,list(log2(OlahSamples_Genes+1)))})

#Only take genes/transcripts which are shared across dataset.
SharedGenes <- lapply(MergedIndivCorrected$Gene,rownames) %>% {intersect(.[[1]],intersect(.[[2]],.[[3]]))}
SharedTranscript <- lapply(MergedIndivCorrected$Transcript,rownames) %>% {intersect(.[[1]],intersect(.[[2]],.[[3]]))}
MergedIndivCorrected$Gene <- lapply(MergedIndivCorrected$Gene,function(x) x[SharedGenes,]) %>% {do.call(cbind,.)}
MergedIndivCorrected$Transcript <- lapply(MergedIndivCorrected$Transcript,function(x) x[SharedTranscript,]) %>% {do.call(cbind,.)}
#Remove 0 counts row
MergedIndivCorrected$Gene <- MergedIndivCorrected$Gene[apply(MergedIndivCorrected$Gene,1,function(x) !all(x==0)),]
MergedIndivCorrected$Transcript <- MergedIndivCorrected$Transcript[apply(MergedIndivCorrected$Transcript,1,function(x) !all(x==0)),]

MergedFinalBatchCorrected <- RunBatchCorrection(MergedIndivCorrected$Gene,MergedIndivCorrected$Transcript,c('Galatro','Gosselin','Olah'),expType = F,full = T)
SalmonTPM_Combat_ExpCorrected <- MergedFinalBatchCorrected$MatrixCorrected

#PariseWiseScatterPlots
FinalGalatroCorrected_Gene <- SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged[,sapply(colnames(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged),function(x) substr(x,1,3) == 'GSM')]
FinalGosselinCorrected_Gene <- SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged[,sapply(colnames(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged),function(x) substr(x,1,3) == 'SRR')]
FinalOlahCorrected_Gene <- SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged[,sapply(colnames(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged),function(x) substr(x,1,3) == 'H5K')]

FinalGalatroCorrected_Transcripts <- SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged[,sapply(colnames(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged),function(x) substr(x,1,3) == 'GSM')]
FinalGosselinCorrected_Transcripts <- SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged[,sapply(colnames(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged),function(x) substr(x,1,3) == 'SRR')]
FinalOlahCorrected_Transcripts <- SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged[,sapply(colnames(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged),function(x) substr(x,1,3) == 'H5K')]


PlotPairWiseTPM <- function(params,xsuffix,ysuffix){
  PairDf <- data.frame(med1 = apply(params$dat1,1,function(x) median(x)), med2 = apply(params$dat2,1,function(x) median(x)))
  get_density <- function(x, y, n = 100) {
    dens <- MASS::kde2d(x = x, y = y, n = n)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }
  
  PairDf$density <- get_density(PairDf$med1,PairDf$med2)
  PlotObj <- ggplot(PairDf, aes(x=med1,y=med2,colour=density)) +
    geom_point(size=0.5) +
    geom_smooth(method=lm,   # Add linear regression line
                se=FALSE) + stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
    theme_bw() + labs(x = paste(params$sample[1],xsuffix),y=paste(params$sample[2],ysuffix)) + 
    ggtitle(paste0(params$titlePre,params$sample[2],' VS ',params$sample[1])) +  scale_color_viridis()
  return(PlotObj)
}
library(ggplot2)
library(devtools)
source_gist("524eade46135f6348140",filename = "ggplot_smooth_func.R")
library(parallel)
library(viridis)
CorrectedPairs_Gene <- list(list(dat1=FinalGalatroCorrected_Gene,dat2=FinalGosselinCorrected_Gene,sample=c('Galatro','Gosselin'),titlePre = 'Corrected Gene - \n'),
                       list(dat1=FinalGalatroCorrected_Gene,dat2=FinalOlahCorrected_Gene,sample=c('Galatro','Olah'),titlePre = 'Corrected Gene - \n'),
                       list(dat1=FinalGosselinCorrected_Gene,dat2=FinalOlahCorrected_Gene,sample=c('Gosselin','Olah'),titlePre = 'Corrected Gene - \n'))
BatchCorrectedPlots_Gene <- lapply(CorrectedPairs_Gene,function(x) PlotPairWiseTPM(x,xsuffix = 'Corrected Abundance',ysuffix = 'Corrected Abundance'))
UnCorrectedPairs_Gene <- list(list(dat1=log2(GalatroSamples_Genes+1),dat2=log2(GosselinSamples_Genes+1),sample=c('Galatro','Gosselin'),titlePre = 'Pre Correction Gene -\n'),
                       list(dat1=log2(GalatroSamples_Genes+1),dat2=log2(OlahSamples_Genes+1),sample=c('Galatro','Olah'),titlePre = 'Pre Correction Gene-\n'),
                       list(dat1=log2(GosselinSamples_Genes+1),dat2=log2(OlahSamples_Genes+1),sample=c('Gosselin','Olah'),titlePre = 'Pre Correction Gene-\n'))
UnCorrectedPlots_Gene <- lapply(UnCorrectedPairs_Gene,function(x) PlotPairWiseTPM(x,xsuffix='log(TPM+1)',ysuffix='log(TPM+1)'))
egg::ggarrange(plots=c(UnCorrectedPlots_Gene,BatchCorrectedPlots_Gene),ncol=3)


CorrectedPairs_Transcript<- list(list(dat1=FinalGalatroCorrected_Transcript,dat2=FinalGosselinCorrected_Transcript,sample=c('Galatro','Gosselin'),titlePre = 'Corrected Transcript - \n'),
                            list(dat1=FinalGalatroCorrected_Transcript,dat2=FinalOlahCorrected_Transcript,sample=c('Galatro','Olah'),titlePre = 'Corrected Transcript - \n'),
                            list(dat1=FinalGosselinCorrected_Transcript,dat2=FinalOlahCorrected_Transcript,sample=c('Gosselin','Olah'),titlePre = 'Corrected Transcript - \n'))
BatchCorrectedPlots_Transcript <- mclapply(CorrectedPairs_Transcript,function(x) PlotPairWiseTPM(x,xsuffix = 'Corrected Abundance',ysuffix = 'Corrected Abundance'),mc.cores = 3)
UnCorrectedPairs_Transcript <- list(list(dat1=log2(GalatroSamples_Transcripts+1),dat2=log2(GosselinSamples_Transcripts+1),sample=c('Galatro','Gosselin'),titlePre = 'Pre Correction Transcript -\n'),
                              list(dat1=log2(GalatroSamples_Transcripts+1),dat2=log2(OlahSamples_Transcripts+1),sample=c('Galatro','Olah'),titlePre = 'Pre Correction Transcript-\n'),
                              list(dat1=log2(GosselinSamples_Transcripts+1),dat2=log2(OlahSamples_Transcripts+1),sample=c('Gosselin','Olah'),titlePre = 'Pre Correction Transcript-\n'))
UnCorrectedPlots_Transcript <- lapply(UnCorrectedPairs_Transcript,function(x) PlotPairWiseTPM(x,xsuffix='log(TPM+1)',ysuffix='log(TPM+1)'))
egg::ggarrange(plots=c(UnCorrectedPlots_Transcript,BatchCorrectedPlots_Transcript),ncol=3)


WithinSamplePairs <-  list(list(dat1=log2(GalatroSamples_Genes[rownames(FinalGalatroCorrected_Gene),]+1),dat2=FinalGalatroCorrected_Gene,sample=c('Galatro Raw','Galatro Corrected'),titlePre = 'Gene -\n'),
                           list(dat1=log2(GosselinSamples_Genes[rownames(FinalGosselinCorrected_Gene),]+1),dat2=FinalGosselinCorrected_Gene,sample=c('Gosselin Raw','Gosselin Corrected'),titlePre = 'Gene-\n'),
                           list(dat1=log2(OlahSamples_Genes[rownames(FinalOlahCorrected_Gene),]+1),dat2=FinalOlahCorrected_Gene,sample=c('Olah Raw','Olah Post Correction'),titlePre = 'Gene-\n'),
                           list(dat1=log2(GalatroSamples_Transcripts[rownames(FinalGalatroCorrected_Transcript),]+1),dat2=FinalGalatroCorrected_Transcript,sample=c('Galatro Raw','Galatro Corrected'),titlePre = 'Transcript-\n'),
                           list(dat1=log2(GosselinSamples_Transcripts[rownames(FinalGosselinCorrected_Transcript),]+1),dat2=FinalGosselinCorrected_Transcript,sample=c('Gosselin Raw','Gosselin Corrected'),titlePre = 'Transcript-\n'),
                           list(dat1=log2(OlahSamples_Transcripts[rownames(FinalOlahCorrected_Transcript),]+1),dat2=FinalOlahCorrected_Transcript,sample=c('Olah Raw','Olah Corrected'),titlePre = 'Transcript-\n'))

WithinSamplePlots <- lapply(WithinSamplePairs,function(x) PlotPairWiseTPM(x,xsuffix='log(TPM+1)',ysuffix='Abundance'))
egg::ggarrange(plots=WithinSamplePlots,ncol=3)

meanAbundancesGene <- apply(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged,1,function(x) mean(x))
cvGene <- apply(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged,1,function(x) sd(x)/mean(x))
meanAbundancesTranscript <- apply(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged,1,function(x) mean(x))
cvTranscript <- apply(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged,1,function(x) sd(x)/mean(x))

geneLevelCutOffTPM <- quantile(meanAbundancesGene,prob=seq(0.02,0.2,0.02)) %>% round(2)
geneLevelCutOffCV <- 0.3
geneLevelCutOff <- expand.grid(TPM=geneLevelCutOffTPM,CV=geneLevelCutOffCV)

geneLevelFiltCounts <- lapply(1:nrow(geneLevelCutOff), function(i) FilterCountMatrix(countMatrix = SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged,
                                                                                          meanAbundances = meanAbundancesGene,
                                                                                          cv = cvGene,
                                                                                          cvCutOffAbsolute = geneLevelCutOff[i,1],
                                                                                          meanCutOffAbsolute = geneLevelCutOff[i,2],
                                                                                          mappingTable = gtfTables$geneTable) %>%
                                                                                          {ParseBiotypeTable(.$categoricalCounts)} %>%
                                                                                          {table(.$gene_biotype)}) %>% {do.call(rbind,.)} %>% {cbind(geneLevelCutOff,.)}
transcriptLevelFiltCounts <- lapply(1:nrow(geneLevelFiltCounts), function(i) FilterCountMatrix(countMatrix = SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged,
                                                                                          meanAbundances = meanAbundancesTranscript,
                                                                                          cv = cvTranscript,
                                                                                          cvCutOffAbsolute = geneLevelCutOff[i,1],
                                                                                          meanCutOffAbsolute = geneLevelCutOff[i,2],
                                                                                          mappingTable = gtfTables$transcriptTable) %>%
                                                                                          {ParseBiotypeTable(.$categoricalCounts)} %>%
                                                                                          {table(.$source)}) %>% {do.call(rbind,.)} %>% {cbind(geneLevelCutOff,.)}
geneLevelFiltCounts$TPM_Percentile <- seq(2,20,2)
geneLevelFiltCounts <- geneLevelFiltCounts[c(7,1:6)]

transcriptLevelFiltCounts$TPM_Percentile <- seq(2,20,2)
transcriptLevelFiltCounts <- transcriptLevelFiltCounts[c(7,1:6)]

library(gridExtra)
grid.table(round(geneLevelFiltCounts,3))
grid.table(round(transcriptLevelFiltCounts,3))


#Bias in read Distribution
MappingRatePlot <- ggplot() + aes(x=CollectMetadata(SalmonTPM_Gene_Merged,full=F)$mappingRate,y=apply(SalmonTPM_Gene_Merged,2,median),
                                  color=CollectMetadata(SalmonTPM_Gene_Merged,full=F)$dataset) + geom_point() + 
                              labs(color='Dataset',x = 'Mapping Rate',y='Median TPM') + ggtitle('Median TPM VS Mapping Rate ')

TotalReadsPlot <- ggplot() + aes(x=CollectMetadata(SalmonTPM_Gene_Merged,full=F)$numReads,y=apply(SalmonTPM_Gene_Merged,2,median),
                                  color=CollectMetadata(SalmonTPM_Gene_Merged,full=F)$dataset) + geom_point() + 
  labs(color='Dataset',x = 'Num Reads',y='Median TPM') + ggtitle('Median TPM VS Num Reads ')

IntronicRatePlot <- ggplot() + aes(x=CollectMetadata(SalmonTPM_Gene_Merged,full=F)$intronTags/CollectMetadata(SalmonTPM_Gene_Merged,full=F)$totalTags,y=apply(SalmonTPM_Gene_Merged,2,median),
                                  color=CollectMetadata(SalmonTPM_Gene_Merged,full=F)$dataset) + geom_point() + 
  labs(color='Dataset',x = 'Intronic Tags',y='Median TPM') + ggtitle('Median TPM VS Intronic Reads ')
ggarrange(plots=list(MappingRatePlot,TotalReadsPlot,IntronicRatePlot),ncol=3)


# #Remove genes which have 0 expression across all Samples
# SalmonTPM_Gene_Merged <- SalmonTPM_Gene_Merged[!apply(SalmonTPM_Gene_Merged,1,function(x) all(x==0)),]
# SalmonTPM_Transcript_Merged <- SalmonTPM_Transcript_Merged[!apply(SalmonTPM_Transcript_Merged,1,function(x) all(x==0)),]
# 
# SalmonTPM_Combat <- RunCombat(SalmonTPM_Gene_Merged,SalmonTPM_Transcript_Merged,Samples = c('Galatro','Gosselin'),expType = F,full=F)
# SalmonTPM_Combat_ExpCorrected <- RunCombat(SalmonTPM_Gene_Merged,SalmonTPM_Transcript_Merged,Samples = c('Galatro','Gosselin'),expType = T,full=F)
# 
# save(SalmonTPM_Combat,file='../Count_Data/Batch_Corrected/Salmon_TPM_Combat_Par.rda')
# save(SalmonTPM_Combat_ExpCorrected,file='../Count_Data/Batch_Corrected/SalmonTPM_Combat_Par_ExpCorrected.rda')
# 
# #Look at Biotypes under different cutoffs
# meanAbundancesGene <- apply(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged,1,function(x) mean(x))
# cvGene <- apply(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged,1,function(x) sd(x) / mean(x))
# #Get Transcript Level Stats
# meanAbundancesTranscript <- apply(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged,1,function(x) mean(x))
# cvTranscript <- apply(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged,1,function(x) sd(x) / mean(x))
# 
# geneLevelCutOffTPM <- quantile(meanAbundancesGene,prob=seq(0.1,0.3,0.1)) %>% round(2)
# geneLevelCutOffCV <- quantile(cvGene,prob=seq(0.1,0.3,0.1)) %>% round(2)
# geneLevelCutOff <- expand.grid(TPM=geneLevelCutOffTPM,CV=geneLevelCutOffCV)
# 
# geneLevelFiltCounts <- lapply(1:nrow(geneLevelCutOff), function(i) FilterCountMatrix(countMatrix = SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged,
#                                                                                           meanAbundances = meanAbundancesGene,
#                                                                                           cv = cvGene,
#                                                                                           cvCutOffAbsolute = geneLevelCutOff[i,1],
#                                                                                           meanCutOffAbsolute = geneLevelCutOff[i,2],
#                                                                                           mappingTable = gtfTables$geneTable) %>%
#                                                                                           {ParseBiotypeTable(.$categoricalCounts)} %>%
#                                                                                           {table(.$gene_biotype)}) %>% {do.call(rbind,.)} %>% {cbind(geneLevelCutOff,.)}
# transcriptLevelFiltCounts <- lapply(1:nrow(geneLevelFiltCounts), function(i) FilterCountMatrix(countMatrix = SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged,
#                                                                                           meanAbundances = meanAbundancesTranscript,
#                                                                                           cv = cvTranscript,
#                                                                                           cvCutOffAbsolute = geneLevelCutOff[i,1],
#                                                                                           meanCutOffAbsolute = geneLevelCutOff[i,2],
#                                                                                           mappingTable = gtfTables$transcriptTable) %>%
#                                                                                           {ParseBiotypeTable(.$categoricalCounts)} %>%
#                                                                                           {table(.$source)}) %>% {do.call(rbind,.)} %>% {cbind(geneLevelCutOff,.)}
# library(gridExtra)
# grid.table(round(geneLevelFiltCounts,3))
# grid.table(round(transcriptLevelFiltCounts,3))
# 
# 
# #Plot resulting boxplot before batch correction according to mapping rate. 
# TPMVsReads <- ggplot(stack(as.data.frame(SalmonTPM_Gene_Merged)) %>% 
#                dplyr::left_join(CollectMetadata(SalmonTPM_Gene_Merged,full=F),by=c('ind' = 'Sample_Name'))) + 
#   geom_boxplot(aes(x = ind, y = values,fill=numReads),outlier.shape = NA) + scale_x_discrete(name = "Sample") +
#   scale_y_continuous(name = "TPM",limits = c(0,4)) + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.6,size =15)) + scale_colour_gradientn(colours = rainbow(7))
# 
# source('pca_analysis.R')
# library(egg)
# #Plot PCA including expType and not including expType
# PCANoExpType <- CalcPCA(SalmonTPM_Combat$SalmonTPM_Gene_Combat_Merged)
# PCAPlotNoExpType_Instrument <- autoplot(PCANoExpType$PCA, data = PCANoExpType$Df, colour = 'Instrument',size=4,shape=F) + 
#   ggtitle('PCA of Non-Exp Type Corrected - Instrument')
# PCAPlotNoExpType_expType <- autoplot(PCANoExpType$PCA, data = PCANoExpType$Df, colour = 'expType',size=4,shape=F) + 
#                                    ggtitle('PCA of Non-Exp Type Corrected - Exp Type')
# PCAPlotNoExpType <- ggarrange(PCAPlotNoExpType_Instrument,PCAPlotNoExpType_expType,ncol=2)
# 
# PCAExpType <- CalcPCA(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged)
# PCAPlotExpType_Instrument <- autoplot(PCAExpType$PCA, data = PCAExpType$Df, colour = 'Instrument',size=4,shape=F) + ggtitle('PCA of Exp Type Corrected - Instrument')
# PCAPlotExpType_ExpType <- autoplot(PCAExpType$PCA, data = PCAExpType$Df, colour = 'expType',size=4,shape=F) + ggtitle('PCA of Exp Type Corrected - ExpType')
# PCAPlotExpType <- ggarrange(PCAPlotExpType_Instrument,PCAPlotExpType_ExpType,ncol=2)
# 
# #Boxplot of batch corrected  TPM
# TPMExpCorrected<- ggplot(stack(as.data.frame(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged)) %>% 
#                        dplyr::left_join(CollectMetadata(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged,full=F),by=c('ind' = 'Sample_Name'))) + 
#   geom_boxplot(aes(x = ind, y = values,fill=expType),outlier.shape = NA) + scale_x_discrete(name = "Sample") +
#   scale_y_continuous(name = "log(TPM+1)",limits = c(-2,4)) + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.6,size =15)) + scale_colour_gradientn(colours = rainbow(7)) + ggtitle('With Exp Type Correction')
# TPMNonExpCorrected<- ggplot(stack(as.data.frame(SalmonTPM_Combat$SalmonTPM_Gene_Combat_Merged)) %>% 
#                            dplyr::left_join(CollectMetadata(SalmonTPM_Combat$SalmonTPM_Gene_Combat_Merged,full=F),by=c('ind' = 'Sample_Name'))) + 
#   geom_boxplot(aes(x = ind, y = values,fill=expType),outlier.shape = NA) + scale_x_discrete(name = "Sample") +
#   scale_y_continuous(name = "log(TPM+1)",limits = c(-2,4)) + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.6,size =15)) + scale_colour_gradientn(colours = rainbow(7)) + ggtitle('No Exp Type Correction')
# TPMExpCorVsNonExpCor <-  ggarrange(TPMExpCorrected,TPMNonExpCorrected,nrow=2)
# 
# #Scatterplot of median exp of genes Galatro vs. Gosselin before batch correction 
# GalatroSamples <- log(SalmonTPM_Gene_Merged[,sapply(colnames(SalmonTPM_Gene_Merged),function(x) substr(x,1,3) == 'GSM')]+1)
# GosselinSamples <- log(SalmonTPM_Gene_Merged[,sapply(colnames(SalmonTPM_Gene_Merged),function(x) substr(x,1,3) == 'SRR')]+1)
# GalatroVsGosselin <- data.frame(Gosselin = apply(GosselinSamples,1,function(x) median(x)), Galatro = apply(GalatroSamples,1,function(x) median(x)))
# 
# library(devtools)
# source_gist("524eade46135f6348140",filename = "ggplot_smooth_func.R")
# GalatroVsGosselinPlot <- ggplot(GalatroVsGosselin, aes(x=Gosselin,y=Galatro)) +
#   geom_point() +   
#   geom_smooth(method=lm,   # Add linear regression line
#               se=FALSE) + stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) + 
#   theme_bw() + labs(x = 'Gosselin log(TPM+1)',y='Galatro log(TPM+1)') + ggtitle('Median Exp - Galatro vs Gosselin') 
# 
# # GalatroVsRandom <- data.frame(Galatro = GalatroVsGosselin$Galatro,Random = GalatroVsGosselin$Gosselin[sample(1:nrow(GalatroVsGosselin),size=nrow(GalatroVsGosselin),replace = F)])
# # GalatroVsRandomPlot <-ggplot(GalatroVsRandom, aes(x=Random,y=Galatro)) +
# #   geom_point() +   
# #   geom_smooth(method=lm,   # Add linear regression line
# #               se=FALSE,formula = formula) + stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) + 
# #   theme_bw() + labs(x = 'Gosselin Permuted log(TPM+1)',y='Galatro log(TPM+1)') + ggtitle('Median Exp - Galatro vs Gosselin') 
# # ggarrange(GalatroVsGosselinPlot,GalatroVsRandomPlot,ncol=2)
# 
# #Scatterplot of median exp of genes Galatro vs. Gosselin after batch correction. 
# GalatroSamplesCorrected <- SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged[,sapply(colnames(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged),function(x) substr(x,1,3) == 'GSM')]
# GosselinSamplesCorrected <- SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged[,sapply(colnames(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged),function(x) substr(x,1,3) == 'SRR')]
# GalatroVsGosselinCorrected <- data.frame(Gosselin = apply(GosselinSamplesCorrected,1,function(x) median(x)), Galatro = apply(GalatroSamplesCorrected,1,function(x) median(x)))
# GalatroVsGosselinPlotCorrected <- ggplot(GalatroVsGosselinCorrected, aes(x=Gosselin,y=Galatro)) +
#   geom_point() +   
#   geom_smooth(method=lm,   # Add linear regression line
#               se=FALSE) + stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) + 
#   theme_bw() + labs(x = 'Gosselin log(TPM+1)',y='Galatro log(TPM+1)') + ggtitle('Batch Corrected Median Exp - Galatro vs Gosselin') 
# 
# # GalatroVsRandomCorrected <- data.frame(Galatro = GalatroVsGosselinCorrected$Galatro,Random = GalatroVsGosselinCorrected$Gosselin[sample(1:nrow(GalatroVsGosselinCorrected),size=nrow(GalatroVsGosselinCorrected),replace = F)])
# # GalatroVsRandomPlotCorrected <- ggplot(GalatroVsRandomCorrected, aes(x=Random,y=Galatro)) +
# #   geom_point() +   
# #   geom_smooth(method=lm,   # Add linear regression line
# #               se=FALSE,formula = formula) + stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) + 
# #   theme_bw() + labs(x = 'Gosselin Permuted log(TPM+1)',y='Galatro log(TPM+1)') + ggtitle('Batch Corrected Median Exp - Galatro vs Gosselin') 
# ggarrange(GalatroVsGosselinPlot,GalatroVsGosselinPlotCorrected,ncol=2)
# 
# ####################################################COMPARE WITH WHOLE BRAIN###############################################################
# 
# #Boxplot of Brain Samples
# load('../Count_Data/Galatro/SalmonTPM_Gene_WholeBrain.rda')
# BulkBrainMetadata <- CollectMetadata(SalmonTPM_Gene_WholeBrain)
# BulkBrainTPM <- ggplot(stack(as.data.frame(SalmonTPM_Gene_WholeBrain)) %>% 
#                        dplyr::left_join(CollectMetadata(SalmonTPM_Gene_WholeBrain),by=c('ind' = 'Sample_Name')) %>% dplyr::mutate(percIntergenic = intergenicTags/totalTags)) + 
#   geom_boxplot(aes(x = ind, y = values,fill=exonTags > 4e6 & percIntergenic < 0.25),outlier.shape = NA) + scale_x_discrete(name = "Sample") +
#   scale_y_continuous(name = "TPM",limits = c(0,4)) + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.6,size =15)) + scale_colour_gradientn(colours = rainbow(7))
# BulkBrainPCA <- CalcPCA(SalmonTPM_Gene_WholeBrain)
# BulkBrainPCA$Df$percIntergenic <- BulkBrainPCA$Df$intergenicTags / BulkBrainPCA$Df$totalTags
# BulkBrainPCAPlot <- autoplot(BulkBrainPCA$PCA, data = BulkBrainPCA$Df,colour='percIntergenic',size=4,shape=F) + ggtitle('PCA of Bulk Brain') 
# 
# #Remove bad Samples, and compare with Galatro and Gosselin
# SalmonTPM_Gene_WholeBrain <- SalmonTPM_Gene_WholeBrain[,BulkBrainMetadata %>%
#                                                          dplyr::mutate(percIntergenic = intergenicTags/totalTags) %>%
#                                                          dplyr::filter(percIntergenic < 0.25) %>% dplyr::select(Sample_Name) %>% t() %>% as.vector()]
# SalmonTPM_Gene_WholeBrain <- SalmonTPM_Gene_WholeBrain[rownames(GalatroSamples),]
# WholeBrainVsGalatro <- data.frame(Galatro = apply(GalatroSamples,1,function(x) median(x)),WholeBrain = apply(log(SalmonTPM_Gene_WholeBrain+1),1,function(x) median(x)))
# WholeBrainVsGosselin <- data.frame(Gosselin = apply(GosselinSamples,1,function(x) median(x)),WholeBrain = apply(log(SalmonTPM_Gene_WholeBrain+1),1,function(x) median(x)))
# WholeBrainVsGalatroPlot <- ggplot(WholeBrainVsGalatro, aes(x=WholeBrain,y=Galatro)) +
#   geom_point() +   
#   geom_smooth(method=lm,   # Add linear regression line
#               se=FALSE) + stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) + 
#   theme_bw() + labs(x = 'Whole Brain log(TPM+1)',y='Galatro log(TPM+1)') + ggtitle('Median Exp - Galatro vs WholeBrain') 
# 
# WholeBrainVsGosselinPlot <- ggplot(WholeBrainVsGosselin, aes(x=WholeBrain,y=Gosselin)) +
#   geom_point() +   
#   geom_smooth(method=lm,   # Add linear regression line
#               se=FALSE) + stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) + 
#   theme_bw() + labs(x = 'Whole Brain log(TPM+1)',y='Gosselin log(TPM+1)') + ggtitle('Median Exp - Gosselin vs WholeBrain') 
# 
# ggarrange(WholeBrainVsGalatroPlot,WholeBrainVsGosselinPlot,ncol=2)
# 
# #PCA plot comparing within Galatro
# PCAWholeBrainVsGalatro <- CalcPCA(cbind(log(SalmonTPM_Gene_WholeBrain+1),GalatroSamples))
# PCAWholeBrainVsGalatro$Df$Dataset <- ifelse(PCAWholeBrainVsGalatro$Df$readLength == 202 & PCAWholeBrainVsGalatro$Df$BulkBrain == 0,'GalatroIllumina',ifelse(PCAWholeBrainVsGalatro$Df$readLength == 252 & PCAWholeBrainVsGalatro$Df$BulkBrain == 0,'GalatroTakara','Whole Brain'))
# autoplot(PCAWholeBrainVsGalatro$PCA, data = PCAWholeBrainVsGalatro$Df, colour = 'Dataset',size=4,shape=F) + 
#   ggtitle('PCA of Tissue Type within Galatro')
# 
# #PCA plot adjusting according to Batch1 and Batch2
# WholeBrainandGalatro <- cbind(log(SalmonTPM_Gene_WholeBrain+1),GalatroSamples)
# WholeBrainandGalatro <- WholeBrainandGalatro[!apply(WholeBrainandGalatro,1,function(x) all(x==0)),]
# BatchCorrectedWholeBrainGalatro <- RunCombat(WholeBrainandGalatro,Samples='Galatro')
# PCABatchCorrected <- CalcPCA(BatchCorrectedWholeBrainGalatro$SalmonTPM_Gene_Combat_Merged)
# PCABatchCorrected$Df$Dataset <- ifelse(PCAWholeBrainVsGalatro$Df$readLength == 202 & PCAWholeBrainVsGalatro$Df$BulkBrain == 0,'GalatroIllumina',ifelse(PCAWholeBrainVsGalatro$Df$readLength == 252 & PCAWholeBrainVsGalatro$Df$BulkBrain == 0,'GalatroTakara','Whole Brain'))
# autoplot(PCABatchCorrected$PCA, data = PCABatchCorrected$Df, colour = 'Dataset',size=4,shape=F) + 
#   ggtitle('PCA of Tissue Type within Galatro - Batch Corrected')

