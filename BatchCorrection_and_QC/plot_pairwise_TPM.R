#Helper functin to extract a dataset from a merged matrix.
ExtractDataset <- function(countMatrix,Sample){
  if(Sample=='Galatro'){
    countMatrix[,sapply(colnames(countMatrix),function(x) substr(x,1,3) == 'GSM')]
  }else if (Sample=='Gosselin'){
    countMatrix[,sapply(colnames(countMatrix),function(x) substr(x,1,3) == 'SRR')]
  }else if (Sample=='Olah'){
    countMatrix[,sapply(colnames(countMatrix),function(x) substr(x,1,3) == 'H5K')]
  }
}

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
SavePairwiseTPMPlots <- function(){
  load('../../CodeImages/TPMFiltering_0.75.RData')
  load('../../Count_Data/Batch_Corrected/SalmonTPM_Combat_ExpCorrected.rda')
  library(ggplot2)
  library(devtools)
  source_gist("524eade46135f6348140",filename = "ggplot_smooth_func.R")
  library(parallel)
  library(viridis)
  
  #Get All Corrected Genes
  FinalGalatroCorrected_Gene <- ExtractDataset(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged,'Galatro')
  FinalGosselinCorrected_Gene <- ExtractDataset(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged,'Gosselin')
  FinalOlahCorrected_Gene <- ExtractDataset(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged,'Olah')
  
  FinalGalatroCorrected_Transcripts <- ExtractDataset(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged,'Galatro')
  FinalGosselinCorrected_Transcripts <- ExtractDataset(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged,'Gosselin')
  FinalOlahCorrected_Transcripts <- ExtractDataset(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged,'Olah')
  
  #Within Sample comparison, batch corrected vs not.
  WithinSamplePairs <-  list(list(dat1=log2(GalatroSamples_Genes[rownames(FinalGalatroCorrected_Gene),]+1),dat2=FinalGalatroCorrected_Gene,sample=c('Galatro Raw','Galatro Corrected'),titlePre = 'Gene -\n'),
                             list(dat1=log2(GosselinSamples_Genes[rownames(FinalGosselinCorrected_Gene),]+1),dat2=FinalGosselinCorrected_Gene,sample=c('Gosselin Raw','Gosselin Corrected'),titlePre = 'Gene-\n'),
                             list(dat1=log2(OlahSamples_Genes[rownames(FinalOlahCorrected_Gene),]+1),dat2=FinalOlahCorrected_Gene,sample=c('Olah Raw','Olah Corrected'),titlePre = 'Gene-\n'),
                             list(dat1=log2(GalatroSamples_Transcripts[rownames(FinalGalatroCorrected_Transcripts),]+1),dat2=FinalGalatroCorrected_Transcripts,sample=c('Galatro Raw','Galatro Corrected'),titlePre = 'Transcript-\n'),
                             list(dat1=log2(GosselinSamples_Transcripts[rownames(FinalGosselinCorrected_Transcripts),]+1),dat2=FinalGosselinCorrected_Transcripts,sample=c('Gosselin Raw','Gosselin Corrected'),titlePre = 'Transcript-\n'),
                             list(dat1=log2(OlahSamples_Transcripts[rownames(FinalOlahCorrected_Transcripts),]+1),dat2=FinalOlahCorrected_Transcripts,sample=c('Olah Raw','Olah Corrected'),titlePre = 'Transcript-\n'))
  
  WithinSamplePlots <- lapply(WithinSamplePairs,function(x) PlotPairWiseTPM(x,xsuffix='log(TPM+1)',ysuffix='Abundance'))
  tiff(filename = '../../Figures/WithinSample_TPM.tiff',width = 1000,height = 550)
  egg::ggarrange(plots=WithinSamplePlots,ncol=3)
  dev.off()
  
  #Inter Sample comparison, Gene Level. 
  CorrectedPairs_Gene <- list(list(dat1=FinalGalatroCorrected_Gene,dat2=FinalGosselinCorrected_Gene,sample=c('Galatro','Gosselin'),titlePre = 'Corrected Gene - \n'),
                              list(dat1=FinalGalatroCorrected_Gene,dat2=FinalOlahCorrected_Gene,sample=c('Galatro','Olah'),titlePre = 'Corrected Gene - \n'),
                              list(dat1=FinalGosselinCorrected_Gene,dat2=FinalOlahCorrected_Gene,sample=c('Gosselin','Olah'),titlePre = 'Corrected Gene - \n'))
  BatchCorrectedPlots_Gene <- lapply(CorrectedPairs_Gene,function(x) PlotPairWiseTPM(x,xsuffix = 'Corrected Abundance',ysuffix = 'Corrected Abundance'))
  UnCorrectedPairs_Gene <- list(list(dat1=log2(GalatroSamples_Genes+1),dat2=log2(GosselinSamples_Genes+1),sample=c('Galatro','Gosselin'),titlePre = 'Pre Correction Gene -\n'),
                                list(dat1=log2(GalatroSamples_Genes+1),dat2=log2(OlahSamples_Genes+1),sample=c('Galatro','Olah'),titlePre = 'Pre Correction Gene-\n'),
                                list(dat1=log2(GosselinSamples_Genes+1),dat2=log2(OlahSamples_Genes+1),sample=c('Gosselin','Olah'),titlePre = 'Pre Correction Gene-\n'))
  UnCorrectedPlots_Gene <- lapply(UnCorrectedPairs_Gene,function(x) PlotPairWiseTPM(x,xsuffix='log(TPM+1)',ysuffix='log(TPM+1)'))
  tiff(filename = '../../Figures/InterSample_Gene_TPM.tiff',width = 1000,height = 550)
  egg::ggarrange(plots=c(UnCorrectedPlots_Gene,BatchCorrectedPlots_Gene),ncol=3)
  dev.off()
  
  #Inter Sample Comparisons, Transcript Level
  CorrectedPairs_Transcripts<- list(list(dat1=FinalGalatroCorrected_Transcripts,dat2=FinalGosselinCorrected_Transcripts,sample=c('Galatro','Gosselin'),titlePre = 'Corrected Transcript - \n'),
                                    list(dat1=FinalGalatroCorrected_Transcripts,dat2=FinalOlahCorrected_Transcripts,sample=c('Galatro','Olah'),titlePre = 'Corrected Transcript - \n'),
                                    list(dat1=FinalGosselinCorrected_Transcripts,dat2=FinalOlahCorrected_Transcripts,sample=c('Gosselin','Olah'),titlePre = 'Corrected Transcript - \n'))
  BatchCorrectedPlots_Transcripts <- lapply(CorrectedPairs_Transcripts,function(x) PlotPairWiseTPM(x,xsuffix = 'Corrected Abundance',ysuffix = 'Corrected Abundance'))
  UnCorrectedPairs_Transcripts <- list(list(dat1=log2(GalatroSamples_Transcripts+1),dat2=log2(GosselinSamples_Transcripts+1),sample=c('Galatro','Gosselin'),titlePre = 'Pre Correction Transcript -\n'),
                                       list(dat1=log2(GalatroSamples_Transcripts+1),dat2=log2(OlahSamples_Transcripts+1),sample=c('Galatro','Olah'),titlePre = 'Pre Correction Transcript-\n'),
                                       list(dat1=log2(GosselinSamples_Transcripts+1),dat2=log2(OlahSamples_Transcripts+1),sample=c('Gosselin','Olah'),titlePre = 'Pre Correction Transcript-\n'))
  UnCorrectedPlots_Transcripts <- lapply(UnCorrectedPairs_Transcripts,function(x) PlotPairWiseTPM(x,xsuffix='log(TPM+1)',ysuffix='log(TPM+1)'))
  tiff(filename = '../../Figures/InterSample_Transcript_TPM.tiff',width = 1000,height = 550)
  egg::ggarrange(plots=c(UnCorrectedPlots_Transcripts,BatchCorrectedPlots_Transcripts),ncol=3)
  dev.off()
  
  #Compare Raw Microglia Expression with Raw Brain expression.
  BrainComparison <- list(list(dat1=log2(TPM_WholeBrain_Gene[intersect(rownames(TPM_WholeBrain_Gene),rownames(GalatroSamples_Genes)),]+1),
                               dat2=GalatroSamples_Genes[intersect(rownames(TPM_WholeBrain_Gene),rownames(GalatroSamples_Genes)),],sample=c('Brain','Galatro Raw'),titlePre = 'Gene -\n'),
                          list(dat1=log2(TPM_WholeBrain_Gene[intersect(rownames(TPM_WholeBrain_Gene),rownames(GosselinSamples_Genes)),]+1),
                               dat2=GosselinSamples_Genes[intersect(rownames(TPM_WholeBrain_Gene),rownames(GosselinSamples_Genes)),],sample=c('Brain','Gosselin Corrected'),titlePre = 'Gene-\n'),
                          list(dat1=log2(TPM_WholeBrain_Gene[intersect(rownames(TPM_WholeBrain_Gene),rownames(OlahSamples_Genes)),]+1),
                               dat2=OlahSamples_Genes[intersect(rownames(TPM_WholeBrain_Gene),rownames(OlahSamples_Genes)),],sample=c('Brain','Olah Corrected'),titlePre = 'Gene-\n'),
                          list(dat1=log2(TPM_WholeBrain_Transcript[intersect(rownames(TPM_WholeBrain_Transcript),rownames(GalatroSamples_Transcripts)),]+1),
                               dat2=GalatroSamples_Transcripts[intersect(rownames(TPM_WholeBrain_Transcript),rownames(GalatroSamples_Transcripts)),],sample=c('Brain','Galatro Corrected'),titlePre = 'Transcript-\n'),
                          list(dat1=log2(TPM_WholeBrain_Transcript[intersect(rownames(TPM_WholeBrain_Transcript),rownames(GosselinSamples_Transcripts)),]+1),
                               dat2=GosselinSamples_Transcripts[intersect(rownames(TPM_WholeBrain_Transcript),rownames(GosselinSamples_Transcripts)),],sample=c('Brain','Gosselin Corrected'),titlePre = 'Transcript-\n'),
                          list(dat1=log2(TPM_WholeBrain_Transcript[intersect(rownames(TPM_WholeBrain_Transcript),rownames(OlahSamples_Transcripts)),]+1),
                               dat2=OlahSamples_Transcripts[intersect(rownames(TPM_WholeBrain_Transcript),rownames(OlahSamples_Transcripts)),],sample=c('Brain','Olah Corrected'),titlePre = 'Transcript-\n'))
  BrainComparisonPlots <- lapply(BrainComparison,function(x) PlotPairWiseTPM(x,xsuffix='log(TPM+1)',ysuffix='log(TPM+1)'))
  tiff(filename = '../../Figures/MicrogliaVSBrain_TPM.tiff',width = 1000,height = 550)
  egg::ggarrange(plots=BrainComparisonPlots,ncol=3)
  dev.off()
  # save(list = ls(environment()),file='../../CodeImages/PairWisePlots.RData')
  
}
