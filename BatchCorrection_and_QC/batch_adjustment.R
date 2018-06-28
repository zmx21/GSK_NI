#Batch correction, 
#input should be count matrix, expType specifies whether batch correction should be performed on 
#expType of Gosselin. Full specifies whether batch correction should be performed on dataset alone. 
#output is Combat and PC corrected count matrix
RunCombat <- function(Gene_Matrix,Transcript_Matrix = NULL,Samples,expType = F,full=F){
  library(dplyr)
  #Batch Adjustment using Combat. Input is count Matrix (where rows are genes and columns are samples)
  #Metadata is data.frame, with mapping between batch and sample. 
  BatchAdjustCombat <- function(countMatrix,metadata){
    library(sva)
    #Add no other covariates to Combat batch correct
    modcombat <- model.matrix(~1,data=metadata)
    #Do batch correction using EB of Combat, based on batch as provided in metadata
    combat_edata = sva::ComBat(dat=countMatrix, batch=metadata$batch, mod=modcombat,prior.plots = F,par.prior = T)
    return(combat_edata)
  }
  
  #Load batch Information of each sample, Galatro et al
  runTable_Galatro <- rbind(read.table(file = '/local/data/public/zmx21/zmx21_private/GSK/Galatro/SraRunTable.txt',header = T,sep = '\t',stringsAsFactors = F),
                            read.table(file = '/local/data/public/zmx21/zmx21_private/GSK/Galatro_Brain/SraRunTable.txt',header = T,sep = '\t',stringsAsFactors = F))
  #No batches in Olah
  OlahBatches <- data.frame(Sample_Name = colnames(Gene_Matrix)[sapply(colnames(Gene_Matrix),function(x) grepl('H5KN',x))])
  OlahBatches$batch <- rep('Olah',nrow(OlahBatches))
  
  
  #Load batch Information of each sample, Galatro et al
  runTable_Gosselin <- read.table(file = '/local/data/public/zmx21/zmx21_private/GSK/Gosselin/SraRunTable_Parsed.txt',header = T,sep = '\t')
  runTable_Gosselin <- dplyr::mutate(runTable_Gosselin,expType = ifelse(grepl('ExVivo',Library_Name),'ExVivo','InVitro'))
  
  #Choose which samples to include, based on provided Samples character vector.
  GalatroSamples <- which(colnames(Gene_Matrix)%in%runTable_Galatro$Sample_Name)
  GosselinSamples <- which(colnames(Gene_Matrix)%in%runTable_Gosselin$Sample_Name)
  OlahSamples <- which(colnames(Gene_Matrix)%in%OlahBatches$Sample_Name)
  SamplesIndices <- list(Galatro = GalatroSamples,Gosselin = GosselinSamples, Olah = OlahSamples)

  
  if(!full){
    GalatroBatches <- dplyr::filter(runTable_Galatro,Sample_Name%in%colnames(Gene_Matrix)) %>%
      dplyr::select(AvgSpotLen,Sample_Name,age,gender) %>%
      dplyr::mutate(batch=ifelse(AvgSpotLen==202 | AvgSpotLen==199,'GalatroIllumina','GalatroTakara')) %>%
      dplyr::distinct(Sample_Name,.keep_all=T) %>% dplyr::select(Sample_Name,batch)
    #Based on whether expType should be considered as a batch, merge different data frames. 
    if(!expType){
      GosselinBatches <- dplyr::filter(runTable_Gosselin,Sample_Name%in%colnames(Gene_Matrix)) %>%
        dplyr::select(AvgSpotLen,Sample_Name,gender,Instrument) %>%
        dplyr::mutate(batch=ifelse(Instrument=='Illumina HiSeq 4000','GosselinIllumina','GosselinNextSeq')) %>%
        dplyr::distinct(Sample_Name,.keep_all=T) %>% dplyr::select(Sample_Name,batch)
    }else{
      GosselinBatches <- dplyr::filter(runTable_Gosselin,Sample_Name%in%colnames(Gene_Matrix)) %>%
        dplyr::select(AvgSpotLen,Sample_Name,gender,Instrument,batch = expType) %>%
        dplyr::distinct(Sample_Name,.keep_all=T) %>% dplyr::select(Sample_Name,batch)
    }
    Batches <- list(Galatro = GalatroBatches, Gosselin = GosselinBatches, Olah = OlahBatches)
    
    #Make sure batch df and samples are in same order
    SamplesIndices <- unlist(SamplesIndices[Samples])
    Batches <- do.call(rbind,Batches[Samples])
    rownames(Batches) <- Batches$Sample_Name
    Batches <- Batches[colnames(Gene_Matrix[,SamplesIndices]),]
  }else{
    GalatroBatches <- dplyr::filter(runTable_Galatro,Sample_Name%in%colnames(Gene_Matrix)) %>%
      dplyr::select(AvgSpotLen,Sample_Name,age,gender) %>%
      dplyr::mutate(batch='Galatro') %>%
      dplyr::distinct(Sample_Name,.keep_all=T) %>% dplyr::select(Sample_Name,batch)
    GosselinBatches <- dplyr::filter(runTable_Gosselin,Sample_Name%in%colnames(Gene_Matrix)) %>%
      dplyr::select(AvgSpotLen,Sample_Name,gender,Instrument) %>% dplyr::mutate(batch='Gosselin') %>%
      dplyr::distinct(Sample_Name,.keep_all=T) %>% dplyr::select(Sample_Name,batch)
    Batches <- list(Galatro = GalatroBatches, Gosselin = GosselinBatches, Olah = OlahBatches)
    SamplesIndices <- list(Galatro = GalatroSamples,Gosselin = GosselinSamples, Olah = OlahSamples)
    
    #Make sure batch df and samples are in same order
    SamplesIndices <- unlist(SamplesIndices[Samples])
    Batches <- do.call(rbind,Batches[Samples])
    rownames(Batches) <- Batches$Sample_Name
    Batches <- Batches[colnames(Gene_Matrix[,SamplesIndices]),]
  }
  #Run batch correction, on all datasets which were specified in Samples
  SalmonTPM_Gene_Combat_Merged <- BatchAdjustCombat(Gene_Matrix[,SamplesIndices],Batches)
  if(!is.null(Transcript_Matrix)){
    SalmonTPM_Transcript_Combat_Merged <- BatchAdjustCombat(Transcript_Matrix[,SamplesIndices],Batches)
  }else{
    SalmonTPM_Transcript_Combat_Merged <- NULL
  }
  return(list(SalmonTPM_Gene_Combat_Merged=SalmonTPM_Gene_Combat_Merged,
              SalmonTPM_Transcript_Combat_Merged=SalmonTPM_Transcript_Combat_Merged,
              Batches=Batches))
}

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


#Preprocess data to be passed to Combat
RunBatchCorrection <- function(GeneMatrix,TranscriptMatrix,Sample,expType=F,full=F,plots=F,figName=NULL){
  #Filter any rows with all zero counts
  GeneMatrix <- GeneMatrix[apply(GeneMatrix,1,function(x) !all(x==0)),]
  TranscriptMatrix <- TranscriptMatrix[apply(TranscriptMatrix,1,function(x) !all(x==0)),]
  
  MatrixCorrected <- RunCombat(GeneMatrix,TranscriptMatrix,Samples = Sample,expType = expType,full= full)
  PCA_PostCorrection_Gene <- CalcPCA(MatrixCorrected$SalmonTPM_Gene_Combat_Merged)
  PCA_PostCorrection_Gene$Df$batch <- MatrixCorrected$Batches$batch
  
  PCA_PostCorrection_Transcript <- CalcPCA(MatrixCorrected$SalmonTPM_Transcript_Combat_Merged)
  PCA_PostCorrection_Transcript$Df$batch <- MatrixCorrected$Batches$batch
  
  PCA_PreCorrection_Gene <- CalcPCA(GeneMatrix)
  PCA_PreCorrection_Gene$Df$batch <- MatrixCorrected$Batches$batch
  
  PCA_PreCorrection_Transcript <- CalcPCA(TranscriptMatrix)
  PCA_PreCorrection_Transcript$Df$batch <- MatrixCorrected$Batches$batch
  
  
  if(plots){
    PCA_PostCorrectionPlot <- list(autoplot(PCA_PostCorrection_Gene$PCA,data=PCA_PostCorrection_Gene$Df,colour='batch',size=4,shape=F) + ggtitle('Post Batch Correction - Gene'),
                                   autoplot(PCA_PostCorrection_Transcript$PCA,data=PCA_PostCorrection_Transcript$Df,colour='batch',size=4,shape=F) + ggtitle('Post Batch Correction - Transcript'))
    
    PCA_PreCorrectionPlot <- list(autoplot(PCA_PreCorrection_Gene$PCA,data=PCA_PreCorrection_Gene$Df,colour='batch',size=4,shape=F) + ggtitle('Pre Batch Correction - Gene'),
                                  autoplot(PCA_PreCorrection_Transcript$PCA,data=PCA_PreCorrection_Transcript$Df,colour='batch',size=4,shape=F) + ggtitle('Pre Batch Correction - Transcript'))
    library(egg)
    tiff(file=paste0('../../Figures/Batch_Correction/',figName,'.tiff'),width=800,height=500)
    ggarrange(plots = c(PCA_PreCorrectionPlot,PCA_PostCorrectionPlot),ncol=2)
    dev.off()
    return(list(PCA_Plots = list(PCA_PreCorrectionPlot,PCA_PostCorrectionPlot),MatrixCorrected=MatrixCorrected))
  }
  
  return(list(MatrixCorrected=MatrixCorrected))
}

RunBatchForAllData <- function(plots){
  source('pca_analysis.R')
  load('../../Count_Data/TPM_Filtered/TPM_Microglia_Gene_Merged.rda')
  load('../../Count_Data/TPM_Filtered/TPM_Microglia_Transcript_Merged.rda')
  GalatroSamples_Genes <- ExtractDataset(TPM_Microglia_Gene_Merged,'Galatro')
  GosselinSamples_Genes <- ExtractDataset(TPM_Microglia_Gene_Merged,'Gosselin')
  OlahSamples_Genes <- ExtractDataset(TPM_Microglia_Gene_Merged,'Olah')
  GalatroSamples_Transcripts <- ExtractDataset(TPM_Microglia_Transcript_Merged,'Galatro')
  GosselinSamples_Transcripts <- ExtractDataset(TPM_Microglia_Transcript_Merged,'Gosselin')
  OlahSamples_Transcripts <- ExtractDataset(TPM_Microglia_Transcript_Merged,'Olah')
  
  #Batch correction of Galatro lib type
  GalatroCorrected <- RunBatchCorrection(log2(GalatroSamples_Genes+1),log2(GalatroSamples_Transcripts+1),Sample='Galatro',expType = F,full=F,plots = plots,figName = 'GalatroLibCorrection')
  #Batch Correction for Gosselin ExpType
  Gosselin_ExpType_Corrected <- RunBatchCorrection(log2(GosselinSamples_Genes+1),log2(GosselinSamples_Transcripts+1),Sample='Gosselin',expType = T,full=F,plots = plots,figName = 'GosselinExpTypeCorrection')
  #Batch Correction for Gosselin Instrument Type
  Gosselin_Instrument_Corrected <- RunBatchCorrection(Gosselin_ExpType_Corrected$MatrixCorrected$SalmonTPM_Gene_Combat_Merged,
                                                      Gosselin_ExpType_Corrected$MatrixCorrected$SalmonTPM_Transcript_Combat_Merged,
                                                      Sample='Gosselin',expType = F,full=F,plots = plots,figName = 'GosselinInstrumentCorrection')
  #Merge all individual datasets after individually correction
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
  
  #Global batch correction.
  MergedFinalBatchCorrected <- RunBatchCorrection(MergedIndivCorrected$Gene,MergedIndivCorrected$Transcript,c('Galatro','Gosselin','Olah'),expType = F,full = T,plots = T,figName = 'GlobalCorrection')
  SalmonTPM_Combat_ExpCorrected <- MergedFinalBatchCorrected$MatrixCorrected
  save(SalmonTPM_Combat_ExpCorrected,file='../../Count_Data/Batch_Corrected/SalmonTPM_Combat_ExpCorrected.rda')
  save(list = ls(environment()),file='../../CodeImages/BatchCorrection.RData')
  
}





########################################PC Regress out##########################################
# PC_Correction <- function(countMatrix,metadata){
#   PC_Estimate = function(countMatrix,metadata){
#     # determine number of principal components to adjust for
#     mod <- model.matrix(~1,data=metadata)
#     n.pc <- sva::num.sv(countMatrix, mod,method='be')
#     return(n.pc)
#   }
#   ComputePCLoadings <- function(countMatrix){
#     usv <- svd(scale(t(countMatrix)))
#     return(usv$u)
#   }
#   PC_Correct <- function(countMatrix, loadings, n.pc){
#     dat <- t(countMatrix)
#     n.pc <- c(1:n.pc) #Range of PC to consider
#     print(paste("removing", n.pc, "PCs", nrow(dat)))
#     # use residuals from top n.pc principal components
#     dat.adjusted <- lm(dat ~ loadings[,n.pc])$residuals
#     return(t(dat.adjusted))
#   }
#   n.pc <- PC_Estimate(countMatrix,metadata) #Estimate # of PC to regress out
#   loadings <- ComputePCLoadings(countMatrix) #Compute loadings of each PC
#   #Retur result of batch correct PC, based on 1stPC and 1st + 2nd PC
#   return(list(OnePC=PC_Correct(countMatrix,loadings,1),TwoPC=PC_Correct(countMatrix,loadings,2)))
# }
# SalmonTPM_Gene_PCAdj <- PC_Correction(log(SalmonTPM_Gene_Filt+1),runTable)
# SalmonTPM_Transcript_PCAdj <- PC_Correction(log(SalmonTPM_Transcript_Filt+1),runTable)
# 
# #Compare PCA of methods, after batch correction. 
# 
# source('pca_analysis.R')
# No_Adj <- CalcPCA(log(SalmonTPM_Gene_Filt+1),tables)
# Gene_Combat <- CalcPCA(SalmonTPM_Gene_Combat,tables)
# Gene_PCAdj_One <- CalcPCA(SalmonTPM_Gene_PCAdj$OnePC,tables)
# Gene_PCAdj_Two <- CalcPCA(SalmonTPM_Gene_PCAdj$TwoPC,tables)
# 
# library(ggplot2)
# library(egg)
# p1 = autoplot(No_Adj$PCA, data = No_Adj$Df, colour = 'readLength',size=4,shape=F) + ggtitle('PCA - No Batch Correction')
# p2 = autoplot(Gene_Combat$PCA, data = Gene_Combat$Df, colour = 'readLength',size=4,shape=F) + ggtitle('PCA - Combat')
# p3 = autoplot(Gene_PCAdj_One$PCA, data = Gene_PCAdj_One$Df, colour = 'readLength',size=4,shape=F) + ggtitle('PCA - 1st PCA Regression')
# p4 = autoplot(Gene_PCAdj_Two$PCA, data = Gene_PCAdj_Two$Df, colour = 'readLength',size=4,shape=F) + ggtitle('PCA - 1st and 2nd PCA Regression')
# ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)

