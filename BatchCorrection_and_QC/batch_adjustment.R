##########################################################################################################
#Functions for batch correction, for within study and inter study biases.
##########################################################################################################

#Subfunction to run batch correction
#input should be count matrix
#expType specifies whether batch correction should be performed on expType of Gosselin. 
#Full specifies whether batch correction should be performed on dataset alone. 
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
  runTable_Galatro <- rbind(read.table(file = '../../Galatro/SraRunTable.txt',header = T,sep = '\t',stringsAsFactors = F))
  #No batches in Olah
  OlahBatches <- data.frame(Sample_Name = colnames(Gene_Matrix)[sapply(colnames(Gene_Matrix),function(x) grepl('H5KN',x))])
  OlahBatches$batch <- rep('Olah',nrow(OlahBatches))
  
  
  #Load batch Information of each sample, Galatro et al
  runTable_Gosselin <- read.table(file = '../../Gosselin/SraRunTable_Parsed.txt',header = T,sep = '\t')
  runTable_Gosselin <- dplyr::mutate(runTable_Gosselin,expType = ifelse(grepl('ExVivo',Library_Name),'ExVivo','InVitro'))
  
  #Choose which samples to include, based on provided Samples character vector.
  GalatroSamples <- which(colnames(Gene_Matrix)%in%runTable_Galatro$Sample_Name)
  GosselinSamples <- which(colnames(Gene_Matrix)%in%runTable_Gosselin$Sample_Name)
  OlahSamples <- which(colnames(Gene_Matrix)%in%OlahBatches$Sample_Name)
  SamplesIndices <- list(Galatro = GalatroSamples,Gosselin = GosselinSamples, Olah = OlahSamples)

  
  if(!full){
    GalatroBatches <- dplyr::filter(runTable_Galatro,Sample_Name%in%colnames(Gene_Matrix)) %>%
      dplyr::select(AvgSpotLen,Sample_Name,age,gender) %>%
      dplyr::mutate(batch=ifelse(AvgSpotLen==202 | AvgSpotLen==199,'Illumina\nLibrary','Takara\nLibrary')) %>%
      dplyr::distinct(Sample_Name,.keep_all=T) %>% dplyr::select(Sample_Name,batch)
    #Based on whether expType should be considered as a batch, merge different data frames. 
    if(!expType){
      GosselinBatches <- dplyr::filter(runTable_Gosselin,Sample_Name%in%colnames(Gene_Matrix)) %>%
        dplyr::select(AvgSpotLen,Sample_Name,gender,Instrument) %>%
        dplyr::mutate(batch=ifelse(Instrument=='Illumina HiSeq 4000','Illumina\nInstrument','NextSeq\nInstrument')) %>%
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
  
  #Run batch correciton, get pre and post correction PCA 
  MatrixCorrected <- RunCombat(GeneMatrix,TranscriptMatrix,Samples = Sample,expType = expType,full= full)
  PCA_PostCorrection_Gene <- CalcPCA(MatrixCorrected$SalmonTPM_Gene_Combat_Merged)
  PCA_PostCorrection_Gene$Df$batch <- MatrixCorrected$Batches$batch
  
  PCA_PostCorrection_Transcript <- CalcPCA(MatrixCorrected$SalmonTPM_Transcript_Combat_Merged)
  PCA_PostCorrection_Transcript$Df$batch <- MatrixCorrected$Batches$batch
  
  PCA_PreCorrection_Gene <- CalcPCA(GeneMatrix)
  PCA_PreCorrection_Gene$Df$batch <- MatrixCorrected$Batches$batch
  
  PCA_PreCorrection_Transcript <- CalcPCA(TranscriptMatrix)
  PCA_PreCorrection_Transcript$Df$batch <- MatrixCorrected$Batches$batch
  
  #PCA Plots
  if(plots){
    library(ggfortify)
    PCA_PostCorrection_Gene_Sum <- summary(PCA_PostCorrection_Gene$PCA)$importance[2,1:2]
    PCA_PostCorrection_Transcript_Sum <- summary(PCA_PostCorrection_Transcript$PCA)$importance[2,1:2]
    
    PCA_PostCorrectionPlot <- list(ggplot2::autoplot(PCA_PostCorrection_Gene$PCA,data=PCA_PostCorrection_Gene$Df,colour='batch',size=4,shape=F) + ggtitle('Post Batch Correction - Gene') +
                                     xlab(paste0('PC1 (',round(PCA_PostCorrection_Gene_Sum[1]*100,1),'%)')) + labs(color='Batch')+theme(legend.spacing.y = unit(2,'cm')) +
                                            ylab(paste0('PC2 (',round(PCA_PostCorrection_Gene_Sum[2]*100,1),'%)')),
                                   ggplot2::autoplot(PCA_PostCorrection_Transcript$PCA,data=PCA_PostCorrection_Transcript$Df,colour='batch',size=4,shape=F) + ggtitle('Post Batch Correction - Transcript')+
                                     xlab(paste0('PC1 (',round(PCA_PostCorrection_Transcript_Sum[1]*100,1),'%)')) + labs(color='Batch')+theme(legend.spacing.y = unit(2,'cm')) +
                                     ylab(paste0('PC2 (',round(PCA_PostCorrection_Transcript_Sum[2]*100,1),'%)')))
    PCA_PreCorrection_Gene_Sum <- summary(PCA_PreCorrection_Gene$PCA)$importance[2,1:2]
    PCA_PreCorrection_Transcript_Sum <- summary(PCA_PreCorrection_Transcript$PCA)$importance[2,1:2]
    
    PCA_PreCorrectionPlot <- list(ggplot2::autoplot(PCA_PreCorrection_Gene$PCA,data=PCA_PreCorrection_Gene$Df,colour='batch',size=4,shape=F) + ggtitle('Pre Batch Correction - Gene')+
                                    xlab(paste0('PC1 (',round(PCA_PreCorrection_Gene_Sum[1]*100,1),'%)')) + labs(color='Batch')+theme(legend.spacing.y = unit(2,'cm')) +
                                    ylab(paste0('PC2 (',round(PCA_PreCorrection_Gene_Sum[2]*100,1),'%)')),
                                  ggplot2::autoplot(PCA_PreCorrection_Transcript$PCA,data=PCA_PreCorrection_Transcript$Df,colour='batch',size=4,shape=F) + ggtitle('Pre Batch Correction - Transcript')+
                                    xlab(paste0('PC1 (',round(PCA_PreCorrection_Transcript_Sum[1]*100,1),'%)')) + labs(color='Batch')+ theme(legend.spacing.y = unit(2,'cm')) +
                                    ylab(paste0('PC2 (',round(PCA_PreCorrection_Transcript_Sum[2]*100,1),'%)')))
    library(egg)
    tiff(file=paste0('../../Figures/Batch_Correction/',figName,'.tiff'),width=800,height=500)
    egg::ggarrange(plots = c(PCA_PreCorrectionPlot,PCA_PostCorrectionPlot),ncol=2)
    dev.off()
    
    
    library(ggpubr)
    ggpubr::ggexport(ggpubr::ggarrange(plotlist = list(PCA_PreCorrectionPlot[[1]] + ggtitle('A)'),
                                            PCA_PostCorrectionPlot[[1]] + ggtitle('B)'),
                                            PCA_PreCorrectionPlot[[2]] + ggtitle('C)'),
                                            PCA_PostCorrectionPlot[[2]] + ggtitle('D)'))),
                     filename = paste0('../../FinalFigures/BatchCorrection/',figName,'.pdf'),height = 6,width = 10)
    return(list(PCA_Plots = list(PCA_PreCorrectionPlot,PCA_PostCorrectionPlot),MatrixCorrected=MatrixCorrected))
  }
  
  return(list(MatrixCorrected=MatrixCorrected))
}

#Main function, to run batch correction for within and between study.  
RunBatchForAllData <- function(plots=F){
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
  MergedFinalBatchCorrected <- RunBatchCorrection(MergedIndivCorrected$Gene,MergedIndivCorrected$Transcript,c('Galatro','Gosselin','Olah'),expType = F,full = T,plots = plots,figName = 'GlobalCorrection')
  SalmonTPM_Combat_ExpCorrected <- MergedFinalBatchCorrected$MatrixCorrected
  
  GalatroBatches <- GalatroCorrected$MatrixCorrected$Batches
  GosselinBatches <- dplyr::left_join(Gosselin_ExpType_Corrected$MatrixCorrected$Batches,
                                      Gosselin_Instrument_Corrected$MatrixCorrected$Batches,by=c('Sample_Name'='Sample_Name'))
  GosselinBatches$DetailedBatch <- sapply(1:nrow(GosselinBatches),function(i) paste(GosselinBatches$batch.x[i],GosselinBatches$batch.y[i],sep = '\n'))
  GosselinBatches$DetailedBatch <- sapply(GosselinBatches$DetailedBatch,function(x) gsub(x=x,pattern = "\nInstrument",replacement = ""))
  GosselinBatches <- GosselinBatches %>% dplyr::select(Sample_Name,batch=DetailedBatch)
  
  #Same batch in Olah dataset
  OlahBatches <- data.frame(Sample_Name=colnames(OlahSamples_Genes),batch=rep('Olah',ncol(OlahSamples_Genes)))
  
  #Use short form sample names. 
  colnames(TPM_Microglia_Gene_Merged) <-  sapply(colnames(TPM_Microglia_Gene_Merged),function(x) ifelse(grepl(x=x,pattern = 'H5'),paste0('OLA_',substr(x,11,12),'_',
                                                                                                                                      substr(x,20,22)),
                                                                                                                 ifelse(grepl(x=x,pattern = 'GSM'),
                                                                                                                        paste0('GAL',substr(x,9,10)),
                                                                                                                        paste0('GOS',substr(x,9,10)))))
  allDetailedBatches <- rbind(GalatroBatches,GosselinBatches,OlahBatches)
  allDetailedBatches$Sample_Name <- sapply(allDetailedBatches$Sample_Name,function(x) ifelse(grepl(x=x,pattern = 'H5'),paste0('OLA_',substr(x,11,12),'_',
                                                                                                                       substr(x,20,22)),
                                                                                                  ifelse(grepl(x=x,pattern = 'GSM'),
                                                                                                         paste0('GAL',substr(x,9,10)),
                                                                                                         paste0('GOS',substr(x,9,10)))))
  
  studyNameBatch <- SalmonTPM_Combat_ExpCorrected$Batches
  studyNameBatch$Sample_Name <- sapply(studyNameBatch$Sample_Name,function(x) ifelse(grepl(x=x,pattern = 'H5'),paste0('OLA_',substr(x,11,12),'_',
                                                                                                                          substr(x,20,22)),
                                                                                         ifelse(grepl(x=x,pattern = 'GSM'),
                                                                                                paste0('GAL',substr(x,9,10)),
                                                                                                paste0('GOS',substr(x,9,10)))))
  pre_batch <- stack(as.data.frame(TPM_Microglia_Gene_Merged)) %>% dplyr::left_join(allDetailedBatches,by=c('ind'='Sample_Name')) %>%
    dplyr::left_join(studyNameBatch %>% dplyr::select(Sample_Name,StudyName=batch),by=c('ind'='Sample_Name')) 
  pre_batch <- pre_batch %>% dplyr::arrange(batch)
  
  #TPM plots at batch correction steps
  library(forcats)
  p1 <- ggplot(pre_batch) +
    geom_boxplot(aes(x = fct_inorder(ind), y = values,fill=batch),outlier.shape = NA) +
    scale_x_discrete(name = "Sample",labels = wrap_format(5)) +
    scale_y_continuous(name = "TPM",limits = quantile(pre_batch$values, c(0.1, 0.9))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.6,size =8)) +
    labs(fill='Batch') + facet_grid(.~StudyName,scales = "free", space = "free") + ggtitle('A)')
  colnames(MergedIndivCorrected$Gene) <- sapply(colnames(MergedIndivCorrected$Gene),function(x) ifelse(grepl(x=x,pattern = 'H5'),paste0('OLA_',substr(x,11,12),'_',
                                                                                                                                        substr(x,20,22)),
                                                                                                       ifelse(grepl(x=x,pattern = 'GSM'),
                                                                                                              paste0('GAL',substr(x,9,10)),
                                                                                                              paste0('GOS',substr(x,9,10)))))
  post_indiv_batch <- stack(as.data.frame(MergedIndivCorrected$Gene)) %>% dplyr::left_join(allDetailedBatches,by=c('ind'='Sample_Name')) %>%
    dplyr::left_join(studyNameBatch %>% dplyr::select(Sample_Name,StudyName=batch),by=c('ind'='Sample_Name')) 
  post_indiv_batch <- post_indiv_batch %>% dplyr::arrange(batch)
  
  p2 <- ggplot(post_indiv_batch) +
    geom_boxplot(aes(x = fct_inorder(ind), y = values,fill=batch),outlier.shape = NA) +
    scale_x_discrete(name = "Sample",labels = wrap_format(5)) +
    scale_y_continuous(name = "Corrected Abundance",limits = quantile(post_indiv_batch$values, c(0.1, 0.9))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.6,size =8)) +
    labs(fill='Batch') + facet_grid(.~StudyName,scales = "free", space = "free") + ggtitle('B)')
  
  
  colnames(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged) <- sapply(colnames(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged),function(x) ifelse(grepl(x=x,pattern = 'H5'),paste0('OLA_',substr(x,11,12),'_',
                                                                                                                                            substr(x,20,22)),
                                                                                                           ifelse(grepl(x=x,pattern = 'GSM'),
                                                                                                                  paste0('GAL',substr(x,9,10)),
                                                                                                                  paste0('GOS',substr(x,9,10)))))
  post_batch <- stack(as.data.frame(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged)) %>% dplyr::left_join(allDetailedBatches,by=c('ind'='Sample_Name')) %>%
    dplyr::left_join(studyNameBatch %>% dplyr::select(Sample_Name,StudyName=batch),by=c('ind'='Sample_Name')) 
  post_batch <- post_batch %>% dplyr::arrange(batch)
  p3 <- ggplot(post_batch) +
    geom_boxplot(aes(x = fct_inorder(ind), y = values,fill=batch),outlier.shape = NA) +
    scale_x_discrete(name = "Sample",labels = wrap_format(5)) +
    scale_y_continuous(name = "Corrected Abundance",limits = quantile(post_batch$values, c(0.1, 0.9))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.6,size =8)) +
    labs(fill='Batch') + facet_grid(.~StudyName,scales = "free", space = "free") + ggtitle('C)')
  
  ggpubr::ggexport(ggpubr::ggarrange(p1),filename = '../../FinalFigures/Supplementary/TPM_PreBatch.pdf',width = 12,height = 5)
  ggpubr::ggexport(ggpubr::ggarrange(p2),filename = '../../FinalFigures/Supplementary/TPM_Batch.pdf',width = 12,height = 5)
  ggpubr::ggexport(ggpubr::ggarrange(p3),filename = '../../FinalFigures/Supplementary/TPM_PostBatch.pdf',width = 12,height = 5)
  
}