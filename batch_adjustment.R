################################################################################################
#Batch correction, 
#input should be count matrix, expType specifies whether batch correction should be performed on 
#expType of Gosselin. Full specifies whether batch correction should be performed on dataset alone. 
#output is Combat and PC corrected count matrix
########################################COMBAT##########################################
RunCombat <- function(SalmonTPM_Gene_Merged,SalmonTPM_Transcript_Merged = NULL,Samples,expType = F,full=F){
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
                            read.table(file = '/local/data/public/zmx21/zmx21_private/GSK/Galatro_Brain//SraRunTable.txt',header = T,sep = '\t',stringsAsFactors = F))
  #No batches in Olah
  OlahBatches <- data.frame(Sample_Name = colnames(SalmonTPM_Gene_Merged)[sapply(colnames(SalmonTPM_Gene_Merged),function(x) grepl('H5KN',x))])
  OlahBatches$batch <- rep('Olah',nrow(OlahBatches))
  
  
  #Load batch Information of each sample, Galatro et al
  runTable_Gosselin <- read.table(file = '/local/data/public/zmx21/zmx21_private/GSK/Gosselin/SraRunTable_Parsed.txt',header = T,sep = '\t')
  runTable_Gosselin <- dplyr::mutate(runTable_Gosselin,expType = ifelse(grepl('ExVivo',Library_Name),'ExVivo','InVitro'))
  
  #Choose which samples to include, based on provided Samples character vector.
  GalatroSamples <- which(colnames(SalmonTPM_Gene_Merged)%in%runTable_Galatro$Sample_Name)
  GosselinSamples <- which(colnames(SalmonTPM_Gene_Merged)%in%runTable_Gosselin$Sample_Name)
  OlahSamples <- which(colnames(SalmonTPM_Gene_Merged)%in%OlahBatches$Sample_Name)
  SamplesIndices <- list(Galatro = GalatroSamples,Gosselin = GosselinSamples, Olah = OlahSamples)

  
  if(!full){
    GalatroBatches <- dplyr::filter(runTable_Galatro,Sample_Name%in%colnames(SalmonTPM_Gene_Merged)) %>%
      dplyr::select(AvgSpotLen,Sample_Name,age,gender) %>%
      dplyr::mutate(batch=ifelse(AvgSpotLen==202 | AvgSpotLen==199,'GalatroIllumina','GalatroTakara')) %>%
      dplyr::distinct(Sample_Name,.keep_all=T) %>% dplyr::select(Sample_Name,batch)
    #Based on whether expType should be considered as a batch, merge different data frames. 
    if(!expType){
      GosselinBatches <- dplyr::filter(runTable_Gosselin,Sample_Name%in%colnames(SalmonTPM_Gene_Merged)) %>%
        dplyr::select(AvgSpotLen,Sample_Name,gender,Instrument) %>%
        dplyr::mutate(batch=ifelse(Instrument=='Illumina HiSeq 4000','GosselinIllumina','GosselinNextSeq')) %>%
        dplyr::distinct(Sample_Name,.keep_all=T) %>% dplyr::select(Sample_Name,batch)
    }else{
      GosselinBatches <- dplyr::filter(runTable_Gosselin,Sample_Name%in%colnames(SalmonTPM_Gene_Merged)) %>%
        dplyr::select(AvgSpotLen,Sample_Name,gender,Instrument,batch = expType) %>%
        dplyr::distinct(Sample_Name,.keep_all=T) %>% dplyr::select(Sample_Name,batch)
    }
    Batches <- list(Galatro = GalatroBatches, Gosselin = GosselinBatches, Olah = OlahBatches)
    
    #Make sure batch df and samples are in same order
    SamplesIndices <- unlist(SamplesIndices[Samples])
    Batches <- do.call(rbind,Batches[Samples])
    rownames(Batches) <- Batches$Sample_Name
    Batches <- Batches[colnames(SalmonTPM_Gene_Merged[,SamplesIndices]),]
  }else{
    GalatroBatches <- dplyr::filter(runTable_Galatro,Sample_Name%in%colnames(SalmonTPM_Gene_Merged)) %>%
      dplyr::select(AvgSpotLen,Sample_Name,age,gender) %>%
      dplyr::mutate(batch='Galatro') %>%
      dplyr::distinct(Sample_Name,.keep_all=T) %>% dplyr::select(Sample_Name,batch)
    GosselinBatches <- dplyr::filter(runTable_Gosselin,Sample_Name%in%colnames(SalmonTPM_Gene_Merged)) %>%
      dplyr::select(AvgSpotLen,Sample_Name,gender,Instrument) %>% dplyr::mutate(batch='Gosselin') %>%
      dplyr::distinct(Sample_Name,.keep_all=T) %>% dplyr::select(Sample_Name,batch)
    Batches <- list(Galatro = GalatroBatches, Gosselin = GosselinBatches, Olah = OlahBatches)
    SamplesIndices <- list(Galatro = GalatroSamples,Gosselin = GosselinSamples, Olah = OlahSamples)
    
    #Make sure batch df and samples are in same order
    SamplesIndices <- unlist(SamplesIndices[Samples])
    Batches <- do.call(rbind,Batches[Samples])
    rownames(Batches) <- Batches$Sample_Name
    Batches <- Batches[colnames(SalmonTPM_Gene_Merged[,SamplesIndices]),]
  }
  #Run batch correction, on all datasets which were specified in Samples
  SalmonTPM_Gene_Combat_Merged <- BatchAdjustCombat(SalmonTPM_Gene_Merged[,SamplesIndices],Batches)
  if(!is.null(SalmonTPM_Transcript_Merged)){
    SalmonTPM_Transcript_Combat_Merged <- BatchAdjustCombat(SalmonTPM_Transcript_Merged[,SamplesIndices],Batches)
  }else{
    SalmonTPM_Transcript_Combat_Merged <- NULL
  }
  return(list(SalmonTPM_Gene_Combat_Merged=SalmonTPM_Gene_Combat_Merged,
              SalmonTPM_Transcript_Combat_Merged=SalmonTPM_Transcript_Combat_Merged,
              Batches=Batches))
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

