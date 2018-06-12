################################################################################################
#Batch correction, 
#input should be count matrix, 
#output is Combat and PC corrected count matrix
################################################################################################

library(sva)
#Load count matrices, which have been filtered by TPM and CV
load(file='../Count_Data/SalmonTPM_Gene_Microglia_Filt.rda')
load(file='../Count_Data/SalmonTPM_Transcript_Microglia_Filt.rda')

########################################COMBAT##########################################
#Batch Adjustment using Combat. Input is count Matrix (where rows are genes and columns are samples)
#Metadata is data.frame, with mapping between batch and sample. 
BatchAdjustCombat <- function(countMatrix,metadata){
  #Add no other covariates to Combat batch correct
  modcombat <- model.matrix(~1,data=metadata)
  #Do batch correction using EB of Combat, based on batch as provided in metadata
  combat_edata = sva::ComBat(dat=countMatrix, batch=metadata$batch, mod=modcombat, par.prior=F,prior.plots = F)
  return(combat_edata)
}

#Load Information of each sample
runTable <- read.table(file = '/local/data/public/zmx21/zmx21_private/GSK/Galatro/SraRunTable.txt',header = T,sep = '\t')
runTable <- dplyr::filter(runTable,Sample_Name%in%colnames(SalmonTPM_Gene_Filt)) %>% 
            dplyr::select(AvgSpotLen,Sample_Name,age,gender) %>%
            dplyr::mutate(batch=ifelse(AvgSpotLen==202,1,2)) %>%
            dplyr::distinct(Sample_Name,.keep_all=T)
rownames(runTable) <- runTable$Sample_Name
runTable <- runTable[colnames(SalmonTPM_Gene_Filt),] #Same order as count matrix

#Do batch Adjustment
SalmonTPM_Gene_Combat <- BatchAdjustCombat(log(SalmonTPM_Gene_Filt+1),runTable)
SalmonTPM_Transcript_Combat <- BatchAdjustCombat(log(SalmonTPM_Transcript_Filt+1),runTable)

########################################PC REgress out##########################################
PC_Correction <- function(countMatrix,metadata){
  PC_Estimate = function(countMatrix,metadata){
    # determine number of principal components to adjust for
    mod <- model.matrix(~1,data=metadata)
    n.pc <- sva::num.sv(countMatrix, mod,method='be')
    return(n.pc)
  }
  ComputePCLoadings <- function(countMatrix){
    usv <- svd(scale(t(countMatrix)))
    return(usv$u)
  }
  PC_Correct <- function(countMatrix, loadings, n.pc){
    dat <- t(countMatrix)
    n.pc <- c(1:n.pc) #Range of PC to consider
    print(paste("removing", n.pc, "PCs", nrow(dat)))
    # use residuals from top n.pc principal components
    dat.adjusted <- lm(dat ~ loadings[,n.pc])$residuals
    return(t(dat.adjusted))
  }
  n.pc <- PC_Estimate(countMatrix,metadata) #Estimate # of PC to regress out
  loadings <- ComputePCLoadings(countMatrix) #Compute loadings of each PC
  #Retur result of batch correct PC, based on 1stPC and 1st + 2nd PC
  return(list(OnePC=PC_Correct(countMatrix,loadings,1),TwoPC=PC_Correct(countMatrix,loadings,2)))
}
SalmonTPM_Gene_PCAdj <- PC_Correction(log(SalmonTPM_Gene_Filt+1),runTable)
SalmonTPM_Transcript_PCAdj <- PC_Correction(log(SalmonTPM_Transcript_Filt+1),runTable)

#Compare PCA of methods, after batch correction. 

source('pca_analysis.R')
No_Adj <- CalcPCA(log(SalmonTPM_Gene_Filt+1),tables)
Gene_Combat <- CalcPCA(SalmonTPM_Gene_Combat,tables)
Gene_PCAdj_One <- CalcPCA(SalmonTPM_Gene_PCAdj$OnePC,tables)
Gene_PCAdj_Two <- CalcPCA(SalmonTPM_Gene_PCAdj$TwoPC,tables)

library(ggplot2)
library(egg)
p1 = autoplot(No_Adj$PCA, data = No_Adj$Df, colour = 'readLength',size=4,shape=F) + ggtitle('PCA - No Batch Correction')
p2 = autoplot(Gene_Combat$PCA, data = Gene_Combat$Df, colour = 'readLength',size=4,shape=F) + ggtitle('PCA - Combat')
p3 = autoplot(Gene_PCAdj_One$PCA, data = Gene_PCAdj_One$Df, colour = 'readLength',size=4,shape=F) + ggtitle('PCA - 1st PCA Regression')
p4 = autoplot(Gene_PCAdj_Two$PCA, data = Gene_PCAdj_Two$Df, colour = 'readLength',size=4,shape=F) + ggtitle('PCA - 1st and 2nd PCA Regression')
ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)

