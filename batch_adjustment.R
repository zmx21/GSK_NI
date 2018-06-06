library(sva)
load(file='../Count_Data/SalmonTPM_Gene_Microglia_Filt.rda')
load(file='../Count_Data/SalmonTPM_Transcript_Microglia_Filt.rda')


BatchAdjustCombat <- function(countMatrix,metadata){
  modcombat <- model.matrix(~1,data=metadata)
  combat_edata = sva::ComBat(dat=countMatrix, batch=metadata$batch, mod=modcombat, par.prior=T,prior.plots = F)
  return(combat_edata)
}


runTable <- read.table(file = '/local/data/public/zmx21/zmx21_private/GSK/Galatro/SraRunTable.txt',header = T,sep = '\t')
runTable <- dplyr::filter(runTable,Sample_Name%in%colnames(SalmonTPM_Gene_Filt)) %>% 
            dplyr::select(AvgSpotLen,Sample_Name,age,gender) %>%
            dplyr::mutate(batch=ifelse(AvgSpotLen==202,1,2)) %>%
            dplyr::distinct(Sample_Name,.keep_all=T)
rownames(runTable) <- runTable$Sample_Name
runTable <- runTable[colnames(SalmonTPM_Gene_Filt),] #Same order as count matrix

SalmonTPM_Gene_Combat <- BatchAdjustCombat(log(SalmonTPM_Gene_Filt+1),runTable)
# SalmonTPM_Transcript_Combat <- BatchAdjustCombat(SalmonTPM_Transcript_Filt,runTable)

PC_Correction <- function(countMatrix,metadata){
  PC_Estimate = function(countMatrix,metadata){
    ## determine number of principal components to adjust for
    mod <- model.matrix(~1,data=metadata)
    n.pc <- sva::num.sv(countMatrix, mod,method="be")
    return(n.pc)
  }
  ComputePCLoadings <- function(countMatrix){
    usv <- svd(scale(t(countMatrix)))
    return(usv$u)
  }
  PC_Correct <- function(countMatrix, loadings, n.pc){
    dat <- t(countMatrix)
    n.pc <- c(1:n.pc)
    print(paste("removing", n.pc, "PCs", nrow(dat)))
    ## use residuals from top n.pc principal components
    dat.adjusted <- lm(dat ~ loadings[,n.pc])$residuals
    return(t(dat.adjusted))
  }
  n.pc <- PC_Estimate(countMatrix,metadata)
  loadings <- ComputePCLoadings(countMatrix)
  return(list(OnePC=PC_Correct(countMatrix,loadings,1),AllPC=PC_Correct(countMatrix,loadings,2)))
}
SalmonTPM_Gene_PCAdj <- PC_Correction(log(SalmonTPM_Gene_Filt+1),runTable)
# SalmonTPM_Transcript_PCAdj <- PC_Correction(SalmonTPM_Transcript_Filt,runTable)

source('pca_analysis.R')
No_Adj <- CalcPCA(log(SalmonTPM_Gene_Filt+1),tables)
Gene_Combat <- CalcPCA(SalmonTPM_Gene_Combat,tables)
Gene_PCAdj_One <- CalcPCA(SalmonTPM_Gene_PCAdj$OnePC,tables)
Gene_PCAdj_Two <- CalcPCA(SalmonTPM_Gene_PCAdj$AllPC,tables)


library(ggplot2)
library(egg)
p1 = autoplot(No_Adj$PCA, data = No_Adj$Df, colour = 'readLength',size=4,shape=F) + ggtitle('PCA - No Batch Correction') 
p2 = autoplot(Gene_Combat$PCA, data = Gene_Combat$Df, colour = 'readLength',size=4,shape=F) + ggtitle('PCA - Combat') 
p3 = autoplot(Gene_PCAdj_One$PCA, data = Gene_PCAdj_One$Df, colour = 'readLength',size=4,shape=F) + ggtitle('PCA - First PCA Regression') 
p4 = autoplot(Gene_PCAdj_Two$PCA, data = Gene_PCAdj_Two$Df, colour = 'readLength',size=4,shape=F) + ggtitle('PCA - First and Second PCA Regression') 
ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)

