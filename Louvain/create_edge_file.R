CreateEdgeFiles <- function(plots=F){
  a <- Sys.time()
  library(Hmisc)
  library(parallel)
  library(plyr)

  if(!'MicrogliaGeneCodingCorr' %in% ls()){
    # load('../../Count_Data/CV_Filtered/MicrogliaGeneCVFiltered.rda')
    # MicrogliaGeneCodingCorr <- rcorr(t(MicrogliaGeneCVFiltered$coding),type = 'pearson')
    print('loading coding corr')
    load('../../Count_Data/Correlation_Matrices/MicrogliaGeneCodingCorr.rda')
    
  }
  if(!'MicrogliaGeneAllCorr' %in% ls()){

    # MicrogliaGeneAll <- rbind(MicrogliaGeneCVFiltered$coding,MicrogliaGeneCVFiltered$noncoding)
    # MicrogliaGeneAllCorr <- rcorr(t(MicrogliaGeneAll),type = 'pearson')
    
    print('loading all corr')
    load('../../Count_Data/Correlation_Matrices/MicrogliaGeneAllCorr.rda')
  }
  if(!'RandomMicrogliaGeneCodingCorr' %in% ls()){
    #For random test, shuffle node labels while keeping the same network. 
    # load('../../Count_Data/CV_Filtered/MicrogliaGeneCVFiltered.rda')
    # set.seed(199)
    # permutation <- sample(1:nrow(MicrogliaGeneCVFiltered$coding),size = nrow(MicrogliaGeneCVFiltered$coding),replace = F)
    # RandomMicrogliaGeneCodingCorr <- rcorr(t(MicrogliaGeneCVFiltered$coding),type = 'pearson')
    # colnames(RandomMicrogliaGeneCodingCorr$P) <- colnames(RandomMicrogliaGeneCodingCorr$P)[permutation]
    # colnames(RandomMicrogliaGeneCodingCorr$r) <- colnames(RandomMicrogliaGeneCodingCorr$r)[permutation]
    # rownames(RandomMicrogliaGeneCodingCorr$P) <- rownames(RandomMicrogliaGeneCodingCorr$P)[permutation]
    # rownames(RandomMicrogliaGeneCodingCorr$r) <- rownames(RandomMicrogliaGeneCodingCorr$r)[permutation]
  }
  if(!'BrainGeneCodingCorr' %in% ls()){
    # load('../../Count_Data/Correlation_Matrices/BrainGeneCodingCorr.rda')

    # load('../../Count_Data/CV_Filtered/BrainGeneCVFiltered.rda')
    # BrainGeneCodingCorr <- rcorr(t(BrainGeneCVFiltered$coding),type = 'pearson')
  }
  if(!'BrainGeneAllCorr' %in% ls()){
    # load('../../Count_Data/Correlation_Matrices/BrainGeneAllCorr.rda')

    # BrainGeneAll <- rbind(BrainGeneCVFiltered$coding,BrainGeneCVFiltered$noncoding)
    # BrainGeneAllCorr <- rcorr(t(BrainGeneAll),type = 'pearson')
  }
  if(plots){
    #Observe distribution of CV and p-value. Microglia
    tiff(filename = '../../Figures/Correlation_Network/Corr_PVal_Dist.tiff',width = 600,height=400)
    par(mfrow=c(2,2))
    hist(as.vector(MicrogliaGeneCodingCorr$r),xlab='Correlation Coeff',main=paste0('Coding Network -\nCorrelation',' (N=',as.character(length(as.vector(MicrogliaGeneCodingCorr$r))),')'))
    hist(as.vector(MicrogliaGeneCodingCorr$P),xlab='P Value',main=paste0('Coding Network -\nP-Value',' (N=',as.character(length(as.vector(MicrogliaGeneCodingCorr$r))),')'))
    hist(as.vector(MicrogliaGeneAllCorr$r),xlab='Correlation Coeff',main=paste0('Coding and Non Coding Network -\nCorrelation',' (N=',as.character(length(as.vector(MicrogliaGeneAllCorr$r))),')'))
    hist(as.vector(MicrogliaGeneAllCorr$P),xlab='P Value',main=paste0('Coding and Non Coding Network -\nP-Value',' (N=',as.character(length(as.vector(MicrogliaGeneAllCorr$r))),')'))
    dev.off()
    
    #Distribution of Correlation values, after cutoff. Microglia
    tiff(filename = '../../Figures/Correlation_Network/Corr_Filtered.tiff',width = 300,height=400)
    par(mfrow=c(2,1))
    hist(as.vector(MicrogliaGeneCodingCorr$r[which(MicrogliaGeneCodingCorr$P<pValCutOff)]),xlab='Correlation Coeff',
         main=paste0('Coding Network -\nCorrelation Filtered',' (N=',as.character(length(as.vector(MicrogliaGeneCodingCorr$r[which(MicrogliaGeneCodingCorr$P<pValCutOff)]))),')'))
    hist(as.vector(MicrogliaGeneAllCorr$r[MicrogliaGeneAllCorr$P<pValCutOff]),xlab='Correlation Coeff',
         main=paste0('Coding and Non Coding Network -\nCorrelation Filtered',' (N=',as.character(length(as.vector(MicrogliaGeneAllCorr$r[which(MicrogliaGeneAllCorr$P<pValCutOff)]))),')'))
    dev.off()
    return()
    
    #Observe distribution of CV and p-value. Brain
    tiff(filename = '../../Figures/Correlation_Network/Corr_PVal_Dist_Brain.tiff',width = 600,height=400)
    par(mfrow=c(2,2))
    hist(as.vector(BrainGeneCodingCorr$r),xlab='Correlation Coeff',main=paste0('Coding Network -\nCorrelation',' (N=',as.character(length(as.vector(BrainGeneCodingCorr$r))),')'))
    hist(as.vector(BrainGeneCodingCorr$P),xlab='P Value',main=paste0('Coding Network -\nP-Value',' (N=',as.character(length(as.vector(BrainGeneCodingCorr$r))),')'))
    hist(as.vector(BrainGeneAllCorr$r),xlab='Correlation Coeff',main=paste0('Coding and Non Coding Network -\nCorrelation',' (N=',as.character(length(as.vector(BrainGeneAllCorr$r))),')'))
    hist(as.vector(BrainGeneAllCorr$P),xlab='P Value',main=paste0('Coding and Non Coding Network -\nP-Value',' (N=',as.character(length(as.vector(BrainGeneAllCorr$r))),')'))
    dev.off()
    
    #Distribution of Correlation values, after cutoff. Brain
    tiff(filename = '../../Figures/Correlation_Network/Corr_Filtered_Brain.tiff',width = 300,height=400)
    par(mfrow=c(2,1))
    hist(as.vector(BrainGeneCodingCorr$r[which(BrainGeneCodingCorr$P<pValCutOff)]),xlab='Correlation Coeff',
         main=paste0('Coding Network -\nCorrelation Filtered',' (N=',as.character(length(as.vector(BrainGeneCodingCorr$r[which(BrainGeneCodingCorr$P<pValCutOff)]))),')'))
    hist(as.vector(BrainGeneAllCorr$r[BrainGeneAllCorr$P<pValCutOff]),xlab='Correlation Coeff',
         main=paste0('Coding and Non Coding Network -\nCorrelation Filtered',' (N=',as.character(length(as.vector(BrainGeneAllCorr$r[which(BrainGeneAllCorr$P<pValCutOff)]))),')'))
    dev.off()
    return()
  }
  
  WriteFile <- function(corrList,corCutOff,pValCutOff,absoluteValue,path){
    #P value multiple correction adjustment, only taking upper diagonal
    adjPValue <- p.adjust(corrList$P[upper.tri(corrList$P,diag = F)],method = 'BH')
    adjMatrixFullP <- matrix(Inf,nrow = nrow(corrList$P),ncol = ncol(corrList$P))
    adjMatrixFullP[upper.tri(adjMatrixFullP,diag = F)] <- adjPValue

    if(absoluteValue){
      adjMatrixFull <- ifelse(abs(corrList$r)>corCutOff & adjMatrixFullP < pValCutOff,1,0)
    }else{
      adjMatrixFull <- ifelse(corrList$r>corCutOff & adjMatrixFullP < pValCutOff,1,0)
    }
    adjMatrixFull[lower.tri(adjMatrixFull,diag = T)] <- 0 #only consider top triangular
    significantIndex <- which(adjMatrixFull == 1,arr.ind = T,useNames = F)
    geneNames<- rownames(corrList$r)
    print(paste('writing: ',path))
    outFile <- file(path,'w')
    for(i in 1:nrow(significantIndex)){
      x <- significantIndex[i,]
      currentRow <- c(geneNames[x[1]],
                      geneNames[x[2]],
                      as.character(abs(round(corrList$r[x[1],x[2]],3))))
      write(currentRow,outFile,append = T,ncolumns = 4)
    }
    close.connection(outFile)
  }
  #microglia coding genes
  WriteFile(MicrogliaGeneCodingCorr,corCutOff=0.25,pValCutOff = 0.05,absoluteValue=T,'../../Louvain_Edge_List/Pearson/CodingGenesEdgeListMicroglia_Pearson_pval0p05_cor0p25_abs.txt')
  
  #microglia all genes
  WriteFile(MicrogliaGeneAllCorr,corCutOff=0.25,pValCutOff = 0.05,absoluteValue=T,'../../Louvain_Edge_List/Pearson/AllGenesEdgeListMicroglia_Pearson_pval0p05_cor0p25_abs.txt')
  
  #microglia coding genes
  WriteFile(MicrogliaGeneCodingCorr,corCutOff=0.25,pValCutOff = 0.05,absoluteValue=F,'../../Louvain_Edge_List/Pearson/CodingGenesEdgeListMicroglia_Pearson_pval0p05_cor0p25_noabs.txt')
  
  #microglia all genes
  WriteFile(MicrogliaGeneAllCorr,corCutOff=0.25,pValCutOff = 0.05,absoluteValue=F,'../../Louvain_Edge_List/Pearson/AllGenesEdgeListMicroglia_Pearson_pval0p05_cor0p25_noabs.txt')
  
  #microglia coding genes
  sortedCor <- sort(MicrogliaGeneCodingCorr$r[upper.tri(MicrogliaGeneCodingCorr$r,diag = F)],decreasing=F)
  top1milCutOff <- sortedCor[length(sortedCor) - 1e6]
  WriteFile(MicrogliaGeneCodingCorr,corCutOff=top1milCutOff,pValCutOff = Inf,absoluteValue=F,'../../Louvain_Edge_List/Pearson/CodingGenesEdgeListMicroglia_Pearson_top1milpos.txt')
  
  #microglia all genes
  sortedCor <- sort(MicrogliaGeneAllCorr$r[upper.tri(MicrogliaGeneAllCorr$r,diag = F)],decreasing=F)
  top1milCutOff <- sortedCor[length(sortedCor) - 1e6]
  WriteFile(MicrogliaGeneAllCorr,corCutOff=top1milCutOff,pValCutOff = Inf,absoluteValue=F,'../../Louvain_Edge_List/Pearson/AllGenesEdgeListMicroglia_Pearson_top1milpos.txt')
  
}
