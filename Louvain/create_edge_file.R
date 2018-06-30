CreateEdgeFiles <- function(plots){
  a <- Sys.time()
  library(Hmisc)
  library(parallel)
  library(plyr)

  if(!'MicrogliaGeneCodingCorr' %in% ls()){
    load('../../Count_Data/Correlation_Matrices/MicrogliaGeneCodingCorr.rda')
    # load('../../Count_Data/CV_Filtered/MicrogliaGeneCVFiltered.rda')
    # MicrogliaGeneCodingCorr <- rcorr(t(MicrogliaGeneCVFiltered$coding),type = 'pearson')
  }
  if(!'MicrogliaGeneAllCorr' %in% ls()){
    load('../../Count_Data/Correlation_Matrices/MicrogliaGeneAllCorr.rda')
    
    # MicrogliaGeneAll <- rbind(MicrogliaGeneCVFiltered$coding,MicrogliaGeneCVFiltered$noncoding)
    # MicrogliaGeneAllCorr <- rcorr(t(MicrogliaGeneAll),type = 'pearson')
  }
  if(!'BrainGeneCodingCorr' %in% ls()){
    load('../../Count_Data/Correlation_Matrices/BrainGeneCodingCorr.rda')

    # load('../../Count_Data/CV_Filtered/BrainGeneCVFiltered.rda')
    # BrainGeneCodingCorr <- rcorr(t(BrainGeneCVFiltered$coding),type = 'pearson')
  }
  if(!'BrainGeneAllCorr' %in% ls()){
    load('../../Count_Data/Correlation_Matrices/BrainGeneAllCorr.rda')

    # BrainGeneAll <- rbind(BrainGeneCVFiltered$coding,BrainGeneCVFiltered$noncoding)
    # BrainGeneAllCorr <- rcorr(t(BrainGeneAll),type = 'pearson')
  }
  pValCutOff = 0.05
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
  
  WriteFile <- function(corrList,pValCutOff,path){
    #Convert correlation matrix to a 3 columns. 1st and 2nd columns are nodes, and 3rd column is the weight.
    corrList$P[lower.tri(corrList$P,diag = T)] <- Inf #only consider upper diagonal since symmetric
    significantIndex <- which(corrList$P < pValCutOff,arr.ind = T,useNames = F)
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
  WriteFile(MicrogliaGeneCodingCorr,pValCutOff,'../../Louvain_Edge_List/CodingGenesEdgeListMicroglia.txt')
  print(Sys.time() - a)
  
  a <- Sys.time()
  #microglia all genes
  WriteFile(MicrogliaGeneAllCorr,pValCutOff,'../../Louvain_Edge_List/AllGenesEdgeListMicroglia.txt')
  print(Sys.time() -a)
}
