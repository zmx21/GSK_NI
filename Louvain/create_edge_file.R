CreateEdgeFiles <- function(plots){
  library(Hmisc)
  library(parallel)
  library(plyr)
  load('../../Count_Data/Correlation_Matrices/MicrogliaGeneAllCorr.rda')
  load('../../Count_Data/Correlation_Matrices/MicrogliaGeneCodingCorr.rda')
  
  if(!'MicrogliaGeneCodingCorr' %in% ls()){
    load('../../Count_Data/CV_Filtered/MicrogliaGeneCVFiltered.rda')
    MicrogliaGeneCodingCorr <- rcorr(t(MicrogliaGeneCVFiltered$coding),type = 'pearson')
  }
  if(!'MicrogliaGeneAllCorr' %in% ls()){
    MicrogliaGeneAll <- rbind(MicrogliaGeneCVFiltered$coding,MicrogliaGeneCVFiltered$noncoding)
    MicrogliaGeneAllCorr <- rcorr(t(MicrogliaGeneAll),type = 'pearson')
  }
  if(!'BrainGeneCodingCorr' %in% ls()){
    load('../../Count_Data/CV_Filtered/BrainGeneCVFiltered.rda')
    BrainGeneCodingCorr <- rcorr(t(BrainGeneCVFiltered$coding),type = 'pearson')
  }
  if(!'BrainGeneAllCorr' %in% ls()){
    BrainGeneAll <- rbind(BrainGeneCVFiltered$coding,BrainGeneCVFiltered$noncoding)
    BrainGeneAllCorr <- rcorr(t(BrainGeneAll),type = 'pearson')
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
    
    #Observe distribution of CV and p-value.
    tiff(filename = '../../Figures/Correlation_Network/Corr_PVal_Dist_Brain.tiff',width = 600,height=400)
    par(mfrow=c(2,2))
    hist(as.vector(BrainGeneCodingCorr$r),xlab='Correlation Coeff',main=paste0('Coding Network -\nCorrelation',' (N=',as.character(length(as.vector(BrainGeneCodingCorr$r))),')'))
    hist(as.vector(BrainGeneCodingCorr$P),xlab='P Value',main=paste0('Coding Network -\nP-Value',' (N=',as.character(length(as.vector(BrainGeneCodingCorr$r))),')'))
    hist(as.vector(BrainGeneAllCorr$r),xlab='Correlation Coeff',main=paste0('Coding and Non Coding Network -\nCorrelation',' (N=',as.character(length(as.vector(BrainGeneAllCorr$r))),')'))
    hist(as.vector(BrainGeneAllCorr$P),xlab='P Value',main=paste0('Coding and Non Coding Network -\nP-Value',' (N=',as.character(length(as.vector(BrainGeneAllCorr$r))),')'))
    dev.off()
    
    #Distribution of Correlation values, after cutoff
    tiff(filename = '../../Figures/Correlation_Network/Corr_Filtered_Brain.tiff',width = 300,height=400)
    par(mfrow=c(2,1))
    hist(as.vector(BrainGeneCodingCorr$r[which(BrainGeneCodingCorr$P<pValCutOff)]),xlab='Correlation Coeff',
         main=paste0('Coding Network -\nCorrelation Filtered',' (N=',as.character(length(as.vector(BrainGeneCodingCorr$r[which(BrainGeneCodingCorr$P<pValCutOff)]))),')'))
    hist(as.vector(BrainGeneAllCorr$r[BrainGeneAllCorr$P<pValCutOff]),xlab='Correlation Coeff',
         main=paste0('Coding and Non Coding Network -\nCorrelation Filtered',' (N=',as.character(length(as.vector(BrainGeneAllCorr$r[which(BrainGeneAllCorr$P<pValCutOff)]))),')'))
    dev.off()
    return()
  }
  #Covert correlation matrix to a 3 columns. 1st and 2nd columns are nodes, and 3rd column is the weight.
  CodingSignificantIndex <- which(MicrogliaGeneCodingCorr$P < pValCutOff,arr.ind = T,useNames = F)
  CodingGeneNames <- rownames(MicrogliaGeneCodingCorr$r)
  codingOutFile <- file('../../Louvain_Edge_List/CodingEdgeList.txt','w')
  print('writing coding file')
  for(i in 1:nrow(CodingSignificantIndex)){
    x <- CodingSignificantIndex[i,]
    currentRow <- c(x[1],x[2],abs(round(MicrogliaGeneCodingCorr$r[x[1],x[2]],3)))
    write(currentRow,codingOutFile,append = T,ncolumns = 4)
  }
  close.connection(codingOutFile)
  save(CodingGeneNames,file='../../Louvain_Edge_List/CodingGeneNames.rda')

  #Covert correlation matrix to a 3 columns. 1st and 2nd columns are nodes, and 3rd column is the weight.
  AllSignificantIndex <- which(MicrogliaGeneAllCorr$P < pValCutOff,arr.ind = T,useNames = F)
  AllGeneNames <- rownames(MicrogliaGeneAllCorr$r)
  AllOutFile <- file('../../Louvain_Edge_List/AllEdgeList.txt','w')
  print('writing all file')
  for(i in 1:nrow(AllSignificantIndex)){
    x <- AllSignificantIndex[i,]
    currentRow <- c(x[1],x[2],abs(round(MicrogliaGeneAllCorr$r[x[1],x[2]],3)))
    write(currentRow,AllOutFile,append = T,ncolumns = 4)
  }
  close.connection(AllOutFile)
  save(AllGeneNames,file='../../Louvain_Edge_List/AllGeneNames.rda')

}

#Filter according to p-value, observe resulting distribution of correlation values.
# tiff(filename = '../../Figures/Correlation_Network/Corr_PVal_Dist.tiff',width = 600,height=400)
# par(mfrow=c(1,2))
# hist(as.vector(MicrogliaGeneCodingCorr$r),main=paste0('Filtered Coding Network -\nCorrelation',' (N=',as.character(length(as.vector(MicrogliaGeneCodingCorr$r))),')'))
# hist(as.vector(MicrogliaGeneAllCorr$r),main=paste0('Filtered Coding and Non Coding Network -\nCorrelation',' (N=',as.character(length(as.vector(MicrogliaGeneAllCorr$r))),')'))
#
# dev.off()