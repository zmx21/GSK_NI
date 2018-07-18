#use for transcript level, slower but cannot load whole matrix into RAM.
a <- Sys.time()
library(parallel)
WriteEdgeList <- function(expMatrix,pValCutOff,corCutOff,path){
  print(paste('writing: ',path))
  numGenes <- nrow(expMatrix)
  geneNames <- rownames(expMatrix)
  numSamples <- ncol(expMatrix)
  #Write data in columns of the correlation matrix, to reduce RAM usage.
  #Loop through columns, take rows larger than the column index (lower triangular)
  #Write correlation values above p-val cutoff to file.
  cor.p.values <- function(r, n) {
    df <- n - 2
    t <- sqrt(df) * r / sqrt(1 - r^2)
    p <- pt(t, df)
    return(2 * pmin(p, 1 - p))
  }
  RunCorrelationForPair <- function(expMatrix,pValCutOff,corCutOff,i,j,numSamples){
    corResult <- cor(expMatrix[i,],expMatrix[j,],method = 'pearson')
    if(pValCutOff==Inf){
      weightVal <- ifelse(abs(corResult) > corCutOff,round(abs(corResult),3),c(NA))
    }else{
      weightVal <- ifelse(abs(corResult) > corCutOff & cor.p.values(corResult,numSamples) < pValCutOff,round(abs(corResult),3),c(NA))
    }
    return(weightVal)
  }
  for(j in seq(1,numGenes-1,by=1)){
    currentRows <- seq(j+1,numGenes,1)
    correlationResults <- unlist(mclapply(currentRows,function(i) RunCorrelationForPair(expMatrix,pValCutOff,corCutOff,i,j,numSamples),mc.cores = 20))
    overCutOffRows <- is.na(correlationResults)
    correlationResults <- data.frame(Node1=rep(geneNames[j],length(currentRows[!overCutOffRows])),
                                     Node2=geneNames[currentRows[!overCutOffRows]],
                                     W = correlationResults[!overCutOffRows])
    data.table::fwrite(correlationResults,file = path,append = T,row.names = F,col.names = F,sep = ' ',quote = F)
  }
}

# a <- Sys.time()
load('../../Count_Data/CV_Filtered/MicrogliaTranscriptCVFiltered.rda')
# #Microglia coding transcripts
# WriteEdgeList(MicrogliaTranscriptCVFiltered$coding,corCutOff = 0.2,pValCutOff=Inf,'../../Louvain_Edge_List/CodingTranscriptsEdgeListMicroglia.txt')
# print(Sys.time() - a)
# 
# 
a <- Sys.time()
# #Microglia all transcripts
WriteEdgeList(rbind(MicrogliaTranscriptCVFiltered$coding,MicrogliaTranscriptCVFiltered$noncoding),corCutOff = 0.2,pValCutOff=Inf,'../../Louvain_Edge_List/AllTranscriptsEdgeListMicroglia.txt')
print(Sys.time() - a)
