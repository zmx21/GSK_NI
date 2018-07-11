#use for transcript level, slower but cannot load whole matrix into RAM.
a <- Sys.time()
library(parallel)
#Construct edge list, which is 3 columns. 1st and 2nd columns are nodes, and 3rd column is the weight.
WriteEdgeList <- function(expMatrix,pValCutOff,path){
  print(paste('writing: ',path))
  numGenes <- nrow(expMatrix)
  geneNames <- rownames(expMatrix)
  numCores <- min(c(numGenes/2000,40))
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
  RunCorrelationForPair <- function(expMatrix,pValCutOff,i,j,numSamples){
    corResult <- cor(expMatrix[i,],expMatrix[j,],method = 'pearson')
    weightVal <- ifelse(cor.p.values(corResult,numSamples) < pValCutOff,round(abs(corResult),3),c(NA))
    return(weightVal)
  }
  for(j in seq(1,numGenes-1,by=1)){
    currentRows <- seq(j+1,numGenes,1)
    correlationResults <- unlist(mclapply(currentRows,function(i) RunCorrelationForPair(expMatrix,pValCutOff,i,j,numSamples),mc.cores = numCores))
    overCutOffRows <- is.na(correlationResults)
    correlationResults <- data.frame(Node1=rep(geneNames[j],length(currentRows[!overCutOffRows])),
                                     Node2=geneNames[currentRows[!overCutOffRows]],
                                     W = correlationResults[!overCutOffRows])
    data.table::fwrite(correlationResults,file = path,append = T,row.names = F,col.names = F,sep = ' ',quote = F)
  }
}

pValCutOff = 0.05
# a <- Sys.time()
load('../../Count_Data/CV_Filtered/MicrogliaTranscriptCVFiltered.rda')
# #Microglia coding transcripts
WriteEdgeList(MicrogliaTranscriptCVFiltered$coding,pValCutOff,'../../Louvain_Edge_List/CodingTranscriptsEdgeListMicroglia.txt')
print(Sys.time() - a)
# 
# 
# a <- Sys.time()
# #Microglia all transcripts
WriteEdgeList(rbind(MicrogliaTranscriptCVFiltered$coding,MicrogliaTranscriptCVFiltered$noncoding),pValCutOff,'../../Louvain_Edge_List/AllTranscriptsEdgeListMicroglia.txt')
print(Sys.time() - a)
