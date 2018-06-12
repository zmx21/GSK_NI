CalcCoeffVariation <- function(countMatrix){
  #Remove genes which are not expressed in all samples
  allZero <- apply(countMatrix,1,function(x) all(x==0))
  countMatrix <- countMatrix[!allZero,]
  coeffVar <- apply(countMatrix,1,function(x) sd(x)/mean(x))
  coeffMean <- apply(countMatrix,1,function(x) mean(x))
  
  return(list(coeffVar=coeffVar,coeffMean = coeffMean))
}
coeffVarGene <- CalcCoeffVariation(SalmonTPM$geneLevel$abundance)$coeffVar
coeffVarTranscript <- CalcCoeffVariation(SalmonTPM$transcriptLevel$abundance)$coeffVar
meanGene <- CalcCoeffVariation(SalmonTPM$geneLevel$abundance)$coeffMean
meanTranscript <- CalcCoeffVariation(SalmonTPM$transcriptLevel$abundance)$coeffMean

par(mfrow=c(1,2))
hist(coeffVarTranscript,breaks = 50,xlab = 'CV',ylab='Count',main='Transcript Level CV')
hist(coeffVarGene,breaks = 50,xlab = 'CV',ylab='Count',main='Gene Level CV')