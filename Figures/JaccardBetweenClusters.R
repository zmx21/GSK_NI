##########################################################################################################
# Fig 3.2.2.1 
##########################################################################################################

library(dplyr)
library(parallel)
library(dplyr)
load('../../Count_Data/Final_Results/allPASCALResults.rda')
allPairwiseComparison <- expand.grid(Biotype=unique(allPASCALResults$Biotype),
                               AggregationType=unique(allPASCALResults$AggregationType),
                               Method=unique(allPASCALResults$Method)) %>% dplyr::filter(!(Biotype=='all' & AggregationType=='transcript'))
CalcJaccard <- function(Gene1,Gene2){
  return(length(intersect(Gene1,Gene2))/length(union(Gene1,Gene2)))
}
#Find maximum jaccard index for each combination
FindJaccardPairs <- function(df1,df2){
  dfList <- list(df1,df2)
  #df with least number of clusters
  minDf <- dfList[[which.min(sapply(dfList,nrow))]]
  maxDf <- dfList[[which.max(sapply(dfList,nrow))]]
  
  minDf <- minDf[sample(1:nrow(minDf),size = min(c(nrow(minDf),5000)),replace = F),]
  #Calculate Jaccard index in parallel, for each module 
  jaccardIndex <- unlist(mclapply(1:nrow(minDf),function(i){
    allJaccard <- sapply(1:nrow(maxDf),function(j) CalcJaccard(minDf$Genes[[i]],maxDf$Genes[[j]]))
    jaccardIndex <- max(allJaccard)
    return(jaccardIndex)
  },mc.cores = 55))
  return(jaccardIndex)
}
#min max jaccard index of each pairs of methods
ComparisonMatrix <- matrix(NA,nrow = nrow(allPairwiseComparison),ncol=nrow(allPairwiseComparison))
for(i in 1:nrow(allPairwiseComparison)){
  for(j in (i+1):nrow(allPairwiseComparison)){
    if(j > ncol(ComparisonMatrix)){
      break
    }
    print(c(i,j))
    df1 <- allPASCALResults %>% dplyr::filter(Biotype == allPairwiseComparison$Biotype[i] & 
                           AggregationType == allPairwiseComparison$AggregationType[i] & 
                           Method == allPairwiseComparison$Method[i])
    df2 <- allPASCALResults %>% dplyr::filter(Biotype == allPairwiseComparison$Biotype[j] & 
                           AggregationType == allPairwiseComparison$AggregationType[j] & 
                           Method == allPairwiseComparison$Method[j])
    ComparisonMatrix[i,j] <- mean(FindJaccardPairs(df1,df2))
  }
}
save(ComparisonMatrix,file='../../Count_Data/JaccardAcrossMethods.rda')

#Plot comparision matrix as heat matrix.
load('../../Count_Data/JaccardAcrossMethods.rda')
rownames(ComparisonMatrix) <- c('CPG','APG','CPT','CWG','AWG','CWT')
colnames(ComparisonMatrix) <- c('CPG','APG','CPT','CWG','AWG','CWT')
diag(ComparisonMatrix) <- 1
ComparisonMatrix[lower.tri(ComparisonMatrix)] <- t(ComparisonMatrix)[lower.tri(ComparisonMatrix)]
library(corrplot)
corrplot(ComparisonMatrix, method = "circle",p.mat = ComparisonMatrix, insig = "p-value", sig.level = -1,
         is.corr = F,type = 'lower',tl.srt = 45,cl.pos = 'r',diag = T,
         tl.col = sapply(colnames(ComparisonMatrix),function(x) ifelse(substr(start = 2,stop = 2,x=x)=='P','blue','red')))
mtext(side = 4,text = 'Jaccard Index')
