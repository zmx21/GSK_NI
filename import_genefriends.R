# load('../Count_Data/Correlation_Matrices/Combat_Cor.rda')
library(dplyr)
library(pbmcapply)
ImportGeneFriends <- function(geneIDs=c(),importAll=F){
  ImportTxtFile <- function(path,gene,geneIDs){
    mergedDf <- tryCatch({
      df <- read.table(paste0(path,gene,'.txt'),stringsAsFactors = F)
      colnames(df) <- 'Cor'
      dplyr::left_join(data.frame(Sample=geneIDs,stringsAsFactors = F),
                       data.frame(Sample=rownames(df),Cor=df$Cor,stringsAsFactors = F),by="Sample")
    }, warning = function(w) {
      data.frame(Sample=geneIDs,Cor=rep(NA,length(geneIDs)),stringsAsFactors = F)
    }, error = function(e) {
      data.frame(Sample=geneIDs,Cor=rep(NA,length(geneIDs)),stringsAsFactors = F)
    }, finally = {
    })
    return(mergedDf)
  }
  path <- '/local/data/public/zmx21/zmx21_private/GSK/Genes10Perc/'
  if(importAll){
    geneIDs <- sapply(dir(path),function(x) unlist(strsplit(x,'\\.'))[1])
    allDf <- pbmclapply(1:length(geneIDs),function(i) ImportTxtFile(path,geneIDs[i],geneIDs)$Cor,mc.cores = 5)
  }else{
    allDf <- pbmclapply(1:length(geneIDs),function(i) ImportTxtFile(path,geneIDs[i],geneIDs)$Cor,mc.cores = 5)
  }
  corMatrix <- do.call(rbind,allDf)
  rownames(corMatrix) <- geneIDs
  colnames(corMatrix) <- geneIDs
  return(corMatrix)
}
# geneFriendsCorMatrix <- ImportGeneFriends(rownames(SalmonCor_Gene_Combat))
# save(geneFriendsCorMatrix,file='../Count_Data/Correlation_Matrices/geneFriendsCorMatrix.rda')

# geneFriendsAllCorMatrix <- ImportGeneFriends(importAll = T)
# save(geneFriendsAllCorMatrix,file='../Count_Data/Correlation_Matrices/geneFriendsAllCorMatrix.rda')