######################################################################################################
#Constructs correlation matrix  
#input should be count matrix, 
#output is Combat and PC corrected count matrix
######################################################################################################

########################################Construct Cor Matrix##########################################
#Load Batch corrected matrices
load('../Count_Data/Batch_Corrected/SalmonTPM_Gene_Combat.rda')
load('../Count_Data/Batch_Corrected/SalmonTPM_Gene_PCAdj.rda')
load('../Count_Data/Batch_Corrected/SalmonTPM_Transcript_Combat.rda')
load('../Count_Data/Batch_Corrected/SalmonTPM_Transcript_PCAdj.rda')

#Calculate and save correlation matrix for both Gene level and transcript level
SalmonCor_Gene_Combat <- cor(t(SalmonTPM_Gene_Combat))
SalmonCor_Transcript_Combat <- cor(t(SalmonTPM_Transcript_Combat))
SalmonCor_Gene_PCAdj1 <- cor(t(SalmonTPM_Gene_PCAdj$OnePC))
SalmonCor_Gene_PCAdj2 <- cor(t(SalmonTPM_Gene_PCAdj$TwoPC))
SalmonCor_Transcript_PCAdj1 <- cor(t(SalmonTPM_Transcript_PCAdj$OnePC))
SalmonCor_Transcript_PCAdj2 <- cor(t(SalmonTPM_Transcript_PCAdj$TwoPC))

save(SalmonCor_Gene_Combat,SalmonCor_Transcript_Combat,file = '../Count_Data/Correlation_Matrices/Combat_Cor.rda')
save(SalmonCor_Gene_PCAdj1,SalmonCor_Gene_PCAdj2,SalmonCor_Transcript_PCAdj1,SalmonCor_Transcript_PCAdj2,file = '../Count_Data/Correlation_Matrices/PCAdj_Cor.rda')


###################################Compared with gene friends##########################################
#Import gene friends. Data is stored such that each txt files (named according to the ENSG ID), 
#Containts list of it's co-expressed partner genes. 
load('../Count_Data/Correlation_Matrices/Combat_Cor.rda')
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
#Import Only Common Genes and saveS
geneFriendsCorMatrix <- ImportGeneFriends(rownames(SalmonCor_Gene_Combat))
save(geneFriendsCorMatrix,file='../Count_Data/Correlation_Matrices/geneFriendsCorMatrix.rda')

#Import all genes and save
geneFriendsAllCorMatrix <- ImportGeneFriends(importAll = T)
save(geneFriendsAllCorMatrix,file='../Count_Data/Correlation_Matrices/geneFriendsAllCorMatrix.rda')

load('../Count_Data/Correlation_Matrices/geneFriendsCorMatrix.rda')
load('../Count_Data/Correlation_Matrices/PCAdj_Cor.rda')

#Extract upper diagonal of correlation matricies
upper_GeneFriends <- geneFriendsCorMatrix[upper.tri(geneFriendsCorMatrix)]
upper_Combat <- SalmonCor_Gene_Combat[upper.tri(SalmonCor_Gene_Combat)]
upper_PCAdj1 <- SalmonCor_Gene_PCAdj1[upper.tri(SalmonCor_Gene_PCAdj1)]
upper_PCAdj2 <- SalmonCor_Gene_PCAdj2[upper.tri(SalmonCor_Gene_PCAdj2)]

#Remove NA elements
upper_Combat <- upper_Combat[!is.na(upper_GeneFriends)]
upper_GeneFriends <- upper_GeneFriends[!is.na(upper_GeneFriends)]
upper_Random <- upper_Combat[sample(1:length(upper_Combat),replace = F,size=length(upper_Combat))]
upper_PCAdj1 <- upper_PCAdj1[!is.na(upper_GeneFriends)]
upper_PCAdj2 <- upper_PCAdj2[!is.na(upper_GeneFriends)]

#Sample randomly, number of pairs to reduce compuational time
randSample <- sample(1:length(upper_Combat),size=1000000,replace=F)
compDf <- data.frame(GeneFriends = upper_GeneFriends[randSample],
                     Combat = upper_Combat[randSample],
                     Random = upper_Random[randSample],
                     PCAdj1 = upper_PCAdj1[randSample],
                     PCAdj2 = upper_PCAdj2[randSample])
#Plot correlation, where x is Gene Friends and Y is batch corrected.
p1 <- ggplot(compDf, aes(x=GeneFriends, y=Combat) ) +
  geom_bin2d(aes(fill=log(..count..)),bins = 500)+
  theme_bw() + labs(x = 'Gene Coexpression Database Cor Value',y='Combat Cor Value') + 
  ggtitle('Combat Cor Matrix \n vs. Existing Gene Co-exp Cor Matrix') + xlim(-1,1) + ylim(-1,1)

p2 <- ggplot(compDf, aes(x=GeneFriends, y=Random) ) +
  geom_bin2d(aes(fill=log(..count..)),bins = 500)+
  theme_bw() + labs(x = 'Gene Coexpression Database Cor Value',y='Random Cor Value') + 
  ggtitle('Randomly Permuted Cor Matrix \n vs. Existing Gene Co-exp Cor Matrix') + xlim(-1,1) + ylim(-1,1)
ggarrange(p1,p2,ncol=2)

p3 <- ggplot(compDf, aes(x=GeneFriends, y=PCAdj1) ) +
  geom_bin2d(aes(fill=log(..count..)),bins = 500)+
  theme_bw() + labs(x = 'Gene Coexpression Database Cor Value',y='PCAdj 1 Cor Value') + 
  ggtitle('PCAdj 1 Cor Matrix \n vs. Existing Gene Co-exp Cor Matrix') + xlim(-1,1) + ylim(-1,1)
p4 <- ggplot(compDf, aes(x=GeneFriends, y=PCAdj2) ) +
  geom_bin2d(aes(fill=log(..count..)),bins = 500)+
  theme_bw() + labs(x = 'Gene Coexpression Database Cor Value',y='PCAdj 1 and 2 Cor Value') + 
  ggtitle('PCAdj 1 and 2 Cor Matrix \n vs. Existing Gene Co-exp Cor Matrix') + xlim(-1,1) + ylim(-1,1)
ggarrange(p3,p4,ncol=2)

