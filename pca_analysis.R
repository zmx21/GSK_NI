source('sample_mapping.R')
source('import_salmon.R')
mapping <- GetSampleMapping()
tx2gene <- ImportTx2gene()

SalmonUnmergedTPM <- ImportSalmonCounts('/local/data/public/zmx21/zmx21_private/GSK/Galatro/Salmon_aligned/Salmon_aligned_k19',tx2gene)
SalmonTPM <- ImportSalmonCounts('/local/data/public/zmx21/zmx21_private/GSK/Galatro/Salmon_aligned/Salmon_aligned_merged',tx2gene)

# expDfSalmon <- as.data.frame(t(SalmonUnmergedTPM$geneLevel$abundance),row.names = NULL)
# expDfSalmon$runID <- as.factor(colnames(SalmonUnmergedTPM$geneLevel$abundance))
# expDfSalmon$GSM <- as.factor(sapply(colnames(SalmonUnmergedTPM$geneLevel$abundance),function(x) mapping$GSM [mapping$SRR==x]))
# #Check PCA of multiple runs
# multRunDf <- expDfSalmon[(expDfSalmon$GSM %in% names(which(table(expDfSalmon$GSM) > 1))),]
# rownames(multRunDf) <- sapply(multRunDf$runID,function(x) paste0(unlist(strsplit(as.character(x),''))[9:10],collapse = ''))
# multRunPCA <- prcomp(multRunDf[,which(!colnames(expDfSalmon) %in% c('GSM','runID'))])
# autoplot(multRunPCA, data = multRunDf, colour = 'GSM',size=2,shape = FALSE, label.size = 3)

CollectMetadata <- function(inputMatrix){
  runTable <- read.table(file = '/local/data/public/zmx21/zmx21_private/GSK/Galatro/SraRunTable.txt',header = T,sep = '\t')
  allGSM <- colnames(inputMatrix$abundance)
  readLength <- sapply(allGSM,function(x) unique(subset(runTable,Sample_Name == x)$AvgSpotLen))
  gender <- sapply(allGSM,function(x) unique(subset(runTable,Sample_Name == x)$gender))
  age <- sapply(allGSM,function(x) unique(subset(runTable,Sample_Name == x)$age))
  numRuns <- sapply(allGSM,function(x) nrow(subset(runTable,Sample_Name == x)))
  
  df <- as.data.frame(t(inputMatrix$abundance),row.names = NULL)
  df$GSM <- allGSM; df$readLength <- as.factor(readLength)
  df$gender <- as.factor(gender); df$age <- age; df$numRuns <- as.factor(numRuns)
  return(df)
}
SalmonGeneLevelDf <- CollectMetadata(SalmonTPM$geneLevel)
# SalmonTranscriptLevel <- CollectMetadata(SalmonTPM$transcriptLevel)
SalmonGeneLevePCA <- prcomp(SalmonGeneLevelDf[,which(!colnames(SalmonGeneLevelDf) %in% c('GSM','readLength','gender','age','numRuns'))])
autoplot(SalmonGeneLevePCA, data = SalmonGeneLevelDf, colour = 'readLength',size=3) + ggtitle('PCA of read length')
autoplot(SalmonGeneLevePCA, data = SalmonGeneLevelDf, colour = 'numRuns',size=3) + ggtitle('PCA of number of runs')
autoplot(SalmonGeneLevePCA, data = SalmonGeneLevelDf, colour = 'gender',size=3)+ ggtitle('PCA of gender')
autoplot(SalmonGeneLevePCA, data = SalmonGeneLevelDf, colour = 'age',size=3)+ ggtitle('PCA of age')

