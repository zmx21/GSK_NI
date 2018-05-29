source('sample_mapping.R')
source('import_salmon.R')
mapping <- GetSampleMapping()
tx2gene <- ImportTx2gene()

# SalmonUnmergedTPM <- ImportSalmonCounts('/local/data/public/zmx21/zmx21_private/GSK/Galatro/Salmon_aligned/Salmon_aligned_k19',tx2gene)
SalmonTPM <- ImportSalmonCounts('/local/data/public/zmx21/zmx21_private/GSK/Galatro/Salmon_aligned/Salmon_aligned_merged',tx2gene)
# 
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
  alignmentTable <- read.table(file='/local/data/public/zmx21/zmx21_private/GSK/Galatro/Salmon_aligned/Salmon_aligned_merged/multiqc_Salmon_merged/multiqc_general_stats.txt',header = T)
  
  allGSM <- colnames(inputMatrix)
  readLength <- sapply(allGSM,function(x) unique(subset(runTable,Sample_Name == x)$AvgSpotLen))
  gender <- sapply(allGSM,function(x) unique(subset(runTable,Sample_Name == x)$gender))
  age <- sapply(allGSM,function(x) unique(subset(runTable,Sample_Name == x)$age))
  numRuns <- sapply(allGSM,function(x) nrow(subset(runTable,Sample_Name == x)))
  numReads <- sapply(allGSM,function(x) subset(alignmentTable,Sample==x)$Salmon_num_mapped)
  
  df <- as.data.frame(t(inputMatrix),row.names = NULL)
  df$GSM <- allGSM; df$readLength <- as.factor(readLength)
  df$gender <- as.factor(gender); df$age <- age; df$numRuns <- as.factor(numRuns)
  df$numReads <- numReads
  return(df)
}
#Get metadata for Salmon TPM
SalmonGeneLevelDf <- CollectMetadata(SalmonTPM$geneLevel$abundance)
# # SalmonTranscriptLevel <- CollectMetadata(SalmonTPM$transcriptLevel)
# SalmonGeneLevePCA <- prcomp(SalmonGeneLevelDf[,which(!colnames(SalmonGeneLevelDf) %in% c('GSM','readLength','gender','age','numRuns','numReads'))])
# autoplot(SalmonGeneLevePCA, data = SalmonGeneLevelDf, colour = 'readLength',size=3) + ggtitle('PCA of read length')
# autoplot(SalmonGeneLevePCA, data = SalmonGeneLevelDf, colour = 'numReads',size=3) + ggtitle('PCA of number of reads')
# autoplot(SalmonGeneLevePCA, data = SalmonGeneLevelDf, colour = 'numRuns',size=3) + ggtitle('PCA of number of runs')
# autoplot(SalmonGeneLevePCA, data = SalmonGeneLevelDf, colour = 'gender',size=3)+ ggtitle('PCA of gender')
# autoplot(SalmonGeneLevePCA, data = SalmonGeneLevelDf, colour = 'age',size=3)+ ggtitle('PCA of age')

#Get metadata for Published TPM
# PublishedGeneLevelDf <- CollectMetadata(publishedTPM)
# PublishedGeneLevelPCA <- prcomp(PublishedGeneLevelDf[,which(!colnames(PublishedGeneLevelDf) %in% c('GSM','readLength','gender','age','numRuns','numReads'))])
# autoplot(PublishedGeneLevelPCA, data = PublishedGeneLevelDf, colour = 'readLength',size=3) + ggtitle('PCA of read length')
# autoplot(PublishedGeneLevelPCA, data = PublishedGeneLevelDf, colour = 'numReads',size=3) + ggtitle('PCA of number of reads')
# autoplot(PublishedGeneLevelPCA, data = PublishedGeneLevelDf, colour = 'numRuns',size=3) + ggtitle('PCA of number of runs')
# autoplot(PublishedGeneLevelPCA, data = PublishedGeneLevelDf, colour = 'gender',size=3)+ ggtitle('PCA of gender')
# autoplot(PublishedGeneLevelPCA, data = PublishedGeneLevelDf, colour = 'age',size=3)+ ggtitle('PCA of age')

