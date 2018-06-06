source('sample_mapping.R')
# source('import_salmon.R')
# mapping <- GetSampleMapping()
# tx2gene <- ImportTx2gene()

# SalmonUnmergedTPM <- ImportSalmonCounts('/local/data/public/zmx21/zmx21_private/GSK/Galatro/Salmon_aligned/Salmon_aligned_k19',tx2gene)
# SalmonTPM <- ImportSalmonCounts('/local/data/public/zmx21/zmx21_private/GSK/Galatro/Salmon_aligned/Salmon_aligned_merged',tx2gene)
# 
# expDfSalmon <- as.data.frame(t(SalmonUnmergedTPM$geneLevel$abundance),row.names = NULL)
# expDfSalmon$runID <- as.factor(colnames(SalmonUnmergedTPM$geneLevel$abundance))
# expDfSalmon$GSM <- as.factor(sapply(colnames(SalmonUnmergedTPM$geneLevel$abundance),function(x) mapping$GSM [mapping$SRR==x]))
# #Check PCA of multiple runs
# multRunDf <- expDfSalmon[(expDfSalmon$GSM %in% names(which(table(expDfSalmon$GSM) > 1))),]
# rownames(multRunDf) <- sapply(multRunDf$runID,function(x) paste0(unlist(strsplit(as.character(x),''))[9:10],collapse = ''))
# multRunPCA <- prcomp(multRunDf[,which(!colnames(expDfSalmon) %in% c('GSM','runID'))])
# autoplot(multRunPCA, data = multRunDf, colour = 'GSM',size=2,shape = FALSE, label.size = 3)

# load('../Count_Data/SalmonTPM_Gene_Microglia_Filt.rda')

CollectMetadata <- function(inputMatrix,tables){
  allGSM <- colnames(inputMatrix)
  readLength <- sapply(allGSM,function(x) unique(subset(tables$runTable,Sample_Name == x)$AvgSpotLen))
  gender <- sapply(allGSM,function(x) unique(subset(tables$runTable,Sample_Name == x)$gender))
  age <- sapply(allGSM,function(x) unique(subset(tables$runTable,Sample_Name == x)$age))
  numRuns <- sapply(allGSM,function(x) nrow(subset(tables$runTable,Sample_Name == x)))
  numReads <- sapply(allGSM,function(x) subset(tables$alignmentTable,Sample==x)$Salmon_num_mapped)
  mappingRate <- sapply(allGSM,function(x) subset(tables$alignmentTable,Sample==x)$Salmon_percent_mapped)
  exonReads <- sapply(allGSM,function(x) subset(tables$readDist,Sample==paste0(x,'.read.dist'))$cds_exons_tag_count)
  
  df <- as.data.frame(t(inputMatrix),row.names = NULL)
  df$GSM <- allGSM; df$readLength <- as.factor(readLength)
  df$gender <- as.factor(gender); df$age <- age; df$numRuns <- as.factor(numRuns)
  df$numReads <- numReads; df$mappingRate <- mappingRate; df$exonReads <- exonReads
  return(df)
}
#Get metadata for Salmon TPM
runTable <- read.table(file = '/local/data/public/zmx21/zmx21_private/GSK/Galatro/SraRunTable.txt',header = T,sep = '\t')
alignmentTable <- read.table(file='/local/data/public/zmx21/zmx21_private/GSK/Galatro/Salmon_aligned_merged/multiqc_Salmon_merged/multiqc_general_stats.txt',header = T)
readDist <- read.table(file='/local/data/public/zmx21/zmx21_private/GSK/Galatro/STAR_aligned_merged/multiqc_data/multiqc_rseqc_read_distribution.txt',header=T)
tables <- list(runTable=runTable,alignmentTable=alignmentTable,readDist=readDist)

CalcPCA <- function(countMatrix,tables){
  Df <- CollectMetadata(countMatrix,tables)
  
  cv <- sapply(1:nrow(countMatrix),function(i) sd(countMatrix[i,])/mean(countMatrix[i,]))
  countMatrix <- countMatrix[order(cv,decreasing = F)[1:2],]
  
  MedianExpByGSM <- sapply(Df$GSM,function(x) median(countMatrix[,as.character(x)]))
  SDByGSM <- sapply(Df$GSM,function(x) sd(countMatrix[,as.character(x)]))
  Df <- cbind(Df,data.frame(median_exp=MedianExpByGSM,sd=SDByGSM))
  rownames(Df) <- sapply(Df$GSM,function(x) paste0(unlist(strsplit(as.character(x),''))[9:10],collapse = ''))
  
  PCA <- prcomp(Df[,which(!colnames(Df) %in% c('GSM','readLength','gender','age','numRuns','numReads','median_exp','sd','exonReads'))])
  return(list(PCA = PCA,Df=Df))
}

# badSamples <- c('GSM3081110','GSM3081130','GSM3081125')
# results <- lapply(SalmonTPM_Gene_Filt,function(x) CalcPCA(x,tables))
# SalmonGeneLevelDf <- lapply(results, function(x) x$Df)
# SalmonGeneLevelPCA <- lapply(results,function(x) x$PCA)

# library(ggplot2)
# library(ggfortify)
# autoplot(results$PCA, data = results$Df, colour = 'readLength',size=4,shape=F) + ggtitle('PCA of Read Length') 

# autoplot(SalmonGeneLevelPCA, data = SalmonGeneLevelDf, colour = 'readLength',size=4,shape=F) + ggtitle('PCA of Read Length') 
# autoplot(SalmonGeneLevelPCA, data = SalmonGeneLevelDf, colour = 'numReads',size=3) + ggtitle('PCA of number of reads')
# autoplot(SalmonGeneLevelPCA, data = SalmonGeneLevelDf, colour = 'numRuns',size=3) + ggtitle('PCA of number of runs')
# autoplot(SalmonGeneLevelPCA, data = SalmonGeneLevelDf, colour = 'gender',size=3)+ ggtitle('PCA of gender')
# autoplot(SalmonGeneLevelPCA, data = SalmonGeneLevelDf, colour = 'age',size=3)+ ggtitle('PCA of age')
# autoplot(SalmonGeneLevelPCA, data = SalmonGeneLevelDf, colour = 'median_exp',shape=F,size=4) + ggtitle('PCA of Median Exp') +
#   scale_colour_gradientn(colours = rainbow(7))
# autoplot(SalmonGeneLevelPCA, data = SalmonGeneLevelDf, colour = 'sd',size=4,shape=F) + ggtitle('PCA of SD') + 
#   scale_colour_gradientn(colours = rainbow(7))
# autoplot(SalmonGeneLevelPCA, data = SalmonGeneLevelDf, colour = 'mappingRate',size=4,shape=F) + ggtitle('PCA of Mapping Rate') + 
#   scale_colour_gradientn(colours = rainbow(7))
# autoplot(SalmonGeneLevelPCA, data = SalmonGeneLevelDf, colour = 'exonReads',size=4,shape=F) + ggtitle('PCA of Mapping Rate') +
#   scale_colour_gradientn(colours = rainbow(7))


#Get metadata for Published TPM
# PublishedGeneLevelDf <- CollectMetadata(publishedTPM)
# PublishedGeneLevelPCA <- prcomp(PublishedGeneLevelDf[,which(!colnames(PublishedGeneLevelDf) %in% c('GSM','readLength','gender','age','numRuns','numReads'))])
# autoplot(PublishedGeneLevelPCA, data = PublishedGeneLevelDf, colour = 'readLength',size=3) + ggtitle('PCA of read length')
# autoplot(PublishedGeneLevelPCA, data = PublishedGeneLevelDf, colour = 'numReads',size=3) + ggtitle('PCA of number of reads')
# autoplot(PublishedGeneLevelPCA, data = PublishedGeneLevelDf, colour = 'numRuns',size=3) + ggtitle('PCA of number of runs')
# autoplot(PublishedGeneLevelPCA, data = PublishedGeneLevelDf, colour = 'gender',size=3)+ ggtitle('PCA of gender')
# autoplot(PublishedGeneLevelPCA, data = PublishedGeneLevelDf, colour = 'age',size=3)+ ggtitle('PCA of age')

# p1 <- autoplot(SalmonGeneLevelPCA[[1]], data = SalmonGeneLevelDf[[1]], colour = 'readLength',size=4,shape=F) + ggtitle('TPM>5, CV>0.3') 
# p2 <- autoplot(SalmonGeneLevelPCA[[2]], data = SalmonGeneLevelDf[[2]], colour = 'readLength',size=4,shape=F) + ggtitle('TPM>5,CV>0.5') 
# p3 <- autoplot(SalmonGeneLevelPCA[[3]], data = SalmonGeneLevelDf[[3]], colour = 'readLength',size=4,shape=F) + ggtitle('TPM>10,CV>0.3') 
# p4 <- autoplot(SalmonGeneLevelPCA[[4]], data = SalmonGeneLevelDf[[4]], colour = 'readLength',size=4,shape=F) + ggtitle('TPM>10,CV>0.5') 
# grid.arrange(p1,p2,p3,p4)

