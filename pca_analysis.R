##########################################################################################################
#PCA based on gene expressions. Imports metadata, which could be visualzed as colors on PCA plot.
##########################################################################################################
#Impot required Data
source('sample_mapping.R')
source('import_salmon.R')
mapping <- GetSampleMapping()
tx2gene <- ImportTx2gene()
SalmonTPM <- ImportSalmonCounts('/local/data/public/zmx21/zmx21_private/GSK/Galatro/Salmon_aligned/Salmon_aligned_merged',tx2gene)
load('../Count_Data/SalmonTPM_Gene_Microglia_Filt.rda')

#Check PCA of individual runs
# SalmonUnmergedTPM <- ImportSalmonCounts('/local/data/public/zmx21/zmx21_private/GSK/Galatro/Salmon_aligned/Salmon_aligned_k19',tx2gene)
# expDfSalmon <- as.data.frame(t(SalmonUnmergedTPM$geneLevel$abundance),row.names = NULL)
# expDfSalmon$runID <- as.factor(colnames(SalmonUnmergedTPM$geneLevel$abundance))
# expDfSalmon$GSM <- as.factor(sapply(colnames(SalmonUnmergedTPM$geneLevel$abundance),function(x) mapping$GSM [mapping$SRR==x]))
# multRunDf <- expDfSalmon[(expDfSalmon$GSM %in% names(which(table(expDfSalmon$GSM) > 1))),]
# rownames(multRunDf) <- sapply(multRunDf$runID,function(x) paste0(unlist(strsplit(as.character(x),''))[9:10],collapse = ''))
# multRunPCA <- prcomp(multRunDf[,which(!colnames(expDfSalmon) %in% c('GSM','runID'))])
# autoplot(multRunPCA, data = multRunDf, colour = 'GSM',size=2,shape = FALSE, label.size = 3)

#Collections metadata, from a count matrix (where colnames are sample names)
#also need the tables, which is a list of tables with run info, alignment info, and read distribution info. 
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
#Get metadata (downloaded from SRA), alignment info, and read distribtuion
runTable <- read.table(file = '/local/data/public/zmx21/zmx21_private/GSK/Galatro/SraRunTable.txt',header = T,sep = '\t')
alignmentTable <- read.table(file='/local/data/public/zmx21/zmx21_private/GSK/Galatro/Salmon_aligned_merged/multiqc_Salmon_merged/multiqc_general_stats.txt',header = T)
readDist <- read.table(file='/local/data/public/zmx21/zmx21_private/GSK/Galatro/STAR_aligned_merged/multiqc_data/multiqc_rseqc_read_distribution.txt',header=T)
tables <- list(runTable=runTable,alignmentTable=alignmentTable,readDist=readDist)

CalcPCA <- function(countMatrix,tables){
  Df <- CollectMetadata(countMatrix,tables)
  
  #Extact Median and SD of exp of each gene/transcript. Store in Df
  MedianExpByGSM <- sapply(Df$GSM,function(x) median(countMatrix[,as.character(x)]))
  SDByGSM <- sapply(Df$GSM,function(x) sd(countMatrix[,as.character(x)]))
  Df <- cbind(Df,data.frame(median_exp=MedianExpByGSM,sd=SDByGSM))
  rownames(Df) <- sapply(Df$GSM,function(x) paste0(unlist(strsplit(as.character(x),''))[9:10],collapse = ''))
  
  #Calculates PCA, based on gene/transcript expression values of different samples. 
  PCA <- prcomp(Df[,which(!colnames(Df) %in% c('GSM','readLength','gender','age','numRuns','numReads','median_exp','sd','exonReads'))])
  return(list(PCA = PCA,Df=Df))
}

#Calc PCA and plot results.
results <- CalcPCA(SalmonTPM_Gene_Filt,tables)
SalmonGeneLevelDf <- results$Df
SalmonGeneLevelPCA <- results$PCA

library(ggplot2)
library(ggfortify)
autoplot(results$PCA, data = results$Df, colour = 'readLength',size=4,shape=F) + ggtitle('PCA of Read Length')

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

