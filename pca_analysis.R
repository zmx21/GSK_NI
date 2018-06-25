##########################################################################################################
#PCA based on gene expressions. Imports metadata, which could be visualzed as colors on PCA plot.
##########################################################################################################
#Impot required Data
source('sample_mapping.R')
source('import_salmon.R')
# mapping <- GetSampleMapping()
# tx2gene <- ImportTx2gene()
# SalmonTPM <- ImportSalmonCounts('/local/data/public/zmx21/zmx21_private/GSK/Galatro/Salmon_aligned/Salmon_aligned_merged',tx2gene)
# SalmonTPM <- ImportSalmonCounts('/local/data/public/zmx21/zmx21_private/GSK/Gosselin/salmon/',tx2gene)

# load('../Count_Data/SalmonTPM_Gene_Microglia_Filt.rda')

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
CollectMetadata <- function(inputMatrix,full=F){
  library(dplyr)
  #Get metadata (downloaded from SRA), alignment info, and read distribtuion
  runTable_Galatro <- read.table(file = '/local/data/public/zmx21/zmx21_private/GSK/Galatro/SraRunTable.txt',header = T,sep = '\t',stringsAsFactors = F) %>% 
    dplyr::select(Sample_Name,AvgSpotLen,gender,age) %>%
    {.[!duplicated(.),]}  #Each sample has two runs, so remove duplicated rows
  runTable_Galatro_Brain <- read.table(file = '/local/data/public/zmx21/zmx21_private/GSK/Galatro_Brain//SraRunTable.txt',header = T,sep = '\t',stringsAsFactors = F) %>% 
    dplyr::select(Sample_Name,AvgSpotLen,gender,age) %>%
    {.[!duplicated(.),]}  #Each sample has two runs, so remove duplicated rows
  
  runTable_Galatro$BulkBrain <- rep(0,nrow(runTable_Galatro))
  runTable_Galatro_Brain$BulkBrain <- rep(1,nrow(runTable_Galatro_Brain)) 
  runTable_Galatro <- rbind(runTable_Galatro,runTable_Galatro_Brain)
  
  alignmentTable_Galatro <- rbind(read.table(file='/local/data/public/zmx21/zmx21_private/GSK/Galatro/Salmon_aligned_merged/multiqc_Salmon_merged//multiqc_general_stats.txt',header = T,stringsAsFactors = F),
                                  read.table(file='/local/data/public/zmx21/zmx21_private/GSK/Galatro_Brain/Salmon_aligned_whole_brain/multiqc_data/multiqc_general_stats.txt',header = T,stringsAsFactors = F))
  readDist_Galatro <- rbind(read.table(file='/local/data/public/zmx21/zmx21_private/GSK/Galatro/STAR_aligned_merged/multiqc_data/multiqc_rseqc_read_distribution.txt',header=T,stringsAsFactors = F),
                            read.table(file='/local/data/public/zmx21/zmx21_private/GSK/Galatro_Brain/STAR_aligned_whole_brain/multiqc_data/multiqc_rseqc_read_distribution.txt',header=T,stringsAsFactors = F)) 
  readDist_Galatro$Sample <- sapply(readDist_Galatro$Sample,function(x) unlist(strsplit(x = x,split = '[.]'))[1])
  
  runTable_Gosselin <- read.table(file = '/local/data/public/zmx21/zmx21_private/GSK/Gosselin/SraRunTable_Parsed.txt',header = T,sep = '\t',stringsAsFactors = F) %>% 
    dplyr::select(Sample_Name,AvgSpotLen,gender,age) %>% dplyr::mutate(BulkBrain = 0)
  alignmentTable_Gosselin <- read.table(file='/local/data/public/zmx21/zmx21_private/GSK/Gosselin/salmon/multiqc_data/multiqc_general_stats.txt',header = T,stringsAsFactors = F)
  readDist_Gosselin <- read.table(file='/local/data/public/zmx21/zmx21_private/GSK/Gosselin/QC_Reports/star_multiqc_data/multiqc_rseqc_read_distribution.txt',header=T,stringsAsFactors = F) %>%
    dplyr::mutate(other_intergenic_tag_count=NA,other_intergenic_tag_pct=NA)

  alignmentTable_Olah <- read.table(file='/local/data/public/zmx21/zmx21_private/GSK/Olah/salmon/multiqc_data/multiqc_general_stats.txt',header = T,stringsAsFactors = F)
  readDist_Olah <- read.table(file='/local/data/public/zmx21/zmx21_private/GSK/Olah/QC_Reports/star_multiqc_data/multiqc_rseqc_read_distribution.txt',header=T,stringsAsFactors = F) %>%
    dplyr::mutate(other_intergenic_tag_count=NA,other_intergenic_tag_pct=NA)
  
  
  tables <- list(runTable=rbind(runTable_Galatro,runTable_Gosselin),
                 alignmentTable=rbind(alignmentTable_Galatro,alignmentTable_Gosselin,alignmentTable_Olah),
                 readDist=rbind(readDist_Galatro,readDist_Gosselin,readDist_Olah))
  
  
  allSamples <- colnames(inputMatrix)
  if(full){
    df <- as.data.frame(t(inputMatrix),row.names = NULL)
    df$Sample_Name <- allSamples
  }else{
    df <- data.frame(Sample_Name = allSamples)
  }

  dataset <- sapply(allSamples,function(x) ifelse(grepl('GSM',x),'Galatro',ifelse(grepl('SRR',x),'Gosselin','Olah')))
  df$dataset <- dataset
  #Join all metadata tables, keep information needed. 
  dfMetadata <- dplyr::left_join(data.frame(Sample_Name = df$Sample_Name,stringsAsFactors = F),
                                 dplyr::select(tables$runTable,readLength=AvgSpotLen,gender=gender,age=age,Sample_Name,BulkBrain),by=c('Sample_Name' = 'Sample_Name')) %>% 
    dplyr::left_join(dplyr::select(tables$alignmentTable,Sample_Name = Sample,
                                   numReads=Salmon_num_mapped,mappingRate=Salmon_percent_mapped),
                     by = c('Sample_Name'='Sample_Name')) %>% 
    dplyr::left_join(dplyr::select(tables$readDist,Sample_Name=Sample,exonTags = cds_exons_tag_count,
                                   intronTags = introns_tag_count,totalTags = total_tags,intergenicTags=other_intergenic_tag_count),by=c('Sample_Name'='Sample_Name'))
  df <- cbind(df,dfMetadata)
  #Add expType and instrument for Gosselin Samples
  GosselinAddtlData <- read.table(file = '/local/data/public/zmx21/zmx21_private/GSK/Gosselin/SraRunTable_Parsed.txt',header = T,sep = '\t',stringsAsFactors = F) %>%
    dplyr::select(Sample_Name,Instrument,Library_Name) %>% 
    dplyr::mutate(expType = ifelse(grepl('ExVivo',Library_Name),'ExVivo','InVitro')) %>% 
    dplyr::select(-Library_Name)
  df <- cbind(df,left_join(data.frame(Sample_Name=df$Sample_Name,stringsAsFactors = F),GosselinAddtlData,by=c('Sample_Name' = 'Sample_Name')))
  rownames(df) <- sapply(df$Sample_Name,function(x) paste0(ifelse(substr(x,1,3)=='GSM','L','S'),substr(x,nchar(x)-1,nchar(x))))
  return(df)
}

CalcPCA <- function(countMatrix){
  Df <- CollectMetadata(countMatrix,full=T)
    
  #Calculates PCA, based on gene/transcript expression values of different samples. 
  allGenes <- rownames(countMatrix)
  PCA <- prcomp(Df[,which(colnames(Df) %in% allGenes)])
  return(list(PCA = PCA,Df=Df))
}

#Calc PCA and plot results.
# results <- CalcPCA(SalmonTPM_Gene)
# SalmonGeneLevelDf <- results$Df
# SalmonGeneLevelPCA <- results$PCA

library(ggplot2)
library(ggfortify)
# autoplot(results$PCA, data = results$Df, colour = 'readLength',size=4,shape=F) + ggtitle('PCA of Read Length')
#   scale_colour_gradientn(colours = rainbow(7))
# autoplot(SalmonGeneLevelPCA, data = SalmonGeneLevelDf, colour = 'sd',size=4,shape=F) + ggtitle('PCA of SD') + 
#   scale_colour_gradientn(colours = rainbow(7))
# autoplot(SalmonGeneLevelPCA, data = SalmonGeneLevelDf, colour = 'mappingRate',size=4,shape=F) + ggtitle('PCA of Mapping Rate') + 
#   scale_colour_gradientn(colours = rainbow(7))
# autoplot(SalmonGeneLevelPCA, data = SalmonGeneLevelDf, colour = 'exonReads',size=4,shape=F) + ggtitle('PCA of Mapping Rate') +
#   scale_colour_gradientn(colours = rainbow(7))

