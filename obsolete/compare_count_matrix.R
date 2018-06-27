##########################################################################################################
#Compare any two count matrix, created from different methods.
#Generates a correlation plot. 
##########################################################################################################

########################################Import Published Counts ##########################################
ImportPublishedCounts <- function(){
  source('sample_mapping.R')
  #Get GSM to sample name mapping
  mapping <- GetSampleMapping()
  
  runTable <- read.table(file = '/local/data/public/zmx21/zmx21_private/GSK/Galatro/SraRunTable.txt',header = T,sep = '\t')
  allGSM <- as.character(runTable$Sample_Name)
  
  publishedCounts <- data.matrix(read.table(file='/local/data/public/zmx21/zmx21_private/GSK/Galatro/GSE99074_HumanMicrogliaBrainCounts.txt'))
  #Keep only wanted GSM
  publishedCountsGSM <- sapply(colnames(publishedCounts),function(x) unique(as.character(mapping$GSM[mapping$title==x])))
  publishedCounts <- publishedCounts[,publishedCountsGSM%in%allGSM]
  colnames(publishedCounts) <- publishedCountsGSM[publishedCountsGSM%in%allGSM]
  
  #Get gene lengths
  library(biomaRt)
  biomartDb <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  df <- getBM(attributes=c("ensembl_gene_id","start_position","end_position","external_gene_name"),mart=biomartDb)
  df$gene_length <- abs(df$start_position - df$end_position)
  library(dplyr)
  allGenes <- data.frame(ensembl_gene_id = rownames(publishedCounts),stringsAsFactors = F)
  geneLengths <- dplyr::left_join(allGenes,df)$gene_length/1000

  #Remove genes which are not in db
  scaleFactor <- rowSums(t(publishedCounts))
  publishedCounts <- publishedCounts[!is.na(geneLengths),]
  allGenes <- data.frame(ensembl_gene_id = rownames(publishedCounts),stringsAsFactors = F)
  geneLengths <- dplyr::left_join(allGenes,df)$gene_length

  #Calculate RPKM
  publishedRPKM <- edgeR::rpkm(publishedCounts,geneLengths,normalized.lib.sizes	= T,log=F)
  publishedTPM <- do.call(cbind, lapply(1:ncol(publishedCounts), function(i) {
    rate = log(publishedCounts[,i]) - log(geneLengths)
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))


  colnames(publishedTPM) <- colnames(publishedCounts)
  rownames(publishedTPM) <- allGenes$ensembl_gene_id
  colnames(publishedRPKM) <- colnames(publishedCounts)
  rownames(publishedTPM) <- allGenes$ensembl_gene_id

  publishedTPM <- publishedTPM[!apply(publishedTPM,1,function(x) all(x==0)),]
  publishedCounts <- publishedCounts[!apply(publishedCounts,1,function(x) all(x==0)),]
  publishedRPKM <- publishedRPKM[!apply(publishedRPKM,1,function(x) all(x==0)),]

  return(list(publishedRPKM=publishedRPKM,publishedTPM = publishedTPM,publishedCounts = publishedCounts))
}
ImportPublishedNormalized<- function(){
  source('sample_mapping.R')
  #Get GSM to sample name mapping
  mapping <- GetSampleMapping()
  
  runTable <- read.table(file = '/local/data/public/zmx21/zmx21_private/GSK/Galatro/SraRunTable.txt',header = T,sep = '\t')
  allGSM <- as.character(runTable$Sample_Name)
  
  publishedCounts <- data.matrix(read.table(file='/local/data/public/zmx21/zmx21_private/GSK/Galatro/GSE99074_HumanMicrogliaBrainVoomNormalization.txt'))
  #Keep only wanted GSM
  publishedCountsGSM <- sapply(colnames(publishedCounts),function(x) unique(as.character(mapping$GSM[mapping$title==x])))
  publishedCounts <- publishedCounts[,publishedCountsGSM%in%allGSM]
  colnames(publishedCounts) <- publishedCountsGSM[publishedCountsGSM%in%allGSM]
  
  return(publishedCounts = publishedCounts)
}

########################################Compare Count Matrix ##########################################
source('import_salmon.R')
#Get pairs of genes from the two methods.
GetGenePairs <- function(method1,method2,logScale=F){
  #Choose genes which are in both count matrix
  geneNames <- intersect(rownames(method1),rownames(method2))
  sampleNames <- intersect(colnames(method1),colnames(method2))
  method1 <- method1[geneNames,]
  method2 <- method2[geneNames,]
  
  #Rearrage dataframe, so gene pairs of each sample and their TPM is on the same row.  
  if(logScale){
    expressionPairs <- lapply(sampleNames,function(x) data.frame(method1=log(method1[,x] + 1),method2=log(method2[,x] + 1)))
  }else{
    expressionPairs <- lapply(sampleNames,function(x) data.frame(method1=method1[,x],method2=method2[,x]))
  }
  #Merge gene pairs from all the samples. 
  allPairs <- do.call('rbind',expressionPairs)
  return(allPairs)
}

#Visualize similarity of alignment methods.
PlotMethodCor <- function(allPairs,xlab,ylab,title){
  maxLim <- max(c(allPairs$method1,allPairs$method2))
  ggplot(allPairs, aes(x=method1, y=method2) ) +
    geom_bin2d(aes(fill=log(..count..)),bins = 500)+
    theme_bw() + labs(x = xlab,y=ylab) + ggtitle(title) + xlim(0,maxLim) + ylim(0,maxLim)
}
#Visulization
tx2gene <- ImportTx2gene()
SalmonTPM <- ImportSalmonCounts('/local/data/public/zmx21/zmx21_private/GSK/Galatro/Salmon_aligned/Salmon_aligned_merged',tx2gene)

Salmon_FMD_TPM <- ImportSalmonCounts('/local/data/public/zmx21/zmx21_private/GSK/Galatro/Salmon_aligned/Salmon_aligned_fmd',tx2gene)
Salmon_Quasi_TPM <- ImportSalmonCounts('/local/data/public/zmx21/zmx21_private/GSK/Galatro/Salmon_aligned/Salmon_aligned_k19',tx2gene)
FMD_Against_Quasi <- GetGenePairs(Salmon_FMD_TPM$geneLevel$abundance,Salmon_Quasi_TPM$geneLevel$abundance)
PlotMethodCor(FMD_Against_Quasi,'FMD log(TPM + 1)','Quasi log(TPM + 1)','Salmon - FMD vs Quasi')

SalmonAgainstSTAR <- GetGenePairs(SalmonTPMGeneLevel,publishedTPM,logScale=T)
PlotMethodCor(SalmonAgainstSTAR,'Salmon log(TPM+1)','Published STAR log(TPM+1)','Salmon TPM vs Published Converted TPM')

