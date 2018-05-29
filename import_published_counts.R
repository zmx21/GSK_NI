ImportPublishedCounts<- function(){
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
  
  # #Get gene lengths
  # library(biomaRt)
  # biomartDb <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  # df <- getBM(attributes=c("ensembl_gene_id","start_position","end_position","external_gene_name"),mart=biomartDb)
  # df$gene_length <- abs(df$start_position - df$end_position)
  # library(dplyr)
  # allGenes <- data.frame(ensembl_gene_id = rownames(publishedCounts),stringsAsFactors = F)
  # geneLengths <- dplyr::left_join(allGenes,df)$gene_length/1000
  # 
  # #Remove genes which are not in db
  # scaleFactor <- rowSums(t(publishedCounts))
  # publishedCounts <- publishedCounts[!is.na(geneLengths),]
  # allGenes <- data.frame(ensembl_gene_id = rownames(publishedCounts),stringsAsFactors = F)
  # geneLengths <- dplyr::left_join(allGenes,df)$gene_length
  # 
  # # publishedRPKM <- do.call(cbind, lapply(1:ncol(publishedCounts), function(i) {
  # #   rpk <- publishedCounts[,i] / geneLengths
  # #   rpk / (scaleFactor[i] / 1e6)
  # # }))
  # publishedRPKM <- edgeR::rpkm(publishedCounts,geneLengths,normalized.lib.sizes	= T,log=F)
  # publishedTPM <- do.call(cbind, lapply(1:ncol(publishedCounts), function(i) {
  #   rate = log(publishedCounts[,i]) - log(geneLengths)
  #   denom = log(sum(exp(rate)))
  #   exp(rate - denom + log(1e6))
  # }))
  # 
  # 
  # colnames(publishedTPM) <- colnames(publishedCounts)
  # rownames(publishedTPM) <- allGenes$ensembl_gene_id
  # colnames(publishedRPKM) <- colnames(publishedCounts)
  # rownames(publishedTPM) <- allGenes$ensembl_gene_id
  # 
  # publishedTPM <- publishedTPM[!apply(publishedTPM,1,function(x) all(x==0)),]
  # publishedCounts <- publishedCounts[!apply(publishedCounts,1,function(x) all(x==0)),]
  # publishedRPKM <- publishedRPKM[!apply(publishedRPKM,1,function(x) all(x==0)),]
  # 
  return(publishedCounts)
  # return(list(publishedRPKM=publishedRPKM,publishedTPM = publishedTPM,publishedCounts = publishedCounts))
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