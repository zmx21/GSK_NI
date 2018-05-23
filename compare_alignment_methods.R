source('import_salmon.R')
GetGenePairs <- function(method1,method2){
  geneNames <- intersect(rownames(method1),rownames(method2))
  sampleNames <- intersect(colnames(method1),colnames(method2))
  method1 <- method1[geneNames,]
  method2 <- method2[geneNames,]
  
  expressionPairs <- lapply(sampleNames,function(x) data.frame(method1=log(method1[,x] + 1),method2=log(method2[,x] + 1)))
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
tx2gene <- ImportTx2gene()
SalmonTPM <- ImportSalmonCounts('../Salmon_aligned/Salmon_aligned_merged/',tx2gene)

Salmon_FMD_TPM <- ImportSalmonCounts('../Salmon_aligned/Salmon_aligned_fmd',tx2gene)
Salmon_Quasi_TPM <- ImportSalmonCounts('../Salmon_aligned/Salmon_aligned_k19',tx2gene)
FMD_Against_Quasi <- GetGenePairs(Salmon_FMD_TPM$geneLevel$abundance,Salmon_Quasi_TPM$geneLevel$abundance)
PlotMethodCor(FMD_Against_Quasi,'FMD log(TPM + 1)','Quasi log(TPM + 1)','Salmon - FMD vs Quasi')

source('sample_mapping.R')
mapping <- GetSampleMapping()
allGSM <- colnames(SalmonTPM$geneLevel$abundance)

publishedCounts <- read.table(file='/local/data/public/zmx21/zmx21_private/GSK/Galatro/GSE99074_HumanMicrogliaBrainCounts.txt')
publishedCountsGSM <- sapply(colnames(publishedCounts),function(x) unique(as.character(mapping$GSM[mapping$title==x])))
publishedCounts <- publishedCounts[,publishedCountsGSM%in%allGSM]
colnames(publishedCounts) <- publishedCountsGSM[publishedCountsGSM%in%allGSM]

SalmonTPMGeneLevel <- SalmonTPM$geneLevel$abundance
commonSamples <- intersect(colnames(publishedCounts),colnames(SalmonTPMGeneLevel))
commonGenes <- intersect(rownames(publishedCounts),rownames(SalmonTPMGeneLevel))
publishedCounts <- publishedCounts[commonGenes,commonSamples]
SalmonTPMGeneLevel <- SalmonTPMGeneLevel[commonGenes,commonSamples]

#Borrow gene effective length from Salmon
geneLengths <- SalmonTPM$geneLevel$length[commonGenes,commonSamples]
publishedTPM <- do.call(cbind, lapply(1:ncol(publishedCounts), function(i) {
  rate = log(publishedCounts[,i]) - log(geneLengths[,i])
  denom = log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}))
colnames(publishedTPM) <- colnames(publishedCounts)
SalmonAgainstSTAR <- GetGenePairs(SalmonTPMGeneLevel,publishedTPM)
PlotMethodCor(SalmonAgainstSTAR,'Salmon log(TPM+1)','Published STAR log(TPM+1)','Salmon TPM vs Published Converted TPM')

