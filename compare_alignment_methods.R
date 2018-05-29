source('import_salmon.R')
GetGenePairs <- function(method1,method2,logScale=F){
  geneNames <- intersect(rownames(method1),rownames(method2))
  sampleNames <- intersect(colnames(method1),colnames(method2))
  method1 <- method1[geneNames,]
  method2 <- method2[geneNames,]
  
  if(logScale){
    expressionPairs <- lapply(sampleNames,function(x) data.frame(method1=log(method1[,x] + 1),method2=log(method2[,x] + 1)))
  }else{
    expressionPairs <- lapply(sampleNames,function(x) data.frame(method1=method1[,x],method2=method2[,x]))
  }
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
# tx2gene <- ImportTx2gene()
# SalmonTPM <- ImportSalmonCounts('/local/data/public/zmx21/zmx21_private/GSK/Galatro/Salmon_aligned/Salmon_aligned_merged',tx2gene)

# Salmon_FMD_TPM <- ImportSalmonCounts('/local/data/public/zmx21/zmx21_private/GSK/Galatro/Salmon_aligned/Salmon_aligned_fmd',tx2gene)
# Salmon_Quasi_TPM <- ImportSalmonCounts('/local/data/public/zmx21/zmx21_private/GSK/Galatro/Salmon_aligned/Salmon_aligned_k19',tx2gene)
# FMD_Against_Quasi <- GetGenePairs(Salmon_FMD_TPM$geneLevel$abundance,Salmon_Quasi_TPM$geneLevel$abundance)
# PlotMethodCor(FMD_Against_Quasi,'FMD log(TPM + 1)','Quasi log(TPM + 1)','Salmon - FMD vs Quasi')

# SalmonAgainstSTAR <- GetGenePairs(SalmonTPMGeneLevel,publishedTPM,logScale=T)
# PlotMethodCor(SalmonAgainstSTAR,'Salmon log(TPM+1)','Published STAR log(TPM+1)','Salmon TPM vs Published Converted TPM')

