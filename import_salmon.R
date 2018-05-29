library(GenomicFeatures)
library(tximport)
library(readr)
library(ggplot2)

#Build txdb from GTF file. Used GRCh37, Ensembl 75
ImportTx2gene <- function(){
  txdb <- makeTxDbFromGFF('/local/data/public/zmx21/zmx21_private/GSK/GRCh37_Ensembl75/Homo_sapiens.GRCh37.75.gtf',format = 'gtf')
  k <- keys(txdb, keytype = "TXNAME")
  tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
  return(tx2gene)
}

ImportSalmonCounts <- function(path,tx2gene){
  allDir <- dir(path)
  runDir <- allDir[union(grep('SRR',allDir),grep('GSM',allDir))]
  #Get paths for all samples
  allPaths <- sapply(runDir,function(x) (paste0(path,'/',x,'/quant.sf')))
  
  #TPM counts for transcript
  transcriptLevelMat <- tximport(allPaths, type = "salmon", txOut = T)
  
  #TPM counts for genes
  geneLevelMat <- summarizeToGene(transcriptLevelMat, tx2gene)
  return(list(geneLevel = geneLevelMat,transcriptLevel = transcriptLevelMat))
}
