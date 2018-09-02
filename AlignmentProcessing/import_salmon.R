##########################################################################################################
#Functions to import Salmon abundance using tximport. 
#Input should be path to directory where each subdirectory is a sample, containing a quant.sf
##########################################################################################################

library(GenomicFeatures)
library(tximport)
library(readr)
library(ggplot2)

#Build txdb from GTF file. Used GRCh37, Ensembl 75
ImportTx2gene <- function(){
  #Make database from GTF file
  txdb <- makeTxDbFromGFF('/local/data/public/zmx21/zmx21_private/GSK/GRCh37_Ensembl75/Homo_sapiens.GRCh37.75.gtf',format = 'gtf')
  k <- keys(txdb, keytype = "TXNAME")
  #Get all genes in the datbase as a dataframe.
  tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
  return(tx2gene)
}
#Use tximport to get gene and transcript level counts, from salmon output.
ImportSalmonCounts <- function(path,tx2gene){
  allDir <- dir(path)
  runDir <- allDir[union(union(grep('SRR',allDir),grep('GSM',allDir)),grep('H5KNCADXX',allDir))]
  #Get paths for all samples
  allPaths <- sapply(runDir,function(x) (paste0(path,'/',x,'/quant.sf')))
  
  #TPM counts at transcript (isoform) level
  transcriptLevelMat <- tximport(allPaths, type = "salmon", txOut = T)
  
  #TPM counts at gene level
  geneLevelMat <- summarizeToGene(transcriptLevelMat, tx2gene)
  return(list(geneLevel = geneLevelMat,transcriptLevel = transcriptLevelMat))
}
