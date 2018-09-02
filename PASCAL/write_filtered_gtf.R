#Filter GTF such that only genes in our network are considered
#This is the same approach which the DREAM network consortium used, such that gene background is not the full genome. 
library(data.table)
library(dplyr)
GTFPath <- '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_New/resources/annotation/'
FullGTF <- data.table::fread(paste0(GTFPath,'Homo_sapiens.GRCh37.75.gencodeformat.gtf')) %>% 
  dplyr::filter(V1 %in% sapply(seq(1,22,1),function(x) paste0('chr',as.character(x))))
GeneID <- sapply(FullGTF$V9,function(x){ 
  semiColonSplit <- unlist(strsplit(x,';'))[1];
  geneID <- unlist(strsplit(semiColonSplit,'\"'))[2];
  return(geneID)})
names(GeneID) <- c()

#Load genes which are under consideration.
load('../../Count_Data/CV_Filtered/MicrogliaGeneCVFiltered.rda')
microgliaCodingGenesToInclude <- rownames(MicrogliaGeneCVFiltered$coding)
microgliaAllGenesToInclude <- c(microgliaCodingGenesToInclude,rownames(MicrogliaGeneCVFiltered$noncoding))

microgliaCodingGenesGTF <- FullGTF[GeneID %in% microgliaCodingGenesToInclude,]
microgliaAllGenesGTF <- FullGTF[GeneID %in% microgliaAllGenesToInclude,]


data.table::fwrite(microgliaCodingGenesGTF,file = paste0(GTFPath,'MicrogliaCodingGenes.gtf'),
                   row.names = F,col.names = F,quote = F,sep = '\t')

data.table::fwrite(microgliaAllGenesGTF,file = paste0(GTFPath,'MicrogliaCodinglncRNAGenes.gtf'),
                   row.names = F,col.names = F,quote = F,sep = '\t')

#Generate GTF files
microgliaCodingTranscriptsToInclude <- unlist(qusage::read.gmt(
  '/local/data/public/zmx21/zmx21_private/GSK/Louvain_results/Pearson_Cor0p2/CodingTranscriptsMicroglia_Pearson_cor0p2_abs.gmt')) %>% unique()
microgliaCodingTranscriptGTF <- FullGTF[GeneID %in% microgliaCodingTranscriptsToInclude,]
data.table::fwrite(microgliaCodingTranscriptGTF,file=paste0(GTFPath,'MicrogliaCodingTranscripts.gtf'),
                   row.names = F,col.names = F,quote = F,sep = '\t')

microgliaAllTranscriptsToInclude <- unlist(qusage::read.gmt(
  '/local/data/public/zmx21/zmx21_private/GSK/Louvain_results/Pearson_Cor0p2/AllTranscriptsMicroglia_Pearson_cor0p2_abs.gmt')) %>% unique()
microgliaAllTranscriptGTF <- FullGTF[GeneID %in% microgliaAllTranscriptsToInclude,]
data.table::fwrite(microgliaAllTranscriptGTF,file=paste0(GTFPath,'MicrogliaCodinglncRNATranscripts.gtf'),
                   row.names = F,col.names = F,quote = F,sep = '\t')




