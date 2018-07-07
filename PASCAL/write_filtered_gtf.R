#Filter GTF such that only genes in our network are considered
#This is the same approach which the DREAM network consortium used, such that gene background is not the full genome. 
library(data.table)
library(dplyr)
GTFPath <- '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_Ensembl/resources/annotation/'
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

data.table::fwrite(microgliaAllGenesGTF,file = paste0(GTFPath,'MicrogliaAllGenes.gtf'),
                   row.names = F,col.names = F,quote = F,sep = '\t')
