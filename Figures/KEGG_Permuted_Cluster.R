##########################################################################################################
# Generate permuted modules for Fig 3.2.1.1 
##########################################################################################################
library(qusage)
library(dplyr)
source('../PASCAL/annotate_clusters.R')
set.seed(2)
permuted_pearson_coding <- read.gmt('../../GWAS/PASCAL_New/resources/genesets/Permuted_Pearson_100/permuted_Coding_Pearson.gmt')
permuted_pearson_all <- read.gmt('../../genesets/Permuted_Pearson_100/permuted_All_Pearson.gmt')
permuted_WGCNA_coding <- read.gmt('../../genesets/Permuted_WGCNA_100/permuted_Coding_WGCNA.gmt')
permuted_WGCNA_all <- read.gmt('../../genesets/Permuted_WGCNA_100/permuted_All_WGCNA.gmt')
permuted_WGCNA_coding <- data_frame(Name = names(permuted_WGCNA_coding),
                                      Genes=permuted_WGCNA_coding,Biotype=rep('coding',length(permuted_WGCNA_coding)))
permuted_WGCNA_all <- data_frame(Name = names(permuted_WGCNA_all),
                                    Genes=permuted_WGCNA_all,Biotype=rep('all',length(permuted_WGCNA_all)))
JoinedDfWGCNAPermutedFake <- rbind(permuted_WGCNA_coding,permuted_WGCNA_all)
JoinedDfWGCNAPermutedFake <- JoinedDfWGCNAPermutedFake %>%
  dplyr::mutate(RealName = sapply(Name,function(x) gsub(x=x,pattern = '_perm_.*',replacement = '')))
names(JoinedDfWGCNAPermutedFake$RealName) <- c()
#Merge FDR from real cluster to the permuted clusters
real_WGCNA_coding <- readRDS('../../Count_Data/PASCAL_Results/Annotations/Annotation_JoinedDfMicrogliaWGCNAUnsigned_coding.rds') %>% 
{do.call(rbind,lapply(.,function(x) x$df))}%>% dplyr::distinct(Name,.keep_all=T)%>% dplyr::mutate(KEGG_Overlap=sapply(KEGG_2016.Overlap,function(x) strsplit(x=x,split = '[/]')[[1]][1]))%>%
  dplyr::filter(as.numeric(KEGG_Overlap)>1) 
real_WGCNA_all <- readRDS('../../Count_Data/PASCAL_Results/Annotations/Annotation_JoinedDfMicrogliaWGCNAUnsigned_all.rds')%>% 
{do.call(rbind,lapply(.,function(x) x$df))} %>% dplyr::distinct(Name,.keep_all=T)%>% dplyr::mutate(KEGG_Overlap=sapply(KEGG_2016.Overlap,function(x) strsplit(x=x,split = '[/]')[[1]][1]))%>% 
  dplyr::filter(as.numeric(KEGG_Overlap)>1) 


JoinedDfWGCNAPermutedFake <- JoinedDfWGCNAPermutedFake %>% 
  dplyr::filter(RealName %in% c(real_WGCNA_coding$Name,real_WGCNA_all$Name)) %>%
  dplyr::select(Name=Name,Genes,RealName,Biotype)
uniqueRealNames <- unique(JoinedDfWGCNAPermutedFake$RealName)
JoinedDfWGCNAPermutedFakeTen <- JoinedDfWGCNAPermutedFake[0,]
for(i in 1:length(uniqueRealNames)){
  curDf <- JoinedDfWGCNAPermutedFake %>% dplyr::filter(RealName==uniqueRealNames[i]) %>% {.[1:50,]}
  JoinedDfWGCNAPermutedFakeTen <- rbind(JoinedDfWGCNAPermutedFakeTen,curDf)
}
JoinedDfWGCNAPermutedFakeTen <- dplyr::left_join(JoinedDfWGCNAPermutedFakeTen,rbind(real_WGCNA_coding,real_WGCNA_all),by=c('RealName'='Name')) %>%
  dplyr::select(Name=RealName,adjPvalue,Size=Size,StudyName,Biotype=Biotype.x,Genes=Genes.x) %>% dplyr::mutate(adjPvalue=1)
saveRDS(JoinedDfWGCNAPermutedFakeTen,'../../Count_Data/PASCAL_Results/permuted_fake_new/PermutedWGCNA.rds')


##########################TRANSCRIPT LEVEL########################################################
permuted_WGCNA_transcript_coding <- read.gmt('../../genesets/Permuted_WGCNA_Transcript/permuted_WGCNA_Coding_Transcript.gmt')
permuted_WGCNA_transcript_coding <- data_frame(Name = names(permuted_WGCNA_transcript_coding),
                                    Genes=permuted_WGCNA_transcript_coding,Biotype=rep('coding',length(permuted_WGCNA_transcript_coding)))

JoinedDfWGCNATranscriptPermutedFake <- permuted_WGCNA_transcript_coding
JoinedDfWGCNATranscriptPermutedFake <- JoinedDfWGCNATranscriptPermutedFake %>%
  dplyr::mutate(RealName = sapply(Name,function(x) gsub(x=x,pattern = '_perm_.*',replacement = '')))
names(JoinedDfWGCNATranscriptPermutedFake$RealName) <- c()
#Merge FDR from real cluster to the permuted clusters
real_WGCNATranscript_coding <- readRDS('../../Count_Data/PASCAL_Results/Annotations/Annotation_JoinedDfMicrogliaWGCNATranscripts_coding.rds') %>%
  {do.call(rbind,lapply(.,function(x) x$df))} %>% dplyr::distinct(Name,.keep_all=T)%>% dplyr::mutate(KEGG_Overlap=sapply(KEGG_2016.Overlap,function(x) strsplit(x=x,split = '[/]')[[1]][1])) %>% 
  dplyr::filter(as.numeric(KEGG_Overlap)>1) 

JoinedDfWGCNATranscriptPermutedFake <- JoinedDfWGCNATranscriptPermutedFake %>% 
  dplyr::filter(RealName %in% c(real_WGCNATranscript_coding$Name)) %>%
  dplyr::select(Name=Name,RealName,Genes,Biotype)
uniqueRealNames <- unique(JoinedDfWGCNATranscriptPermutedFake$RealName)
JoinedDfWGCNATranscriptPermutedFakeTen <- JoinedDfWGCNATranscriptPermutedFake[0,]
for(i in 1:length(uniqueRealNames)){
  curDf <- JoinedDfWGCNATranscriptPermutedFake %>% dplyr::filter(RealName==uniqueRealNames[i]) %>% {.[1:50,]}
  JoinedDfWGCNATranscriptPermutedFakeTen <- rbind(JoinedDfWGCNATranscriptPermutedFakeTen,curDf)
}
JoinedDfWGCNATranscriptPermutedFakeTen <- dplyr::left_join(JoinedDfWGCNATranscriptPermutedFakeTen,rbind(real_WGCNATranscript_coding),by=c('RealName'='Name')) %>%
  dplyr::select(Name=RealName,Perm=Name,adjPvalue,Size=Size,StudyName,Biotype=Biotype.x,Genes=Genes.x) %>% dplyr::mutate(adjPvalue=1)
saveRDS(JoinedDfWGCNATranscriptPermutedFakeTen,'../../Count_Data/PASCAL_Results/permuted_fake_new/PermutedWGCNATranscript.rds')

GetTopPercentileOfEachMethod(topPercentile = 1,methodPath = paste0('../../Count_Data/PASCAL_Results/permuted_fake_new/',c('PermutedWGCNA'),'.rds'),
                             method = c('PermutedWGCNA'))
GetTopPercentileOfEachMethod(topPercentile = 1,methodPath = paste0('../../Count_Data/PASCAL_Results/permuted_fake_new/',c('PermutedWGCNATranscript'),'.rds'),
                             method = c('PermutedWGCNATranscript'))
