##########################################################################################################
# Fig 3.2.1.1 
##########################################################################################################
#KEGG P value calculation. Remove those with no hits and those with overlap of 1.
library(dplyr)
#WGCNA Coding Genes
library(dplyr)
WGCNACodingKEGG <- readRDS('../../Count_Data/PASCAL_Results/Annotations/Annotation_JoinedDfMicrogliaWGCNAUnsigned_coding.rds')%>% 
{do.call(rbind,lapply(.,function(x) x$df))} %>%dplyr::distinct(Name,.keep_all=T) %>% dplyr::select(KEGG_P=KEGG_2016.Adjusted.P.value,KEGG_Overlap=KEGG_2016.Overlap) %>%
  dplyr::mutate(Method='WGCNA',Biotype='coding') %>% dplyr::mutate(KEGG_Overlap=sapply(KEGG_Overlap,function(x) strsplit(x=x,split = '[/]')[[1]][1]))
WGCNACodingKEGG$KEGG_Overlap[is.na(WGCNACodingKEGG$KEGG_Overlap)] <- 0
WGCNACodingKEGG$KEGG_P[is.na(WGCNACodingKEGG$KEGG_P)] <- 1
WGCNACodingKEGG <- WGCNACodingKEGG %>% dplyr::filter(as.numeric(KEGG_Overlap)>1) 

#WGCNA All Genes
WGCNAAllKEGG <- readRDS('../../Count_Data/PASCAL_Results/Annotations/Annotation_JoinedDfMicrogliaWGCNAUnsigned_all.rds')%>% 
{do.call(rbind,lapply(.,function(x) x$df))} %>%dplyr::distinct(Name,.keep_all=T) %>% dplyr::select(KEGG_P=KEGG_2016.Adjusted.P.value,KEGG_Overlap=KEGG_2016.Overlap) %>%
  dplyr::mutate(Method='WGCNA',Biotype='all') %>% dplyr::mutate(KEGG_Overlap=sapply(KEGG_Overlap,function(x) strsplit(x=x,split = '[/]')[[1]][1])) 
WGCNAAllKEGG$KEGG_Overlap[is.na(WGCNAAllKEGG$KEGG_Overlap)] <- 0
WGCNAAllKEGG$KEGG_P[is.na(WGCNAAllKEGG$KEGG_P)] <- 1
WGCNAAllKEGG <- WGCNAAllKEGG %>% dplyr::filter(as.numeric(KEGG_Overlap)>1) 

#Pearson Coding Genes
PearsonCodingKEGG <- readRDS('../../Count_Data/PASCAL_Results/Annotations/Annotation_Microglia_Pearson_cor0p2_abs_coding.rds') %>% 
{do.call(rbind,lapply(.,function(x) x$df))} %>%dplyr::distinct(Name,.keep_all=T) %>% dplyr::select(KEGG_P=KEGG_2016.Adjusted.P.value,KEGG_Overlap=KEGG_2016.Overlap) %>%
  dplyr::mutate(Method='Pearson',Biotype='coding') %>% dplyr::mutate(KEGG_Overlap=sapply(KEGG_Overlap,function(x) strsplit(x=x,split = '[/]')[[1]][1]))
PearsonCodingKEGG$KEGG_Overlap[is.na(PearsonCodingKEGG$KEGG_Overlap)] <- 0
PearsonCodingKEGG$KEGG_P[is.na(PearsonCodingKEGG$KEGG_P)] <- 1
PearsonCodingKEGG <- PearsonCodingKEGG %>% dplyr::filter(as.numeric(KEGG_Overlap)>1) 

#Pearson All Genes
PearsonAllKEGG <- readRDS('../../Count_Data/PASCAL_Results/Annotations/Annotation_Microglia_Pearson_cor0p2_abs_all.rds') %>% 
{do.call(rbind,lapply(.,function(x) x$df))} %>%dplyr::distinct(Name,.keep_all=T) %>% dplyr::select(KEGG_P=KEGG_2016.Adjusted.P.value,KEGG_Overlap=KEGG_2016.Overlap) %>%
  dplyr::mutate(Method='Pearson',Biotype='all') %>% dplyr::mutate(KEGG_Overlap=sapply(KEGG_Overlap,function(x) strsplit(x=x,split = '[/]')[[1]][1]))
PearsonAllKEGG$KEGG_Overlap[is.na(PearsonAllKEGG$KEGG_Overlap)] <- 0
PearsonAllKEGG$KEGG_P[is.na(PearsonAllKEGG$KEGG_P)] <- 1
PearsonAllKEGG <- PearsonAllKEGG %>% dplyr::filter(as.numeric(KEGG_Overlap)>1) 

#WGCNA Coding Genes
WGCNARandomCodingKEGG <- readRDS('../../Count_Data/PASCAL_Results/Annotations/Annotation_PermutedWGCNA_coding.rds') %>% 
  dplyr::select(KEGG_P=KEGG_2016.Adjusted.P.value,KEGG_Overlap=KEGG_2016.Overlap) %>%
  dplyr::mutate(Method='WGCNA Permuted',Biotype='coding') %>% dplyr::mutate(KEGG_Overlap=sapply(KEGG_Overlap,function(x) strsplit(x=x,split = '[/]')[[1]][1])) 
WGCNARandomCodingKEGG$KEGG_Overlap[is.na(WGCNARandomCodingKEGG$KEGG_Overlap)] <- 0
WGCNARandomCodingKEGG$KEGG_P[is.na(WGCNARandomCodingKEGG$KEGG_P)] <- 1
WGCNARandomCodingKEGG <- WGCNARandomCodingKEGG %>% dplyr::filter(as.numeric(KEGG_Overlap)>1) 

#Permuted WGCNA All Genes
WGCNARandomAllKEGG <- readRDS('../../Count_Data/PASCAL_Results/Annotations/Annotation_PermutedWGCNA_all.rds') %>% 
  dplyr::select(KEGG_P=KEGG_2016.Adjusted.P.value,KEGG_Overlap=KEGG_2016.Overlap) %>%
  dplyr::mutate(Method='WGCNA Permuted',Biotype='all') %>% dplyr::mutate(KEGG_Overlap=sapply(KEGG_Overlap,function(x) strsplit(x=x,split = '[/]')[[1]][1]))
WGCNARandomAllKEGG$KEGG_Overlap[is.na(WGCNARandomAllKEGG$KEGG_Overlap)] <- 0
WGCNARandomAllKEGG$KEGG_P[is.na(WGCNARandomAllKEGG$KEGG_P)] <- 1
WGCNARandomAllKEGG <- WGCNARandomAllKEGG %>% dplyr::filter(as.numeric(KEGG_Overlap)>1) 

#Permuted WGCNA Coding Transcripts
WGCNATranscriptRandomCodingKEGG <- readRDS('../../Count_Data/PASCAL_Results/Annotations/Annotation_PermutedWGCNATranscript_coding.rds')%>%
  dplyr::distinct(Name,.keep_all=T) %>% dplyr::select(KEGG_P=KEGG_2016.Adjusted.P.value,KEGG_Overlap=KEGG_2016.Overlap) %>%
  dplyr::mutate(Method='WGCNA Permuted',Biotype='coding transcript') %>% dplyr::mutate(KEGG_Overlap=sapply(KEGG_Overlap,function(x) strsplit(x=x,split = '[/]')[[1]][1]))
WGCNATranscriptRandomCodingKEGG$KEGG_Overlap[is.na(WGCNATranscriptRandomCodingKEGG$KEGG_Overlap)] <- 0
WGCNATranscriptRandomCodingKEGG$KEGG_P[is.na(WGCNATranscriptRandomCodingKEGG$KEGG_P)] <- 1
WGCNATranscriptRandomCodingKEGG <- WGCNATranscriptRandomCodingKEGG %>% dplyr::filter(as.numeric(KEGG_Overlap)>1) 

# WGCNA Coding Transcripts
WGCNATranscriptCodingKEGG <- readRDS('../../Count_Data/PASCAL_Results/Annotations/Annotation_JoinedDfMicrogliaWGCNATranscripts_coding.rds')  %>% 
{do.call(rbind,lapply(.,function(x) x$df))} %>% 
  dplyr::select(KEGG_P=KEGG_2016.Adjusted.P.value,KEGG_Overlap=KEGG_2016.Overlap) %>% dplyr::mutate(Method='WGCNA',Biotype='coding transcript') %>% dplyr::mutate(KEGG_Overlap=sapply(KEGG_Overlap,function(x) strsplit(x=x,split = '[/]')[[1]][1]))
WGCNATranscriptCodingKEGG$KEGG_Overlap[is.na(WGCNATranscriptCodingKEGG$KEGG_Overlap)] <- 0
WGCNATranscriptCodingKEGG$KEGG_P[is.na(WGCNATranscriptCodingKEGG$KEGG_P)] <- 1
WGCNATranscriptCodingKEGG <- WGCNATranscriptCodingKEGG %>% dplyr::filter(as.numeric(KEGG_Overlap)>1) 

# Pearson Coding Transcripts
PearsonTranscriptCodingKEGG <- readRDS('../../Count_Data/PASCAL_Results/Annotations/Annotation_JoinedDfMicrogliaPearsonTranscripts_coding.rds')%>%
  dplyr::distinct(Name,.keep_all=T) %>% dplyr::select(KEGG_P=KEGG_2016.Adjusted.P.value,KEGG_Overlap=KEGG_2016.Overlap) %>%
  dplyr::mutate(Method='Pearson',Biotype='coding transcript') %>% dplyr::mutate(KEGG_Overlap=sapply(KEGG_Overlap,function(x) strsplit(x=x,split = '[/]')[[1]][1]))
PearsonTranscriptCodingKEGG$KEGG_Overlap[is.na(PearsonTranscriptCodingKEGG$KEGG_Overlap)] <- 0
PearsonTranscriptCodingKEGG$KEGG_P[is.na(PearsonTranscriptCodingKEGG$KEGG_P)] <- 1
PearsonTranscriptCodingKEGG <- PearsonTranscriptCodingKEGG %>% dplyr::filter(as.numeric(KEGG_Overlap)>1) 

#Permuted Pearson Coding Genes
PearsonRandomCodingKEGG <- readRDS('../../Count_Data/PASCAL_Results/Annotations/Annotation_PermutedPearson_coding.rds') %>% 
  dplyr::select(KEGG_P=Adjusted.P.value,KEGG_Overlap=Overlap) %>%
  dplyr::mutate(Method='Pearson Permuted',Biotype='coding') %>% dplyr::mutate(KEGG_Overlap=sapply(KEGG_Overlap,function(x) strsplit(x=x,split = '[/]')[[1]][1])) 
PearsonRandomCodingKEGG$KEGG_Overlap[is.na(PearsonRandomCodingKEGG$KEGG_Overlap)] <- 0
PearsonRandomCodingKEGG$KEGG_P[is.na(PearsonRandomCodingKEGG$KEGG_P)] <- 1
PearsonRandomCodingKEGG <- PearsonRandomCodingKEGG %>% dplyr::filter(as.numeric(KEGG_Overlap)>1) 

#Permuted Pearson All Genes
PearsonRandomAllKEGG <- readRDS('../../Count_Data/PASCAL_Results/Annotations/Annotation_PermutedPearsonAll_all.rds') %>% 
  dplyr::select(KEGG_P=KEGG_2016.Adjusted.P.value,KEGG_Overlap=KEGG_2016.Overlap) %>%
  dplyr::mutate(Method='Pearson Permuted',Biotype='all') %>% dplyr::mutate(KEGG_Overlap=sapply(KEGG_Overlap,function(x) strsplit(x=x,split = '[/]')[[1]][1])) 
PearsonRandomCodingKEGG$KEGG_Overlap[is.na(PearsonRandomCodingKEGG$KEGG_Overlap)] <- 0
PearsonRandomCodingKEGG$KEGG_P[is.na(PearsonRandomCodingKEGG$KEGG_P)] <- 1
PearsonRandomCodingKEGG <- PearsonRandomCodingKEGG %>% dplyr::filter(as.numeric(KEGG_Overlap)>1) 

PearsonRandomTranscriptKEGG <- readRDS('../../Count_Data/PASCAL_Results/Annotations/Annotation_PermutedPearsonTranscript_coding.rds') %>% 
  dplyr::select(KEGG_P=KEGG_2016.Adjusted.P.value,KEGG_Overlap=KEGG_2016.Overlap) %>%
  dplyr::mutate(Method='Pearson Permuted',Biotype='coding transcript') %>% dplyr::mutate(KEGG_Overlap=sapply(KEGG_Overlap,function(x) strsplit(x=x,split = '[/]')[[1]][1])) 
PearsonRandomTranscriptKEGG$KEGG_Overlap[is.na(PearsonRandomTranscriptKEGG$KEGG_Overlap)] <- 0
PearsonRandomTranscriptKEGG$KEGG_P[is.na(PearsonRandomTranscriptKEGG$KEGG_P)] <- 1
PearsonRandomTranscriptKEGG <- PearsonRandomTranscriptKEGG %>% dplyr::filter(as.numeric(KEGG_Overlap)>1) 

collapsedDf <- rbind(PearsonCodingKEGG,PearsonAllKEGG,WGCNACodingKEGG,WGCNAAllKEGG,
                     WGCNARandomCodingKEGG,WGCNARandomAllKEGG,
                     PearsonTranscriptCodingKEGG,WGCNATranscriptCodingKEGG,WGCNATranscriptRandomCodingKEGG,
                     PearsonRandomCodingKEGG,PearsonRandomAllKEGG,PearsonRandomTranscriptKEGG)
collapsedDf$Method <- factor(collapsedDf$Method,levels=c('Pearson','Pearson Permuted','WGCNA','WGCNA Permuted'))
collapsedDf$Biotype <- factor(collapsedDf$Biotype,levels=c('coding','all','coding transcript'))
library(ggplot2)
library(ggsignif)
library(ggpubr)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 10
cols = gg_color_hue(n)

#Boxplot
group.colors <- c('Pearson'=cols[1],'Pearson Permuted'=cols[2],'WGCNA'=cols[7],'WGCNA Permuted'=cols[6])
p1 <- ggplot(collapsedDf, aes(x=Biotype, y=-1*log10(KEGG_P),fill=Method)) +
  geom_boxplot(outlier.shape = NA) + ylim(0,3.5) + ylab('-log10(KEGG Pathway P-value)') + xlab('Biotype') + 
  theme(text = element_text(size=16),axis.text.x=element_text(angle = 90, vjust = 0.5)) + scale_x_discrete(breaks = c('coding','all','coding transcript'),
                                                                                                           labels=c("Coding \n Genes", "Coding \n and lncRNA \n Genes","Coding \n Transcripts")) + 
  scale_fill_manual('Method',values=group.colors,
                    labels=c("\nPearson\n",
                             '\nPearson \nPermuted\n',
                             '\nWGCNA\n',
                             '\nWGCNA \nPermuted\n'),
                    breaks=c("Pearson",
                             'Pearson Permuted',
                             'WGCNA',
                             'WGCNA Permuted'))
wcTest1 <- compare_means(KEGG_P ~ Method, data = collapsedDf %>% dplyr::filter(Biotype=='coding' & (Method=='Pearson'|Method=='WGCNA')),paired = F)$p.format
wcTest2 <- compare_means(KEGG_P ~ Method, data = collapsedDf %>% dplyr::filter(Biotype=='all'& (Method=='Pearson'|Method=='WGCNA')),paired = F)$p.format
wcTest3 <- compare_means(KEGG_P ~ Method, data = collapsedDf %>% dplyr::filter(Biotype=='coding transcript'& (Method=='Pearson'|Method=='WGCNA')),paired = F)$p.format

wcTest4 <- compare_means(KEGG_P ~ Method, data = collapsedDf %>% dplyr::filter(Biotype=='coding'& (Method=='Pearson'|Method=='Pearson Permuted')),paired = F)$p.format
wcTest5 <- compare_means(KEGG_P ~ Method, data = collapsedDf %>% dplyr::filter(Biotype=='all'& (Method=='Pearson'|Method=='Pearson Permuted')),paired = F)$p.format
wcTest6 <- compare_means(KEGG_P ~ Method, data = collapsedDf %>% dplyr::filter(Biotype=='coding transcript'& (Method=='Pearson'|Method=='Pearson Permuted')),paired = F)$p.format
wcTest7 <- compare_means(KEGG_P ~ Method, data = collapsedDf %>% dplyr::filter(Biotype=='coding'& (Method=='WGCNA'|Method=='WGCNA Permuted')),paired = F)$p.format
wcTest8 <- compare_means(KEGG_P ~ Method, data = collapsedDf %>% dplyr::filter(Biotype=='all'& (Method=='WGCNA'|Method=='WGCNA Permuted')),paired = F)$p.format
wcTest9 <- compare_means(KEGG_P ~ Method, data = collapsedDf %>% dplyr::filter(Biotype=='coding transcript'& (Method=='WGCNA'|Method=='WGCNA Permuted')),paired = F)$p.format

#significance table
testTable <- data.frame(Biotype=c(rep('Coding\nGenes',3),rep('Coding and \n lncRNA Genes',3),rep('Coding\nTranscripts',3)),
                        Comparision=rep(c('Pearson vs \nWGCNA','Pearson vs \n Permuted','WGCNA vs \n Permuted'),3))


testTable$Wilcoxon_P <- c(wcTest1,wcTest4,wcTest7,wcTest2,wcTest5,wcTest8,wcTest3,wcTest6,wcTest9)
colnames(testTable) <- c('Biotype','Method \n Comparison','Mann-Whitney\nU Test\n P-value')
rownames(testTable) <- c()

library(gridExtra)
library(grid)
library(gtable)

g <- tableGrob(testTable)

g <- gtable_add_grob(g,
                     grobs = segmentsGrob( # line across the bottom
                       x0 = unit(0,"npc"),
                       y0 = unit(0,"npc"),
                       x1 = unit(1,"npc"),
                       y1 = unit(0,"npc"),
                       gp = gpar(lwd = 2.0)),
                     t = 4, b = 4, l = 2, r = 4)
g <- gtable_add_grob(g,
                     grobs = segmentsGrob( # line across the bottom
                       x0 = unit(0,"npc"),
                       y0 = unit(0,"npc"),
                       x1 = unit(1,"npc"),
                       y1 = unit(0,"npc"),
                       gp = gpar(lwd = 2.0)),
                     t = 7, b = 7, l = 2, r = 4)


ggpubr::ggexport(ggpubr::ggarrange(p1,g,ncol = 2,widths = c(2,1.5)),filename = '../../FinalFigures/KEGGComparison.pdf',width = 10,height = 7)

Combinations <- expand.grid(Method=unique(collapsedDf$Method),Biotype=unique(collapsedDf$Biotype))
summarizedData <- data.frame(Method=Combinations$Method,Biotype=Combinations$Biotype)
summarizedStats <- data.frame(FractSig=rep(NA,nrow(summarizedData)),NumClusters=rep(NA,nrow(summarizedData)))

for(i in 1:nrow(summarizedStats)){
  curDf <- dplyr::filter(collapsedDf,Method==Combinations$Method[i],Biotype==Combinations$Biotype[i])
  summarizedStats$FractSig[i] <- sum(curDf$KEGG_P<0.05)/nrow(curDf)
  summarizedStats$NumClusters[i] <- nrow(curDf)
}
summarizedData <- cbind(summarizedData,summarizedStats)

p2 <- ggplot(summarizedData, aes(x=Biotype, y=FractSig,fill=Method)) +
  geom_bar(stat = 'identity',position = 'dodge') + ylab('Fraction of Significant KEGG Enrichment') + xlab('Biotype') + 
  theme(text = element_text(size=16),axis.text.x=element_text(angle = 90, vjust = 0.5)) + scale_x_discrete(breaks = c('coding','all','coding transcript'),
                                                                                                           labels=c("Coding \n Genes", "Coding \n and lncRNA \n Genes","Coding \n Transcripts"))
