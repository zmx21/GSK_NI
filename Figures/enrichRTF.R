##########################################################################################################
# Fig 3.4.1.2 and Fig 3.4.2.1
##########################################################################################################
#Get module genes
cpg12Genes <- data.table::fread(file = '../../FinalFigures/CPG12_Genes.txt',header = F)$V1
cpg2Genes <- data.table::fread(file = '../../FinalFigures/CPG2_Genes.txt',header = F)$V1

library(enrichR)
library(dplyr)
#Run enrichR for databases for interest
librariesToRun <- c('WikiPathways_2016','ChEA_2016','Transcription_Factor_PPIs','PPI_Hub_Proteins')
cpg12EnrichR <- enrichr(cpg12Genes,databases = librariesToRun)
cpg2EnrichR <- enrichr(cpg2Genes,databases = librariesToRun)


library(grid)
library(gtable)
library(gridExtra)

#Construct ggtables
t1 <- ttheme_default(core=list(
  bg_params = list(fill=c("#FFFF99","grey90","grey95","grey90","grey95"))
))
WikiPathways <- tableGrob(cpg12EnrichR$WikiPathways_2016[1:5,] %>%
                            dplyr::mutate(Term= sapply(Term,function(x) gsub(x=x,pattern = '_WP.*',''))) %>%
  dplyr::mutate(Genes=sapply(Genes,function(x) gsub(x=x,pattern = ';',replacement = ' ')))%>% dplyr::select('WikiPathways Term'=Term,'Adjusted P-value' = Adjusted.P.value,'Enriched Genes'=Genes),theme = t1)
title <- textGrob("A) WikiPathways 2016", gp = gpar(fontsize = 16))
padding <- unit(0.5,"line")
WikiPathways <- gtable_add_rows(
  WikiPathways, heights = grobHeight(title) + padding, pos = 0
)
WikiPathways <- gtable_add_grob(
  WikiPathways, list(title),
  t = 1, l = 1, r = ncol(WikiPathways)
)


t1 <- ttheme_default(core=list(
  bg_params = list(fill=c("#FFFF99","grey90","grey95","grey90","grey95"))
))
chEA <- tableGrob(cpg12EnrichR$ChEA_2016[1:5,]%>%
  dplyr::mutate(Genes=sapply(Genes,function(x) gsub(x=x,pattern = ';',replacement = ' ')))%>% dplyr::select('ChEA Term'=Term,'Adjusted P-value' = Adjusted.P.value,'Enriched Genes'=Genes),theme = t1)
title <- textGrob("B) ChEA", gp = gpar(fontsize = 16))
padding <- unit(0.5,"line")
chEA <- gtable_add_rows(
  chEA, heights = grobHeight(title) + padding, pos = 0
)
chEA <- gtable_add_grob(
  chEA, list(title),
  t = 1, l = 1, r = ncol(chEA)
)


t1 <- ttheme_default(core=list(
  bg_params = list(fill=c("#FFFF99","grey90","grey95","grey90","grey95"))
))
TF <- tableGrob(cpg12EnrichR$Transcription_Factor_PPIs[1:5,]%>%
  dplyr::mutate(Genes=sapply(Genes,function(x) gsub(x=x,pattern = ';',replacement = ' '))) %>% dplyr::select('TF-PPI Term'=Term,'Adjusted P-value' = Adjusted.P.value,'Enriched Genes'=Genes),theme = t1)
title <- textGrob("C) Transcription Factor PPIs", gp = gpar(fontsize = 16))
padding <- unit(0.5,"line")
TF <- gtable_add_rows(
  TF, heights = grobHeight(title) + padding, pos = 0
)
TF <- gtable_add_grob(
  TF, list(title),
  t = 1, l = 1, r = ncol(TF)
)

t1 <- ttheme_default(core=list(
  bg_params = list(fill=c("#FFFF99","#FFFF99","#FFFF99","#FFFF99","#FFFF99"))
))
PPIHub <- tableGrob(cpg12EnrichR$PPI_Hub_Proteins[1:5,]%>%
  dplyr::mutate(Genes=sapply(Genes,function(x) gsub(x=x,pattern = ';',replacement = ' ')))%>% dplyr::select('PPI-Hub Term'=Term,'Adjusted P-value' = Adjusted.P.value,'Enriched Genes'=Genes),theme = t1)
title <- textGrob("D) PPI Hub Proteins", gp = gpar(fontsize = 16))
padding <- unit(0.5,"line")
PPIHub <- gtable_add_rows(
  PPIHub, heights = grobHeight(title) + padding, pos = 0
)
PPIHub <- gtable_add_grob(
  PPIHub, list(title),
  t = 1, l = 1, r = ncol(PPIHub)
)

library(ggpubr)
ggpubr::ggexport(ggpubr::ggarrange(WikiPathways,chEA,TF,PPIHub,ncol = 1,nrow=4),filename = '../../FinalFigures/TFAnalysis_JakStat.pdf',width = 10,height = 8)


cpg2EnrichR$WikiPathways_2016$Genes[[1]] <- "HLA-DRB5;HLA-DMB;HLA-B;\nHLA-DRA;HLA-DRB1;HLA-DQB1"
t1 <- ttheme_default(core=list(
  bg_params = list(fill=c("#FFFF99","grey90","grey95","grey90","grey95"))
))
WikiPathways <- tableGrob(cpg2EnrichR$WikiPathways_2016[1:5,] %>% dplyr::mutate(Term= sapply(Term,function(x) gsub(x=x,pattern = '_WP.*',''))) %>%
                            dplyr::mutate(Genes=sapply(Genes,function(x) gsub(x=x,pattern = ';',replacement = ' ')))%>% dplyr::select('WikiPathways Term'=Term,'Adjusted P-value' = Adjusted.P.value,'Enriched Genes'=Genes) ,theme = t1)
title <- textGrob("A) WikiPathways 2016", gp = gpar(fontsize = 16))
padding <- unit(0.5,"line")
WikiPathways <- gtable_add_rows(
  WikiPathways, heights = grobHeight(title) + padding, pos = 0
)
WikiPathways <- gtable_add_grob(
  WikiPathways, list(title),
  t = 1, l = 1, r = ncol(WikiPathways)
)

cpg2EnrichR$ChEA_2016$Genes[[1]] <- "HLA-DRB5;C19ORF81;HLA-DMB;\nS100A13;HLA-DRA;PTPRCAP;\nHLA-DRB1;PRRC2A;HLA-DQB1"
cpg2EnrichR$ChEA_2016$Genes[[4]] <- "CEBPA;ORM1;S100A1;\nS100A13;FFAR4;PDSS1;\nCES4A;ENKD1"
chEA <- tableGrob(cpg2EnrichR$ChEA_2016[1:5,]%>%
                    dplyr::mutate(Genes=sapply(Genes,function(x) gsub(x=x,pattern = ';',replacement = ' ')))%>% dplyr::select('ChEA Term'=Term,'Adjusted P-value' = Adjusted.P.value,'Enriched Genes'=Genes))
title <- textGrob("B) ChEA", gp = gpar(fontsize = 16))
padding <- unit(0.5,"line")
chEA <- gtable_add_rows(
  chEA, heights = grobHeight(title) + padding, pos = 0
)
chEA <- gtable_add_grob(
  chEA, list(title),
  t = 1, l = 1, r = ncol(chEA)
)
TF <- tableGrob(cpg2EnrichR$Transcription_Factor_PPIs[1:5,] %>%
                  dplyr::mutate(Genes=sapply(Genes,function(x) gsub(x=x,pattern = ';',replacement = ' ')))%>%dplyr::select('TF-PPI Term'=Term,'Adjusted P-value' = Adjusted.P.value,'Enriched Genes'=Genes))
title <- textGrob("C) Transcription Factor PPIs", gp = gpar(fontsize = 16))
padding <- unit(0.5,"line")
TF <- gtable_add_rows(
  TF, heights = grobHeight(title) + padding, pos = 0
)
TF <- gtable_add_grob(
  TF, list(title),
  t = 1, l = 1, r = ncol(TF)
)
PPIHub <- tableGrob(cpg2EnrichR$PPI_Hub_Proteins[1:5,]%>%
                      dplyr::mutate(Genes=sapply(Genes,function(x) gsub(x=x,pattern = ';',replacement = ' ')))%>%dplyr::select('PPI-Hub Term'=Term,'Adjusted P-value' = Adjusted.P.value,'Enriched Genes'=Genes))
title <- textGrob("D) PPI Hub Proteins", gp = gpar(fontsize = 16))
padding <- unit(0.5,"line")
PPIHub <- gtable_add_rows(
  PPIHub, heights = grobHeight(title) + padding, pos = 0
)
PPIHub <- gtable_add_grob(
  PPIHub, list(title),
  t = 1, l = 1, r = ncol(PPIHub)
)
ggpubr::ggexport(ggpubr::ggarrange(WikiPathways,chEA,TF,PPIHub,ncol = 1,nrow=4),filename = '../../FinalFigures/TFAnalysis_Allograft.pdf',width = 9,height = 10)
