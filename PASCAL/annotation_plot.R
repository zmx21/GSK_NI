library(egg)
load('../../GWAS/PASCAL_results/microglia_gene/ClusterAnnot_AB42_Deming.rda')
librariesToRun <- c('KEGG_2016','TRANSFAC_and_JASPAR_PWMs','Jensen_DISEASES')
AB42DemingDf <- currentStudyResult$df
overlap <- sapply(AB42DemingDf$KEGG_2016.Overlap,function(x) as.numeric(unlist(strsplit(x,'/'))[1])) #/AB42DemingDf$Size
AB42_Deming_Plots <- lapply(librariesToRun,function(x) ggplot(AB42DemingDf) + aes_string(x='adjPvalue',y=paste0(x,'.Adjusted.P.value'),colour='Size') + geom_point() +
  xlab('GWAS adj P Value') + ylab('Annotation P Value') + ggtitle(x))
tiff(filename = '../../Figures/PASCAL/GOAnno/AB42_Deming.tiff',width = 1600,height = 600)
ggarrange(plots=AB42_Deming_Plots,ncol = 3)
dev.off()

load('../../GWAS/PASCAL_results/microglia_gene/ClusterAnnot_tau_Deming.rda')
tau_DemingDf <- currentStudyResult$df
tau_Deming_Plots <- lapply(librariesToRun,function(x) ggplot(tau_DemingDf) + aes_string(x='adjPvalue',y=paste0(x,'.Adjusted.P.value'),colour='Size') + geom_point() +
                              xlab('GWAS adj P Value') + ylab('Annotation P Value') + ggtitle(x) + xlim(0,1))
tiff(filename = '../../Figures/PASCAL/GOAnno/tau_Deming.tiff',width = 1600,height = 600)
ggarrange(plots=tau_Deming_Plots,ncol = 3)
dev.off()

load('../../GWAS/PASCAL_results/microglia_gene/ClusterAnnot_MS_IMSGC_2013.rda')
MS_IMSGC_2013Df <- currentStudyResult$df
MS_IMSGC_2013_Plots <- lapply(librariesToRun,function(x) ggplot(MS_IMSGC_2013Df) + aes_string(x='adjPvalue',y=paste0(x,'.Adjusted.P.value'),colour='Size') + geom_point() +
                             xlab('GWAS adj P Value') + ylab('Annotation P Value') + ggtitle(x)+ xlim(0,1))
tiff(filename = '../../Figures/PASCAL/GOAnno/MS_IMSGC_2013.tiff',width = 1600,height = 600)
ggarrange(plots=MS_IMSGC_2013_Plots,ncol = 3)
dev.off()

load('../../GWAS/PASCAL_results/microglia_gene/ClusterAnnot_ptau_Deming.rda')
ptau_DemingDf <- currentStudyResult$df
ptau_Deming_Plots <- lapply(librariesToRun,function(x) ggplot(ptau_DemingDf) + aes_string(x='adjPvalue',y=paste0(x,'.Adjusted.P.value'),colour='Size') + geom_point() +
                                xlab('GWAS adj P Value') + ylab('Annotation P Value') + ggtitle(x)+ xlim(0,1))
tiff(filename = '../../Figures/PASCAL/GOAnno/ptau_Deming.tiff',width = 1600,height = 600)
ggarrange(plots=ptau_Deming_Plots,ncol = 3)
dev.off()

load('../../GWAS/PASCAL_results/microglia_gene/ClusterAnnot_MS_IMSGC_2011.rda')
MS_IMSGC_2011Df <- currentStudyResult$df
MS_IMSGC_2011_Plots <- lapply(librariesToRun,function(x) ggplot(MS_IMSGC_2011Df) + aes_string(x='adjPvalue',y=paste0(x,'.Adjusted.P.value'),colour='Size') + geom_point() +
                                xlab('GWAS adj P Value') + ylab('Annotation P Value') + ggtitle(x)+ xlim(0,1))
tiff(filename = '../../Figures/PASCAL/GOAnno/MS_IMSGC_2011.tiff',width = 1600,height = 600)
ggarrange(plots=MS_IMSGC_2011_Plots,ncol = 3)
dev.off()

load('../../GWAS/PASCAL_results/microglia_gene/ClusterAnnot_Alzheimer_IGAP_Stage1.rda')
Alzheimer_IGAP_Stage1Df <- currentStudyResult$df
Alzheimer_IGAP_Stage1_Plots <- lapply(librariesToRun,function(x) ggplot(Alzheimer_IGAP_Stage1Df) + aes_string(x='adjPvalue',y=paste0(x,'.Adjusted.P.value'),colour='Size') + geom_point() +
                                xlab('GWAS adj P Value') + ylab('Annotation P Value') + ggtitle(x)+ xlim(0,1))
tiff(filename = '../../Figures/PASCAL/GOAnno/Alzheimer_IGAP_Stage1.tiff',width = 1600,height = 600)
ggarrange(plots=Alzheimer_IGAP_Stage1_Plots,ncol = 3)
dev.off()


