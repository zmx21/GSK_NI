library(data.table)
library(dplyr)
source('annotate_clusters.R')
GWASMetadata <- read.csv('/local/data/public/zmx21/zmx21_private/GSK/GWAS/GWASMetadata.csv',stringsAsFactors = F)

PASCAL_nofilt_Path <- '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_results/Old/KEGG_Ensembl/'
PASCAL_results_nofilt <- dir(PASCAL_nofilt_Path)
PASCAL_results_nofilt <- PASCAL_results_nofilt[grep('sum.genescores',PASCAL_results_nofilt)]
numSigGenes <- sapply(PASCAL_results_nofilt,function(x) data.table::fread(paste0(PASCAL_nofilt_Path,x)) %>% 
                        dplyr::filter(pvalue < 1e-8) %>% nrow)

PASCAL_pathway_path <- '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_results/microglia_all_genes/'
PASCAL_results_pathway <- dir(PASCAL_pathway_path)
PASCAL_results_pathway <- PASCAL_results_pathway[grep('sum.genescores',PASCAL_results_pathway)]
numSigGenesIncluded <- sapply(PASCAL_results_pathway,function(x) data.table::fread(paste0(PASCAL_pathway_path,x)) %>% 
                        dplyr::filter(pvalue < 1e-8) %>% nrow)
sigGenes <- lapply(PASCAL_results_pathway,function(x) data.table::fread(paste0(PASCAL_pathway_path,x)) %>% 
                     dplyr::filter(pvalue < 1e-8) %>% select(gene_id) %>% t() %>% as.vector)

SigStudiesStats <- data.frame(N_Sig_Genes = numSigGenes,N_Sig_Genes_In_Network=numSigGenesIncluded)
SigStudiesStats$StudyName <- sapply(rownames(SigStudiesStats),function(x) unlist(strsplit(x,'.sum'))[1])
SigStudiesStats$StudyName <- sapply(SigStudiesStats$StudyName,function(x) gsub('[.]','_',x))

sigClusters <- do.call(rbind,lapply(GetAnnotationForSignificantClusters(),function(x) x$df))
numSigClusters <- sigClusters %>% dplyr::count(StudyName)
numSigClusters$n[is.na(numSigClusters$n)] <- 0

numSigGenesInSigClusters <- rep(NA,length(sigGenes))
for(i in 1:length(sigGenes)){
  currentStudy <- SigStudiesStats$StudyName[i]
  currentStudySigClusters <- dplyr::filter(sigClusters,StudyName==SigStudiesStats$StudyName[i])
  currentStudyIncludedGenes <- unique(unlist(currentStudySigClusters$Genes))
  numSigGenesInSigClusters[i] <- length(intersect(sigGenes[[i]],currentStudyIncludedGenes))
}
SigStudiesStats$N_SigGenesInSigClusters <- numSigGenesInSigClusters

GWASFullMetadata <- dplyr::left_join(SigStudiesStats,GWASMetadata,by=c('StudyName'='StudyName')) %>% 
  dplyr::left_join(numSigClusters,by=c('StudyName'='StudyName')) %>% rename_('N_Subjects'='N','N_Sig_Clusters'='n')
GWASFullMetadata <- GWASFullMetadata[,c(3,5,1,2,4,6)]
GWASFullMetadata$N_Sig_Clusters[is.na(GWASFullMetadata$N_Sig_Clusters)] <- 0

library(gridExtra)
grid.arrange(tableGrob(GWASFullMetadata))

