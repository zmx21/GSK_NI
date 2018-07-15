library(data.table)
library(dplyr)
# source('annotate_clusters.R')
GetGWASMetdata <- function(sigClusters){
  GWASMetadata <- read.csv('/local/data/public/zmx21/zmx21_private/GSK/GWAS/GWASMetadata.csv',stringsAsFactors = F)
  
  PASCAL_raw_Path <- '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_results/KEGG_Ensembl/'
  PASCAL_results_raw <- dir(PASCAL_raw_Path)
  PASCAL_results_raw <- PASCAL_results_raw[grep('sum.genescores',PASCAL_results_raw)]
  numSigGenes <- sapply(PASCAL_results_raw,function(x) data.table::fread(paste0(PASCAL_raw_Path,x)) %>%
                          dplyr::filter(pvalue < 1e-8) %>% nrow)
  
  PASCAL_cluster_coding_path <- '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_results2/'
  PASCAL_results_coding <- dir(PASCAL_cluster_coding_path)
  PASCAL_results_coding <- PASCAL_results_coding[grep('sum.genescores',PASCAL_results_coding)]
  numSigGenesIncludedCoding <- sapply(PASCAL_results_coding,function(x) data.table::fread(paste0(PASCAL_cluster_coding_path,x)) %>% 
                                  dplyr::filter(pvalue < 1e-8) %>% nrow)
  sigGenesCoding <- lapply(PASCAL_results_coding,function(x) data.table::fread(paste0(PASCAL_cluster_coding_path,x)) %>% 
                       dplyr::filter(pvalue < 1e-8) %>% select(gene_id) %>% t() %>% as.vector)
  
  SigStudiesStats <- data.frame(N_Sig_Genes = numSigGenes,N_Sig_Genes_In_Network=numSigGenesIncluded)
  SigStudiesStats <- data.frame(N_Sig_Genes_In_Network=numSigGenesIncluded)
  SigStudiesStats$StudyName <- sapply(rownames(SigStudiesStats),function(x) unlist(strsplit(x,'.sum'))[1])
  SigStudiesStats$StudyName <- sapply(SigStudiesStats$StudyName,function(x) gsub('[.]','_',x))
  
  # sigClusters <- do.call(rbind,lapply(GetAnnotationForSignificantClusters(),function(x) x$df))
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
  GWASFullMetadata <- GWASFullMetadata[,c(2,1,3)]
  # GWASFullMetadata[[is.na(GWASFullMetadata[,3]),3] <- 0
  
  library(gridExtra)
  grid.arrange(tableGrob(GWASFullMetadata))
  
}

