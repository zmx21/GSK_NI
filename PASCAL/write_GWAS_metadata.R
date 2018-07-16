library(data.table)
library(dplyr)
# source('annotate_clusters.R')
GetGWASMetdata <- function(sigClusters=NULL,clusterLevel=T){
  GWASMetadata <- read.csv('/local/data/public/zmx21/zmx21_private/GSK/GWAS/GWASMetadata.csv',stringsAsFactors = F)
  
  if(!clusterLevel){
    #Get significant genes included in GWAS studies
    PASCAL_raw_Path <- '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_results/KEGG_Ensembl_Coding/'
    PASCAL_results_raw <- dir(PASCAL_raw_Path)
    PASCAL_results_raw <- PASCAL_results_raw[grep('sum.genescores',PASCAL_results_raw)]
    numSigGenesCoding <- sapply(PASCAL_results_raw,function(x) data.table::fread(paste0(PASCAL_raw_Path,x)) %>%
                                  dplyr::filter(pvalue < 5e-8) %>% nrow)
    PASCAL_raw_Path <- '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_results/KEGG_Ensembl_All/'
    PASCAL_results_raw <- dir(PASCAL_raw_Path)
    PASCAL_results_raw <- PASCAL_results_raw[grep('sum.genescores',PASCAL_results_raw)]
    numSigGenesAll <- sapply(PASCAL_results_raw,function(x) data.table::fread(paste0(PASCAL_raw_Path,x)) %>%
                                  dplyr::filter(pvalue < 5e-8) %>% nrow)
    StudyName <- sapply(names(numSigGenesCoding),function(x) unlist(strsplit(x,'.sum'))[1])
    StudyName <- sapply(StudyName,function(x) gsub('[.]','_',x))
    df <- data.frame(StudyName,numSigGenesCoding,numSigGenesAll)
    colnames(df) <- c('Study Name','Sig Coding Genes','Sig Coding + lncRNA Genes')
    rownames(df) <- c()
    df <- left_join(df,GWASMetadata,by=c('Study Name'='StudyName')) %>% dplyr::rename('Num Subjects'=N)
    grid.arrange(tableGrob(df))
    
  }else{
    #Genes which are include in Coding co-exp network 
    PASCAL_cluster_coding_path <- '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_results/Jaccard_Cor0p2/CodingGenesMicroglia_Jaccard_cor0p2_abs/'
    PASCAL_results_coding <- dir(PASCAL_cluster_coding_path)
    PASCAL_results_coding <- PASCAL_results_coding[grep('sum.genescores',PASCAL_results_coding)]
    numSigGenesIncludedCoding <- sapply(PASCAL_results_coding,function(x) data.table::fread(paste0(PASCAL_cluster_coding_path,x)) %>% 
                                          dplyr::filter(pvalue < 5e-8) %>% nrow)
    sigGenesCoding <- lapply(PASCAL_results_coding,function(x) data.table::fread(paste0(PASCAL_cluster_coding_path,x)) %>% 
                               dplyr::filter(pvalue < 5e-8) %>% select(gene_id) %>% t() %>% as.vector)
    
    #Genes which are included in Coding + lncRNA co-exp network.
    PASCAL_cluster_all_path <- '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_results/Jaccard_Cor0p2/AllGenesMicroglia_Jaccard_cor0p2_abs/'
    PASCAL_results_all <- dir(PASCAL_cluster_all_path)
    PASCAL_results_all <- PASCAL_results_all[grep('sum.genescores',PASCAL_results_all)]
    numSigGenesIncludedAll <- sapply(PASCAL_results_all,function(x) data.table::fread(paste0(PASCAL_cluster_all_path,x)) %>% 
                                       dplyr::filter(pvalue < 5e-8) %>% nrow)
    sigGenesAll <- lapply(PASCAL_results_all,function(x) data.table::fread(paste0(PASCAL_cluster_all_path,x)) %>% 
                            dplyr::filter(pvalue < 5e-8) %>% select(gene_id) %>% t() %>% as.vector)
    
    #Store results in DF
    SigStudiesStats <- data.frame(N_Sig_Coding_Genes_In_Network=numSigGenesIncludedCoding,
                                  N_Sig_Coding_and_LncRNA_Genes_In_Network=numSigGenesIncludedAll)
    SigStudiesStats$StudyName <- sapply(rownames(SigStudiesStats),function(x) unlist(strsplit(x,'.sum'))[1])
    SigStudiesStats$StudyName <- sapply(SigStudiesStats$StudyName,function(x) gsub('[.]','_',x))
    
    numSigClusters <- sigClusters %>% dplyr::count(StudyName)
    numSigClusters$n[is.na(numSigClusters$n)] <- 0  #Replace NA with 0
    
    #Find number of significant genes in significant cluster, for coding and coding+lncRNA co-exp network
    numSigGenesInSigClustersCoding <- rep(NA,length(sigGenesCoding))
    numSigGenesInSigClustersAll <- rep(NA,length(sigGenesCoding))
    for(i in 1:length(sigGenesCoding)){
      currentStudy <- SigStudiesStats$StudyName[i]
      
      #consider clusters from coding network
      currentStudySigClusters <- dplyr::filter(sigClusters,StudyName==currentStudy & Biotype=='coding')
      currentStudyIncludedGenes <- unique(unlist(currentStudySigClusters$Genes))
      numSigGenesInSigClustersCoding[i] <- length(intersect(sigGenesCoding[[i]],currentStudyIncludedGenes))
      
      #consider clusters from coding + lncRNA network
      currentStudySigClusters <- dplyr::filter(sigClusters,StudyName==currentStudy & Biotype=='all')
      currentStudyIncludedGenes <- unique(unlist(currentStudySigClusters$Genes))
      numSigGenesInSigClustersAll[i] <- length(intersect(sigGenesAll[[i]],currentStudyIncludedGenes))
      
    }
    SigStudiesStats$N_SigCodingGenesInSigClusters <- numSigGenesInSigClustersCoding
    SigStudiesStats$N_SigAllGenesInSigClusters <- numSigGenesInSigClustersAll
    
    GWASFullMetadata <- dplyr::left_join(SigStudiesStats,GWASMetadata,by=c('StudyName'='StudyName')) %>% 
      dplyr::left_join(numSigClusters,by=c('StudyName'='StudyName')) %>% rename_('N_Subjects'='N','N_Sig_Clusters'='n')
    GWASFullMetadata <- GWASFullMetadata[,c('StudyName','N_Sig_Coding_Genes_In_Network',
                                            'N_SigCodingGenesInSigClusters','N_Sig_Coding_and_LncRNA_Genes_In_Network',
                                            'N_SigAllGenesInSigClusters')]
    colnames(GWASFullMetadata) <- c('Study Name','Sig Coding Genes \n In Network',
                                    'Sig Coding Genes \n In Sig Cluster',
                                    'Sig Coding + lncRNA Genes \n In Network',
                                    'Sig Coding + lncRNA Genes \n In Sig Cluster')
    totalRow <- data.frame('TOTAL',sum(GWASFullMetadata[,2]),
                           sum(GWASFullMetadata[,3]),
                           sum(GWASFullMetadata[,4]),
                           sum(GWASFullMetadata[,5]))
    colnames(totalRow) <- colnames(GWASFullMetadata)
    GWASFullMetadata <- rbind(GWASFullMetadata,totalRow)
    library(gridExtra)
    t1 <- ttheme_default(core=list(
      fg_params=list(fontface=c(rep("plain", nrow(GWASFullMetadata)-1), "bold.italic")),
      bg_params = list(fill=c(rep(c("grey95", "grey90"),
                                  length.out=nrow(GWASFullMetadata)-1), "yellow"))
    ))
    grid.arrange(tableGrob(GWASFullMetadata,theme=t1))
  }
}

