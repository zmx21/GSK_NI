##########################################################################################################
#Create edgelist, genelist, and heatlist for HotNet2 
##########################################################################################################

library(data.table)
library(dplyr)
library(parallel)

CreateEdgeList <- function(PPIDb){
  #Obtain HGNC to Entrez mapping
  EntrezHGNC <- data.table::fread(file='/local/data/public/zmx21/zmx21_private/GSK/GWAS/HINT/ENTREZ_HGNC_ENSEMBL') %>% 
    select(HGNC='HGNC ID',Entrez='Entrez Gene ID',EntrezNCBI = 'Entrez Gene ID(supplied by NCBI)',Name='Approved Symbol') %>% 
    filter(!(is.na(EntrezNCBI) & is.na(Entrez))) %>% mutate(EntrezFull = if_else(is.na(Entrez),EntrezNCBI,Entrez)) %>% select(-Entrez,-EntrezNCBI) %>% rename(Entrez = EntrezFull)
  #Convert PPI to Entrex ID
  PPIDb <- filter(PPIDb,Alias_A!="" & Alias_B!="") %>% select(GeneA=Alias_A,GeneB=Alias_B)
  PPIGeneAEntrez <- mclapply(PPIDb$GeneA,function(x) dplyr::left_join(data.frame(HGNC=unlist(strsplit(x,split = '\\|')),stringsAsFactors = F),EntrezHGNC,by=c('HGNC'='HGNC')) %>% 
                               select(Entrez) %>% filter(!is.na(Entrez)) %>% {unique(.$Entrez)},mc.cores = 60)
  PPIGeneBEntrez <- mclapply(PPIDb$GeneB,function(x) dplyr::left_join(data.frame(HGNC=unlist(strsplit(x,split = '\\|')),stringsAsFactors = F),EntrezHGNC,by=c('HGNC'='HGNC')) %>% 
                               select(Entrez) %>% filter(!is.na(Entrez)) %>% {unique(.$Entrez)},mc.cores = 60)
  #Only include Genes which have mapping.
  GeneAIncl <- which(sapply(PPIGeneAEntrez,function(x) !is.null(x)))
  GeneBIncl <-  which(sapply(PPIGeneBEntrez,function(x) !is.null(x)))
  PairsIncl <- intersect(GeneAIncl,GeneBIncl)
  GeneA <-  PPIGeneAEntrez[PairsIncl]; GeneB = PPIGeneBEntrez[PairsIncl]
  
  #"Expand" the list, by constructing edge list of all pairs.
  edgeList <- data.frame(GeneA = numeric(),GeneB = numeric(),stringsAsFactors = F)
  for(i in 1:length(PairsIncl)){
    currentGeneA <- GeneA[[i]]
    currentGeneB <- GeneB[[i]]
    for(j in 1:length(currentGeneA)){
      for(k in 1:length(currentGeneB)){
        edgeList <- tryCatch({rbind(edgeList,data.frame(GeneA=currentGeneA[j],GeneB=currentGeneB[k],stringsAsFactors = F))},error=function(cond){edgeList})
      }
    }
  }
  edgeList <- dplyr::filter(edgeList,GeneA != GeneB)
  return(edgeList)
}
WriteHotnetFiles <- function(EntrezEdgeList,GWASDb,path){
  system(command = paste0('mkdir ',path)) #Make directory to write filtes
  
  #Get all unique genes
  allgenes <- union(unique(EntrezEdgeList$GeneA),unique(EntrezEdgeList$GeneB))
  #Get index of each gene of the Entrez edge list, based on the order of allgenes.
  GeneAIndex <- sapply(EntrezEdgeList$GeneA,function(x) which(x==allgenes))
  GeneBIndex <- sapply(EntrezEdgeList$GeneB,function(x) which(x==allgenes))
  #Edge list based on index
  indexEdgeList <- data.frame(GeneA=GeneAIndex,GeneB=GeneBIndex)
  #Mapping between index and gene name (Entrez ID)
  geneList <- data.frame(Index = 1:length(allgenes),Gene = allgenes)
  #Mapping between each gene, and it's head score (-log(pval))
  heatList <- dplyr::left_join(geneList,GWASDb,by=c('Gene'='gene_id')) %>% dplyr::mutate(Heat = -1*log(pvalue), Gene = paste0('ENTREZ:',Gene)) %>% dplyr::select(Gene,Heat)
  geneList <- mutate(geneList,Gene=paste0('ENTREZ:',as.character(Gene)))
  #Write files to the directory.
  write.table(indexEdgeList,quote = F,row.names = F,col.names = F,file = paste0(path,'edgelist.txt'))
  write.table(geneList,quote = F,row.names = F,col.names = F,file = paste0(path,'genelist.txt'))
  write.table(heatList,quote = F,row.names = F,col.names = F,file = paste0(path,'heatlist.txt'))
}

#Create network files based on PPI, for all datasets. 
PPIDb <- data.table::fread(file='/local/data/public/zmx21/zmx21_private/GSK/GWAS/HINT/HomoSapiens_lcb_hq.txt')
GWASAlzheimer <- data.table::fread(file='/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_results/KEGG/Alzheimer/IGAP_stage_1_parsed.sum.genescores.txt')
PPIEdgeList <- CreateEdgeList(PPIDb)
AlzheimerEdgeList <- filter(PPIEdgeList,GeneA%in% GWASAlzheimer$gene_id & GeneB %in% GWASAlzheimer$gene_id)  #Genes must be in GWAS
WriteHotnetFiles(AlzheimerEdgeList,GWASAlzheimer,'/local/data/public/zmx21/zmx21_private/GSK/GWAS/HotNet_PPI/Alzheimer_Stage1/')

GWASAlzheimer2 <- data.table::fread(file='/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_results/KEGG/Alzheimer/IGAP_stage_1_2_combined.parsed.sum.genescores.txt')
AlzheimerEdgeList2 <- filter(PPIEdgeList,GeneA%in% GWASAlzheimer2$gene_id & GeneB %in% GWASAlzheimer2$gene_id)  #Genes must be in GWAS
WriteHotnetFiles(AlzheimerEdgeList2,GWASAlzheimer2,'/local/data/public/zmx21/zmx21_private/GSK/GWAS/HotNet_PPI/Alzheimer_Stage2/')

GWASALS <- data.table::fread(file='/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_results/KEGG/ALS/als.meta.all.parsed.sum.genescores.txt')
ALSEdgeList <- filter(PPIEdgeList,GeneA%in% GWASALS$gene_id & GeneB %in% GWASALS$gene_id)  #Genes must be in GWAS
WriteHotnetFiles(ALSEdgeList,GWASALS,'/local/data/public/zmx21/zmx21_private/GSK/GWAS/HotNet_PPI/ALS/')

GWASParkinsons <- data.table::fread(file='/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_results/KEGG/Parkinsons/Parkinsons.parsed.sum.genescores.txt')
ParkinsonEdgeList <- filter(PPIEdgeList,GeneA%in% GWASParkinsons$gene_id & GeneB %in% GWASParkinsons$gene_id)  #Genes must be in GWAS
WriteHotnetFiles(ParkinsonEdgeList,GWASParkinsons,'/local/data/public/zmx21/zmx21_private/GSK/GWAS/HotNet_PPI/Parkinsons/')

GWASMS <- data.table::fread(file='/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_results/KEGG/MS/MS.GWAS.parsed.sum.genescores.txt')
MSEdgeList <- filter(PPIEdgeList,GeneA%in% GWASMS$gene_id & GeneB %in% GWASMS$gene_id)  #Genes must be in GWAS
WriteHotnetFiles(MSEdgeList,GWASMS,'/local/data/public/zmx21/zmx21_private/GSK/GWAS/HotNet_PPI/MS/')

GWASMSImputed <- data.table::fread(file='/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_results/KEGG/MS/MS.imputed.parsed.sum.genescores.txt')
MSImputedEdgeList <- filter(PPIEdgeList,GeneA%in% GWASMSImputed$gene_id & GeneB %in% GWASMSImputed$gene_id)  #Genes must be in GWAS
WriteHotnetFiles(MSImputedEdgeList,GWASMSImputed,'/local/data/public/zmx21/zmx21_private/GSK/GWAS/HotNet_PPI/MS_Imputed/')

