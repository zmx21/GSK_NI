library(qusage)
library(org.Hs.eg.db)
library(dplyr)
biocartaGMT <-qusage::read.gmt(file='/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_Ensembl/resources/genesets/msigdb/msigBIOCARTA_KEGG_REACTOME.gmt')
biocartaGMT <- lapply(1:length(biocartaGMT),function(i) list(Name = names(biocartaGMT)[i],genes = biocartaGMT[[i]],GO = "None"))

#Get Entrez to Ensembl Mapping
db <- org.Hs.egENSEMBL
mapped_genes <- mappedkeys(db)
mapping <- as.data.frame(db[mapped_genes])

#Writes gene modules into a GMT file acceptable by PASCAL
ConstructGMT <- function(geneModules){
  gmtDf <- data.frame(ClusterName = character(length = length(geneModules)),
                      ClusterDescript = character(length = length(geneModules)),
                      Genes = character(length=length(geneModules)),stringsAsFactors = F)
  for(i in 1:length(geneModules)){
    gmtDf$ClusterName[i] <- geneModules[[i]]$Name
    gmtDf$ClusterDescript[i] <- geneModules[[i]]$GO
    gmtDf$Genes[i] <- paste(geneModules[[i]]$genes,collapse = '\t')
  }
  return(gmtDf)
}

biocartaGMTEnsembl <- mclapply(biocartaGMT,function(x) {
  ENSGenes <- dplyr::left_join(data.frame(gene_id = x$genes,stringsAsFactors = F),mapping,by="gene_id") %>% 
    dplyr::filter(!is.na(ensembl_id)) %>% dplyr::select('ensembl_id') %>% t() %>% as.vector()
  return(list(Name=x$Name,GO=x$GO,genes=ENSGenes))
},mc.cores = 60)
biocartaGMTEnsembl <- ConstructGMT(biocartaGMTEnsembl)
write.table(biocartaGMTEnsembl,row.names = F,col.names = F,quote = F,sep = '\t',file='/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_Ensembl/resources/genesets/msigdb/msigBIOCARTA_KEGG_REACTOME_Ensembl.gmt')
