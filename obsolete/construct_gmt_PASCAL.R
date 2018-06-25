##########################################################################################################
#Converts gene modules encoded as a list with ENSG IDs into a gmt file with Entrez ID to be used by PASCAL
##########################################################################################################

library(data.table)
library(dplyr)
source('load_GTF.R')
#Loads data.frame which maps Ensembl Gene ID to Entrez ID, source is file provided by PASCAL
LoadEnsemblEntrezMapping <- function(){
  ensemblEntrezMapping <- data.table::fread(file='/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL/resources/annotation/ucsc/ucsc_known_genes_2013-09-03.txt',
                                            stringsAsFactors = F) %>%
    dplyr::select(EnsemblTranscript = hg19.knownToEnsembl.value,
                  GeneSymbol=hg19.kgXref.geneSymbol,
                  EntrezID=hg19.knownToLocusLink.value)
  ensemblEntrezMapping <- LoadGTF(full = F) %>% dplyr::select(EnsemblTranscript=transcript_id,EnsemblGene=gene_id) %>% 
  {dplyr::left_join(.,ensemblEntrezMapping,by='EnsemblTranscript')} %>% 
    dplyr::filter(!is.na(EntrezID),EntrezID!='n/a') %>%
    dplyr::distinct(EnsemblGene,.keep_all = T) %>%
    dplyr::select(-EnsemblTranscript)
  return(ensemblEntrezMapping)
}

#Writes gene modules into a GMT file acceptable by PASCAL
ConstructGMT <- function(geneModules){
  gmtDf <- data.frame(ClusterName = character(length = length(geneModules)),
                      ClusterDescript = character(length = length(geneModules)),
                      Genes = character(length=length(geneModules)),stringsAsFactors = F)
  for(i in 1:length(geneModules)){
    gmtDf$ClusterName[i] <- as.character(i)
    gmtDf$ClusterDescript[i] <- geneModules[[i]]$GO
    gmtDf$Genes[i] <- paste(testModules[[i]]$genes,collapse = '\t')
  }
  return(gmtDf)
}

#Map modules encoded by Ensembl ID into EntrezID, with same structure. 
EnsemblToEntrez <- function(modules,mapping){
  for(i in 1:length(modules)){
    modules[[i]]$genes <- dplyr::left_join(data.frame(EnsemblGene = modules[[i]]$genes,stringsAsFactors = F),mapping,by='EnsemblGene') %>% 
      dplyr::filter(!is.na(EntrezID)) %>% dplyr::select(EntrezID) %>% t() %>% as.character()
  }
  return(modules)
}
load('../Count_Data/testModules.rda')
# ensemblEntrezMapping <- LoadEnsemblEntrezMapping()

testModules <- EnsemblToEntrez(testModules,ensemblEntrezMapping)

gmtDf <- ConstructGMT(testModules)
write.table(gmtDf,sep = '\t',
            file='/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL/resources/genesets/microglia_modules/modules.gmt',
            quote = F,row.names = F,col.names = F)