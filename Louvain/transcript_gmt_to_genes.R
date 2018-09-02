#Convert modules with transcripts into genes. Remove duplicated genes. 
WriteGeneGmt <- function(GmtPath,outFile){
  #Generate transcript to gene mapping
  library(qusage)
  library(dplyr)
  source('../BatchCorrection_and_QC/load_GTF.R')
  geneGtfTableFull <- LoadGTF(full=T)
  geneGtfTableFull <- geneGtfTableFull %>% dplyr::filter(feature=='transcript')
  transcriptIdToGene <- hashmap::hashmap(keys = geneGtfTableFull$transcript_id,values = geneGtfTableFull$gene_id)
  
  #Read modules
  Gmt <- qusage::read.gmt(GmtPath)
  
  #Convert all transcripts to genes
  GmtAsGene <- vector(mode= 'list',length = length(Gmt))
  transcriptToGeneRatio <- rep(NA,length(Gmt))
  for(i in 1:length(Gmt)){
    currentTranscripts <- Gmt[[i]]
    currentGenes <- unique(transcriptIdToGene[[currentTranscripts]])
    transcriptToGeneRatio[i] <- length(currentTranscripts)/length(currentGenes)
    GmtAsGene[[i]] <- currentGenes
  }
  names(GmtAsGene) <- names(Gmt)
  
  levels <- sapply(names(GmtAsGene),function(x)as.numeric(unlist(strsplit(x,'_'))[5]))
  ratioDf <- data.frame(Level=levels,transcriptToGeneRatio=transcriptToGeneRatio,stringsAsFactors = F)
  # p <- ggplot(ratioDf, aes(x=as.factor(Level), y=transcriptToGeneRatio)) + 
  #   geom_boxplot()
  #Write gmt file with the converted genes.
  for(i in 1:length(GmtAsGene)){
    currentName <- names(GmtAsGene)[i]
    currentGenes <- GmtAsGene[[i]]
    currentLength <- length(currentGenes)
    
    currentString <- paste0(currentName,'\t',currentLength,'\t',paste(currentGenes,collapse = '\t'))
    write(currentString,file=outFile,append = T)
  }
}
GmtPath <- '/local/data/public/zmx21/zmx21_private/GSK/Louvain_results/AllMicrogliaTranscripts/AllMicrogliaTranscripts.gmt'
outFile <- file('/local/data/public/zmx21/zmx21_private/GSK/Louvain_results/AllMicrogliaTranscripts/AllMicrogliaTranscriptsAsGenes.gmt','w')
WriteGeneGmt(GmtPath,outFile)

GmtPath <- '/local/data/public/zmx21/zmx21_private/GSK/Louvain_results/Pearson_Cor0p2/AllPearsonTranscripts/all_transcripts.gmt'
outFile <- file('/local/data/public/zmx21/zmx21_private/GSK/Louvain_results/Pearson_Cor0p2/AllTranscriptsMicroglia_Pearson_cor0p2_abs.gmt','w')
WriteGeneGmt(GmtPath,outFile)

GmtPath <- '/local/data/public/zmx21/zmx21_private/GSK/Louvain_results/Pearson_Cor0p2/CodingPearsonTranscripts/coding_transcripts.gmt'
outFile <- file('/local/data/public/zmx21/zmx21_private/GSK/Louvain_results/Pearson_Cor0p2/CodingTranscriptsMicroglia_Pearson_cor0p2_abs.gmt','w')
WriteGeneGmt(GmtPath,outFile)

GmtPath <- '/local/data/public/zmx21/zmx21_private/GSK/WGCNA_clusters/Transcript_Level/CodingWGCNAUnsigned_Soft4_Size3_DeepSplit2.gmt'
outFile <- file('/local/data/public/zmx21/zmx21_private/GSK/WGCNA_clusters/Transcript_Level/Transcript_As_Genes/CodingWGCNAUnsigned_Soft4_Size3_DeepSplit2.gmt','w')
WriteGeneGmt(GmtPath,outFile)
