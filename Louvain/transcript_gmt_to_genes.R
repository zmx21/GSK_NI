WriteGeneGmt <- function(GmtPath,outFile){
  library(qusage)
  library(dplyr)
  source('../BatchCorrection_and_QC/load_GTF.R')
  geneGtfTableFull <- LoadGTF(full=T)
  geneGtfTableFull <- geneGtfTableFull %>% dplyr::filter(feature=='transcript')
  transcriptIdToGene <- hashmap::hashmap(keys = geneGtfTableFull$transcript_id,values = geneGtfTableFull$gene_id)
  
  Gmt <- qusage::read.gmt(GmtPath)
  
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



#Same gene vs different gene transcripts. 
# geneGtfTableFullIncluded <- geneGtfTableFull %>% dplyr::filter(transcript_id %in% unique(unlist(codingMicrogliaTranscriptGmt)))
# load('../../Count_Data/CV_Filtered/MicrogliaTranscriptCVFiltered.rda')
# transcriptCountMatrix <- MicrogliaTranscriptCVFiltered$coding
# allGenes <- unique(geneGtfTableFullIncluded$gene_id)
# allGenesSamples <- sample(allGenes,size = 1000)
# 
# withinGroupCorr <- c()
# numTranscript <- 0
# for(i in 1:length(allGenesSamples)){
#   currentTranscripts <- dplyr::filter(geneGtfTableFullIncluded,gene_id==allGenesSamples[i]) %>% select(transcript_id) %>% t() %>% as.vector()
#   if(length(currentTranscripts)<2){
#     next
#   }
#   pairWiseCorr <- cor(t(transcriptCountMatrix[currentTranscripts,]))
#   withinGroupCorr <- c(withinGroupCorr,as.vector(pairWiseCorr[upper.tri(pairWiseCorr,diag = F)]))
#   numTranscript <- numTranscript + length(currentTranscripts)
# }
# 
# randomMatrix <- transcriptCountMatrix[sample(rownames(transcriptCountMatrix),size=numTranscript),]
# interGroupCorrMat <- cor(t(randomMatrix))
# interGroupcorr <- as.vector(interGroupCorrMat[upper.tri(interGroupCorrMat,diag = F)])
# 
# 
# plot(density(abs(sample(interGroupcorr,size=length(withinGroupCorr)))),col='red',main='Same Gene Transcripts VS Random Transcripts',xlab='abs(cor)')
# lines(density(abs(withinGroupCorr)),col='blue')
# legend(x='topright',col=c('red','blue'),legend = c('Random','Same Gene'),lty=c(1,1))
