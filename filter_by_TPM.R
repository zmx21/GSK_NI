#Helper functin to extract a dataset from a merged matrix.
ExtractDataset <- function(countMatrix,Sample){
  if(Sample=='Galatro'){
    countMatrix[,sapply(colnames(countMatrix),function(x) substr(x,1,3) == 'GSM')]
  }else if (Sample=='Gosselin'){
    countMatrix[,sapply(colnames(countMatrix),function(x) substr(x,1,3) == 'SRR')]
  }else if (Sample=='Olah'){
    countMatrix[,sapply(colnames(countMatrix),function(x) substr(x,1,3) == 'H5K')]
  }
}

FilterByTPM <- function(plots=T,inputPath,outputPath,percentile){
  library(egg)
  source('load_GTF.R')
  invisible(lapply(paste0(inputPath,dir(inputPath)),load,environment()))
  
  ###################################################MICROGLIA DATA################################################################3
  
  #Filter according to expression of genes
  #The removal criteria is that of the samples are unexpressed
  #(expression level below the specified percentile of the respective sample's distribution).
  sampleQuantilesGenes <- apply(TPM_Microglia_Gene_Merged,2,function(x) quantile(x,percentile))
  sampleQuantilesTranscripts <- apply(TPM_Microglia_Transcript_Merged,2,function(x) quantile(x,percentile))
  genesToKeepMicroglia <- sapply(1:nrow(TPM_Microglia_Gene_Merged),function(i) !all(TPM_Microglia_Gene_Merged[i,] <= sampleQuantilesGenes))
  transcriptsToKeepMicroglia <- sapply(1:nrow(TPM_Microglia_Transcript_Merged),function(i) !all(TPM_Microglia_Transcript_Merged[i,] <= sampleQuantilesTranscripts))
  
  #Store prefilertered data. 
  GalatroSamples_Genes_Prefilt <- ExtractDataset(TPM_Microglia_Gene_Merged,'Galatro')
  GosselinSamples_Genes_Prefilt <- ExtractDataset(TPM_Microglia_Gene_Merged,'Gosselin')
  OlahSamples_Genes_Prefilt <- ExtractDataset(TPM_Microglia_Gene_Merged,'Olah')
  GalatroSamples_Transcripts_Prefilt <- ExtractDataset(TPM_Microglia_Transcript_Merged,'Galatro')
  GosselinSamples_Transcripts_Prefilt <- ExtractDataset(TPM_Microglia_Transcript_Merged,'Gosselin')
  OlahSamples_Transcripts_Prefilt <- ExtractDataset(TPM_Microglia_Transcript_Merged,'Olah')
  TPM_Microglia_Gene_Merged_Prefilt <- TPM_Microglia_Gene_Merged
  TPM_Microglia_Transcript_Merged_Prefilt  <- TPM_Microglia_Transcript_Merged
  
  #Filter matrix according to TPM quantile.
  TPM_Microglia_Gene_Merged <- TPM_Microglia_Gene_Merged[rownames(TPM_Microglia_Gene_Merged)[genesToKeepMicroglia],]
  TPM_Microglia_Transcript_Merged  <- TPM_Microglia_Transcript_Merged[rownames(TPM_Microglia_Transcript_Merged)[transcriptsToKeepMicroglia],]
  
  GalatroSamples_Genes <- ExtractDataset(TPM_Microglia_Gene_Merged,'Galatro')
  GosselinSamples_Genes <- ExtractDataset(TPM_Microglia_Gene_Merged,'Gosselin')
  OlahSamples_Genes <- ExtractDataset(TPM_Microglia_Gene_Merged,'Olah')
  GalatroSamples_Transcripts <- ExtractDataset(TPM_Microglia_Transcript_Merged,'Galatro')
  GosselinSamples_Transcripts <- ExtractDataset(TPM_Microglia_Transcript_Merged,'Gosselin')
  OlahSamples_Transcripts <- ExtractDataset(TPM_Microglia_Transcript_Merged,'Olah')
  

  Filtering_Biotype_Microglia <- lapply(list(TPM_Microglia_Gene_Merged_Prefilt,
                                             TPM_Microglia_Gene_Merged,
                                             TPM_Microglia_Transcript_Merged_Prefilt,
                                             TPM_Microglia_Transcript_Merged),function(x) ExtractBioType(x)) %>% {do.call(rbind,.)} %>% as.data.frame()
  Filtering_Biotype_Microglia$Total <- rowSums(Filtering_Biotype_Microglia)
  Filtering_Biotype_Microglia$Step <- c('Gene - Pre Filt','Gene - Post Filt','Transcript - Pre Filt','Transcript Post Filt')
  Filtering_Biotype_Microglia <- Filtering_Biotype_Microglia[,c(6,1,2,3,4,5)]  #Switch column orders
  
  if(plots){
    tiff(height=600,width=500,filename = '../Figures/TPM_Filtering//MicrgoliaTPMFilteringStats.tiff')
    grid.arrange(top='Microglia TPM Filtering',tableGrob(Filtering_Biotype_Microglia))
    dev.off()
  }

  #TPM Comparison Pre and Post filtering
  PrePostFilt_Gene <- rbind(data.frame(Counts=as.vector(GalatroSamples_Genes_Prefilt)) %>% mutate(Dataset = 'Galatro',Filt='Prefilt'),
                            data.frame(Counts=as.vector(GalatroSamples_Genes)) %>% mutate(Dataset = 'Galatro',Filt='Postfilt'),
                            data.frame(Counts=as.vector(GosselinSamples_Genes_Prefilt)) %>% mutate(Dataset = 'Gosselin',Filt='Prefilt'),
                            data.frame(Counts=as.vector(GosselinSamples_Genes)) %>% mutate(Dataset = 'Gosselin',Filt='Postfilt'),
                            data.frame(Counts=as.vector(OlahSamples_Genes_Prefilt)) %>% mutate(Dataset = 'Olah',Filt='Prefilt'),
                            data.frame(Counts=as.vector(OlahSamples_Genes)) %>% mutate(Dataset = 'Olah',Filt='Postfilt'))
  PrePostFilt_Transcript <- rbind(data.frame(Counts=as.vector(GalatroSamples_Transcripts_Prefilt)) %>% mutate(Dataset = 'Galatro',Filt='Prefilt'),
                                  data.frame(Counts=as.vector(GalatroSamples_Transcripts)) %>% mutate(Dataset = 'Galatro',Filt='Postfilt'),
                                  data.frame(Counts=as.vector(GosselinSamples_Transcripts_Prefilt)) %>% mutate(Dataset = 'Gosselin',Filt='Prefilt'),
                                  data.frame(Counts=as.vector(GosselinSamples_Transcripts)) %>% mutate(Dataset = 'Gosselin',Filt='Postfilt'),
                                  data.frame(Counts=as.vector(OlahSamples_Transcripts_Prefilt)) %>% mutate(Dataset = 'Olah',Filt='Prefilt'),
                                  data.frame(Counts=as.vector(OlahSamples_Transcripts)) %>% mutate(Dataset = 'Olah',Filt='Postfilt'))
  if(plots){
      tiff(height=600,width=1200,filename = '../Figures/TPM_Comparison/PrePostTPMFiltering_MicrogliaTPM.tiff')
      ggarrange(ggplot(PrePostFilt_Gene, aes(x=Dataset, y=Counts, fill=as.factor(Filt))) + 
                geom_boxplot(outlier.shape = NA) + 
                ggtitle('Gene Level') + 
                scale_x_discrete(name = "Dataset") +
                scale_y_continuous(name = "TPM",limits = c(0,10)),ggplot(PrePostFilt_Transcript, aes(x=Dataset, y=Counts, fill=as.factor(Filt))) + 
                geom_boxplot(outlier.shape = NA) + 
                scale_x_discrete(name = "Dataset") + 
                ggtitle('Transcript Level')+
                scale_y_continuous(name = "TPM",limits = c(0,10)),nrow=2)
      dev.off()
  }
  ###################################################BRAIN DATA################################################################3
  #Filter by TPM quantile, same criteria as Microglia data
  sampleQuantilesGenesBrain <- apply(TPM_WholeBrain_Gene,2,function(x) quantile(x,percentile))
  sampleQuantilesTranscriptsBrain <- apply(TPM_WholeBrain_Transcript,2,function(x) quantile(x,percentile))
  genesToKeepBrain <- sapply(1:nrow(TPM_WholeBrain_Gene),function(i) !all(TPM_WholeBrain_Gene[i,] <= sampleQuantilesGenesBrain))
  transcriptsToKeepBrain <- sapply(1:nrow(TPM_WholeBrain_Transcript),function(i) !all(TPM_WholeBrain_Transcript[i,] <= sampleQuantilesTranscriptsBrain))
  
  #Store Prefiltered Data 
  TPM_WholeBrain_Gene_Prefilt <- TPM_WholeBrain_Gene
  TPM_WholeBrain_Transcript_Prefilt <- TPM_WholeBrain_Transcript
  
  #Remove Filtered Genes
  TPM_WholeBrain_Gene <- TPM_WholeBrain_Gene[rownames(TPM_WholeBrain_Gene)[genesToKeepBrain],]
  TPM_WholeBrain_Transcript <- TPM_WholeBrain_Transcript[rownames(TPM_WholeBrain_Transcript)[transcriptsToKeepBrain],]
  
  #Biotype Comparison Pre and Post filtering
  Filtering_Biotype_WholeBrain <- lapply(list(TPM_WholeBrain_Gene_Prefilt,
                                              TPM_WholeBrain_Gene,
                                              TPM_WholeBrain_Transcript_Prefilt,
                                              TPM_WholeBrain_Transcript),function(x) ExtractBioType(x)) %>% {do.call(rbind,.)} %>% as.data.frame()
  Filtering_Biotype_WholeBrain$Total <- rowSums(Filtering_Biotype_WholeBrain)
  Filtering_Biotype_WholeBrain$Step <- c('Gene - Pre Filt','Gene - Post Filt','Transcript - Pre Filt','Transcript Post Filt')
  Filtering_Biotype_WholeBrain <- Filtering_Biotype_WholeBrain[,c(6,1,2,3,4,5)]  #Switch column orders
  
  if(plots){
    tiff(height=600,width=500,filename = '../Figures/TPM_Filtering/WholeBrainTPMFilteringStats.tiff')
    grid.arrange(top='Whole Brain TPM Filtering',tableGrob(Filtering_Biotype_WholeBrain))
    dev.off()
  }
  
  
  #TPM Comparison Pre and Post filtering
  PrePostFilt_Gene_Brain <- rbind(data.frame(Counts=as.vector(TPM_WholeBrain_Gene_Prefilt)) %>% mutate(Dataset = 'Brain',Filt='Prefilt'),
                                  data.frame(Counts=as.vector(TPM_WholeBrain_Gene)) %>% mutate(Dataset = 'Brain',Filt='Postfilt'))
  PrePostFilt_Transcript_Brain <- rbind(data.frame(Counts=as.vector(TPM_WholeBrain_Transcript_Prefilt)) %>% mutate(Dataset = 'Brain',Filt='Prefilt'),
                                        data.frame(Counts=as.vector(TPM_WholeBrain_Transcript)) %>% mutate(Dataset = 'Brain',Filt='Postfilt'))
  if(plots){
      tiff(height=600,width=1200,filename = '../Figures/TPM_Comparison/PrePostTPMFiltering_WholeBrainTPM.tiff')
      ggarrange(ggplot(PrePostFilt_Gene_Brain, aes(x=Dataset, y=Counts, fill=as.factor(Filt))) + 
                geom_boxplot(outlier.shape = NA) + 
                scale_x_discrete(name = "Dataset") +
                ggtitle('Gene Level')+
                scale_y_continuous(name = "TPM",limits = c(0,10)),ggplot(PrePostFilt_Transcript_Brain, aes(x=Dataset, y=Counts, fill=as.factor(Filt))) + 
                geom_boxplot(outlier.shape = NA) + 
                scale_x_discrete(name = "Dataset") +
                ggtitle('Transcript Level')+
                scale_y_continuous(name = "TPM",limits = c(0,10)),nrow = 2)
      dev.off()
  }
  save(TPM_Microglia_Gene_Merged,file=paste0(outputPath,'TPM_Microglia_Gene_Merged.rda'))
  save(TPM_Microglia_Transcript_Merged,file=paste0(outputPath,'TPM_Microglia_Transcript_Merged.rda'))
  save(TPM_WholeBrain_Gene,file=paste0(outputPath,'/TPM_WholeBrain_Gene.rda'))
  save(TPM_WholeBrain_Transcript,file= paste0(outputPath,'/TPM_WholeBrain_Transcript.rda'))

  #save(list = ls(environment()),file=paste0('../CodeImages/TPMFiltering_',as.character(percentile),'.RData'))
}
#FilterByTPM(plots = F,inputPath = '../Count_Data/Read_Filtered/',outputPath = '../Count_Data/TPM_Filtered/',percentile=0.75)













# #Helper functino to extract biotype. 
# FilterCountMatrix <- function(countMatrixInput,meanAbundances,cv,cvCutOffAbsolute,meanCutOffAbsolute,mappingTable=NULL){
#   #Filter for genes according to cutoff
#   belowMeanCutOff <- which(meanAbundances <= meanCutOffAbsolute)
#   belowCVCutoff <- which(cv <= cvCutOffAbsolute | is.nan(cv))
#   rowsToFilter <- union(belowMeanCutOff,belowCVCutoff)
#   
#   filteredCountMatrix <- countMatrixInput[setdiff(1:nrow(countMatrixInput),rowsToFilter),]
#   
#   #Decide whether to return categorial counts.
#   if(is.null(mappingTable)){
#     return(list(countMatrix = filteredCountMatrix,numFiltered = length(rowsToFilter)))
#   }else{
#     type <- ifelse(grepl('ENST',rownames(countMatrixInput)[1]),'transcript_id','gene_id')
#     allFilteredNames <- data.frame(ids = rownames(filteredCountMatrix),stringsAsFactors = F)
#     colnames(allFilteredNames) <- type
#     categoricalCounts <- as.data.frame(dplyr::left_join(allFilteredNames,mappingTable,by=c(type,type)))
#     return(list(countMatrix = filteredCountMatrix,numFiltered = length(rowsToFilter),categoricalCounts = categoricalCounts))
#   }
# }
