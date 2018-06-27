SeperateCodingandNonCoding <- function(inputMatrix){
  library(dplyr)
  source('load_GTF.R')
  bioTypeTable <- ExtractBioType(inputMatrix,table = F)
  codingIDs <- bioTypeTable == 'protein_coding'
  return(list(coding=inputMatrix[codingIDs,],noncoding=inputMatrix[!codingIDs,]))
}
CalcCV <- function(inputMatrix){
  return(apply(inputMatrix,1,function(x) sd(x)/mean(x)))
}
ApplyCVFilter <- function(inputMatrix,percentile){
  cv <- CalcCV(inputMatrix)
  return(inputMatrix[rownames(inputMatrix)[cv > quantile(cv,percentile)],])
}

#Filter genes in batch corrected (or TPM filtered) matrix, according to CV.
FilterByCV <- function(plots=F,inputPath,outputPath,percentile){
  source('load_GTF.R')
  #Load the all data specified in the inputPath list
  for(i in 1:length(inputPath)){
    invisible(lapply(paste0(inputPath[[i]],dir(inputPath[[i]])),load,environment()))
  }

  #Log transform brain data
  TPM_WholeBrain_Gene <- log2(TPM_WholeBrain_Gene+1)
  TPM_WholeBrain_Transcript <- log2(TPM_WholeBrain_Transcript+1)
  
  #Seperate Non-Coding and Protein Coding
  microgliaGene <- SeperateCodingandNonCoding(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged)
  microgliaTranscript <- SeperateCodingandNonCoding(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged)
  brainGene <- SeperateCodingandNonCoding(TPM_WholeBrain_Gene)
  brainTranscript <- SeperateCodingandNonCoding(TPM_WholeBrain_Transcript)
  
  #Apply filter for microglia, according to CV above percentile 
  MicrogliaGeneCVFiltered <- lapply(microgliaGene,function(x) ApplyCVFilter(x,percentile))
  MicrogliaTranscriptCVFiltered <- lapply(microgliaTranscript,function(x) ApplyCVFilter(x,percentile))
  
  #Apply filter for brain, according to CV above percentile 
  BrainGeneCVFiltered <- lapply(microgliaGene,function(x) ApplyCVFilter(x,percentile))
  BrainTranscriptCVFiltered <- lapply(microgliaTranscript,function(x) ApplyCVFilter(x,percentile))
  
  if(plots){
    library(egg)
    library(ggplot2)
    library(viridis)
    get_density <- function(x, y, n = 100) {
      dens <- MASS::kde2d(x = x, y = y, n = n)
      ix <- findInterval(x, dens$x)
      iy <- findInterval(y, dens$y)
      ii <- cbind(ix, iy)
      return(dens$z[ii])
    }
    #TPM and CV, gene level for microglia comparing biotype
    tiff(filename = '../Figures/CV_Filtering/Microglia_TPMVSCV_Biotype.tiff',width = 1200,height=800)
    DfGene  <- rbind(data.frame(mean = apply(microgliaGene$coding,1,function(x) mean(x)),
                            cv = apply(microgliaGene$coding,1,function(x) sd(x)/mean(x)),
                            biotype = rep('Coding',nrow(microgliaGene$coding))),
    data.frame(mean = apply(microgliaGene$noncoding,1,function(x) mean(x)), 
                            cv = apply(microgliaGene$noncoding,1,function(x) sd(x)/mean(x)),
                            biotype = rep('Non-Coding',nrow(microgliaGene$noncoding))))
    
    DfTranscript  <- rbind(data.frame(mean = apply(microgliaTranscript$coding,1,function(x) mean(x)),
                                cv = apply(microgliaTranscript$coding,1,function(x) sd(x)/mean(x)),
                                biotype = rep('Coding',nrow(microgliaTranscript$coding))),
                     data.frame(mean = apply(microgliaTranscript$noncoding,1,function(x) mean(x)), 
                                cv = apply(microgliaTranscript$noncoding,1,function(x) sd(x)/mean(x)),
                                biotype = rep('Non-Coding',nrow(microgliaTranscript$noncoding))))
    
    
    DfGeneCoding <- filter(DfGene,biotype=='Coding') %>% dplyr::select(-biotype)
    DfGeneNonCoding <- filter(DfGene,biotype=='Non-Coding')%>% dplyr::select(-biotype)
    DfGeneCoding$density <- get_density(DfGeneCoding$mean,DfGeneCoding$cv)
    DfGeneNonCoding$density <- get_density(DfGeneNonCoding$mean,DfGeneNonCoding$cv)
    DfTranscriptCoding <- filter(DfTranscript,biotype=='Coding') %>% dplyr::select(-biotype)
    DfTranscriptNonCoding <- filter(DfTranscript,biotype=='Non-Coding')%>% dplyr::select(-biotype)
    DfTranscriptCoding$density <- get_density(DfTranscriptCoding$mean,DfTranscriptCoding$cv)
    DfTranscriptNonCoding$density <- get_density(DfTranscriptNonCoding$mean,DfTranscriptNonCoding$cv)
    
    
    p1 <- ggplot(DfGeneCoding) + aes(x=mean,y = cv,colour=density) + geom_point(size=0.5) + 
      labs(x='Mean TPM',y='CV') + ggtitle('Coding Genes Microglia - mean TPM vs CV') + 
      theme_bw() + ylim(c(0,max(DfGeneNonCoding$cv))) + xlim(c(0,max(DfGeneCoding$mean)))
    
    
    p2 <- ggplot(DfGeneNonCoding) + aes(x=mean,y=cv,colour=density) + geom_point(size=0.5) +
      labs(x='Mean TPM',y='CV') + ggtitle('Noncoding Genes Microglia - mean TPM vs CV') + 
      theme_bw() + ylim(c(0,max(DfGeneNonCoding$cv))) + xlim(c(0,max(DfGeneCoding$mean)))
    
    p3 <- ggplot(DfTranscriptCoding) + aes(x=mean,y = cv,colour=density) + geom_point(size=0.5) + 
      labs(x='Mean TPM',y='CV') + ggtitle('Coding Transcript Microglia - mean TPM vs CV') + 
      theme_bw() + ylim(c(0,max(DfTranscriptNonCoding$cv))) + xlim(c(0,max(DfTranscriptCoding$mean)))
    
    p4 <- ggplot(DfTranscriptNonCoding) + aes(x=mean,y=cv,colour=density) + geom_point(size=0.5) +
      labs(x='Mean TPM',y='CV') + ggtitle('Noncoding Transcript Microglia - mean TPM vs CV') + 
      theme_bw() + ylim(c(0,max(DfTranscriptNonCoding$cv))) + xlim(c(0,max(DfTranscriptCoding$mean)))
    egg::ggarrange(p1,p2,p3,p4,nrow=2,ncol=2)
    
    dev.off()
    
    #Distribution of Microglia CV, according to Biotype
    tiff(filename = '../Figures/CV_Filtering/MicrgoliaCVDist.tiff',width = 1200,height=800)
    par(mfrow=c(2,2))
    hist(CalcCV(microgliaGene$coding),breaks=100,main='CV Gene Micrgolia - Coding',xlab='CV')
    lines(rep(quantile(CalcCV(microgliaGene$coding),percentile),2),y=c(0,100000),col='red',lwd=3)
    hist(CalcCV(microgliaTranscript$coding),breaks=100,main='CV Transcript Microglia - Coding',xlab='CV')
    lines(rep(quantile(CalcCV(microgliaTranscript$coding),percentile),2),y=c(0,100000),col='red',lwd=3)
    hist(CalcCV(microgliaGene$noncoding),breaks=100,main='CV Gene Micrgolia - Non Coding',xlab='CV')
    lines(rep(quantile(CalcCV(microgliaGene$noncoding),percentile),2),y=c(0,100000),col='red',lwd=3)
    hist(CalcCV(microgliaTranscript$noncoding),breaks=100,main='CV Transcript Microglia - Non Coding',xlab='CV')
    lines(rep(quantile(CalcCV(microgliaTranscript$noncoding),percentile),2),y=c(0,100000),col='red',lwd=3)
    dev.off()
    
    #Distribution of Brain CV, according to Biotype
    tiff(filename = '../Figures/CV_Filtering/BrainCVDist.tiff',width = 1200,height=800)
    par(mfrow=c(2,2))
    hist(CalcCV(brainGene$coding),breaks=100,main='CV Gene Micrgolia - Coding',xlab='CV')
    lines(rep(quantile(CalcCV(brainGene$coding),percentile),2),y=c(0,100000),col='red',lwd=3)
    hist(CalcCV(brainTranscript$coding),breaks=100,main='CV Transcript Microglia - Coding',xlab='CV')
    lines(rep(quantile(CalcCV(brainTranscript$coding),percentile),2),y=c(0,100000),col='red',lwd=3)
    hist(CalcCV(brainGene$noncoding),breaks=100,main='CV Gene Micrgolia - Non Coding',xlab='CV')
    lines(rep(quantile(CalcCV(brainGene$noncoding),percentile),2),y=c(0,100000),col='red',lwd=3)
    hist(CalcCV(brainTranscript$noncoding),breaks=100,main='CV Transcript Microglia - Non Coding',xlab='CV')
    lines(rep(quantile(CalcCV(brainTranscript$noncoding),percentile),2),y=c(0,100000),col='red',lwd=3)
    dev.off()
    
    #Pre and post filtering stats, for microglia
    Filtering_Biotype_Microglia <- lapply(list(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged,
                                               do.call(rbind,MicrogliaGeneCVFiltered),
                                               SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged,
                                               do.call(rbind,MicrogliaTranscriptCVFiltered)),function(x) ExtractBioType(x)) %>% {do.call(rbind,.)} %>% as.data.frame()
    Filtering_Biotype_Microglia$Total <- rowSums(Filtering_Biotype_Microglia)
    Filtering_Biotype_Microglia$Step <- c('Gene - Pre Filt','Gene - Post Filt','Transcript - Pre Filt','Transcript Post Filt')
    Filtering_Biotype_Microglia <- Filtering_Biotype_Microglia[,c(6,1,2,3,4,5)]  #Switch column orders
    
    tiff(height=200,width=500,filename = '../Figures/CV_Filtering/MicrgoliaCVFilteringStats.tiff')
    grid.arrange(top='Microglia CV Filtering',tableGrob(Filtering_Biotype_Microglia))
    dev.off()
    
    #Pre and post filtering stats, for brain
    Filtering_Biotype_Brain <- lapply(list(TPM_WholeBrain_Gene,
                                           do.call(rbind,BrainGeneCVFiltered),
                                           TPM_WholeBrain_Transcript,
                                           do.call(rbind,BrainTranscriptCVFiltered)),function(x) ExtractBioType(x)) %>% {do.call(rbind,.)} %>% as.data.frame()
    Filtering_Biotype_Brain$Total <- rowSums(Filtering_Biotype_Brain)
    Filtering_Biotype_Brain$Step <- c('Gene - Pre Filt','Gene - Post Filt','Transcript - Pre Filt','Transcript Post Filt')
    Filtering_Biotype_Brain <- Filtering_Biotype_Brain[,c(6,1,2,3,4,5)]  #Switch column orders
    
    tiff(height=200,width=500,filename = '../Figures/CV_Filtering/BrainCVFilteringStats.tiff')
    grid.arrange(top='Brain CV Filtering',tableGrob(Filtering_Biotype_Brain))
    dev.off()
    
  }
  save(MicrogliaGeneCVFiltered,file='../Count_Data/CV_Filtered/MicrogliaGeneCVFiltered.rda')
  save(MicrogliaTranscriptCVFiltered,file='../Count_Data/CV_Filtered/MicrogliaTranscriptCVFiltered.rda')
  save(BrainGeneCVFiltered,file='../Count_Data/CV_Filtered/BrainGeneCVFiltered.rda')
  save(BrainTranscriptCVFiltered,file = '../Count_Data/CV_Filtered/BrainTranscriptCVFiltered.rda')
  save(list = ls(environment()),file='../CodeImages/CVFiltering.RData')
}
FilterByCV(plots=T,list('../Count_Data/Batch_Corrected/','../Count_Data/TPM_Filtered/'),outputPath = '../Count_Data/CV_Filtered/',0.5)
  