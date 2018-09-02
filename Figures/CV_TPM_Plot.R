##########################################################################################################
#Fig 3.1.0.1
##########################################################################################################

#Extract matrix of coding and non-coding genes within input count matrix.
SeperateCodingandNonCoding <- function(inputMatrix,lncRNAOnly=T){
  library(dplyr)
  source('../BatchCorrection_and_QC/load_GTF.R')
  bioTypeTable <- ExtractBioType(inputMatrix,table = F)
  codingIDs <- bioTypeTable == 'protein_coding'
  if(lncRNAOnly){
    nonCodingIDs <- bioTypeTable == 'lncRNA'
    return(list(coding=inputMatrix[codingIDs,],noncoding=inputMatrix[nonCodingIDs,]))
  }else{
    return(list(coding=inputMatrix[codingIDs,],noncoding=inputMatrix[!codingIDs,]))
  }
}
#CV of each gene across samples
CalcCV <- function(inputMatrix){
  return(apply(inputMatrix,1,function(x) sd(x)/mean(x)))
}
#Remove genes according to CV cutoff.
ApplyCVFilter <- function(inputMatrix,percentile){
  cv <- CalcCV(inputMatrix)
  return(inputMatrix[rownames(inputMatrix)[cv > quantile(cv,percentile)],])
}

#Filter genes in batch corrected (or TPM filtered) matrix, according to CV.
PlotTPMCV <- function(inputPath,percentile){
  # source('load_GTF.R')
  #Load the all data specified in the inputPath list
  for(i in 1:length(inputPath)){
    invisible(lapply(paste0(inputPath[[i]],dir(inputPath[[i]])),load,environment()))
  }

  #Seperate Non-Coding and Protein Coding
  microgliaGene <- SeperateCodingandNonCoding(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged)
  microgliaTranscript <- SeperateCodingandNonCoding(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged)
  
  #Apply filter for microglia, according to CV above percentile 
  MicrogliaGeneCVFiltered <- lapply(microgliaGene,function(x) ApplyCVFilter(x,percentile))
  MicrogliaTranscriptCVFiltered <- lapply(microgliaTranscript,function(x) ApplyCVFilter(x,percentile))
  
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
    DfGene  <- rbind(data.frame(median = apply(microgliaGene$coding,1,function(x) median(x)),
                            cv = apply(microgliaGene$coding,1,function(x) sd(x)/mean(x)),
                            biotype = rep('Coding',nrow(microgliaGene$coding))),
    data.frame(median = apply(microgliaGene$noncoding,1,function(x) median(x)), 
                            cv = apply(microgliaGene$noncoding,1,function(x) sd(x)/mean(x)),
                            biotype = rep('Non-Coding',nrow(microgliaGene$noncoding))))
    
    DfTranscript  <- rbind(data.frame(median = apply(microgliaTranscript$coding,1,function(x) median(x)),
                                cv = apply(microgliaTranscript$coding,1,function(x) sd(x)/mean(x)),
                                biotype = rep('Coding',nrow(microgliaTranscript$coding))),
                     data.frame(median = apply(microgliaTranscript$noncoding,1,function(x) median(x)), 
                                cv = apply(microgliaTranscript$noncoding,1,function(x) sd(x)/mean(x)),
                                biotype = rep('Non-Coding',nrow(microgliaTranscript$noncoding))))
    
    
    DfGeneCoding <- filter(DfGene,biotype=='Coding') %>% dplyr::select(-biotype)
    DfGeneCoding <- DfGeneCoding[sample(1:nrow(DfGeneCoding),size = floor(1 * nrow(DfGeneCoding)),replace = F),]
    
    DfGeneNonCoding <- filter(DfGene,biotype=='Non-Coding')%>% dplyr::select(-biotype)
    DfGeneNonCoding <- DfGeneNonCoding[sample(1:nrow(DfGeneNonCoding),size = floor(1 * nrow(DfGeneNonCoding)),replace = F),]
    
    DfGeneCoding$density <- get_density(DfGeneCoding$median,DfGeneCoding$cv)
    DfGeneNonCoding$density <- get_density(DfGeneNonCoding$median,DfGeneNonCoding$cv)
    
    DfTranscriptCoding <- filter(DfTranscript,biotype=='Coding') %>% dplyr::select(-biotype)
    DfTranscriptCoding <- DfTranscriptCoding[sample(1:nrow(DfTranscriptCoding),size = floor(1 * nrow(DfTranscriptCoding)),replace = F),]
    
    DfTranscriptNonCoding <- filter(DfTranscript,biotype=='Non-Coding')%>% dplyr::select(-biotype)
    DfTranscriptNonCoding <- DfTranscriptNonCoding[sample(1:nrow(DfTranscriptNonCoding),size = floor(1 * nrow(DfTranscriptNonCoding)),replace = F),]
    
    
    DfTranscriptCoding$density <- get_density(DfTranscriptCoding$median,DfTranscriptCoding$cv)
    DfTranscriptNonCoding$density <- get_density(DfTranscriptNonCoding$median,DfTranscriptNonCoding$cv)
    
    
    p1 <- ggplot(DfGeneCoding) + aes(x=cv,y = median,colour=density) + geom_point(size=0.1) + 
      labs(y='Median Batch Corrected Abundance',x='Coefficient of Variation') + ggtitle('') + 
      theme_bw() + ylim(c(0,max(DfGeneNonCoding$median))) + xlim(c(0,max(DfGeneCoding$cv))) + 
      geom_vline(xintercept = quantile(CalcCV(microgliaGene$coding),percentile),colour = 'red') + 
      labs(colour='Density')
    
    
    p2 <- ggplot(DfGeneNonCoding) + aes(x=cv,y=median,colour=density) + geom_point(size=0.1) +
      labs(y='Median Batch Corrected Abundance',x='Coefficient of Variation') + ggtitle('') + 
      theme_bw() + ylim(c(0,max(DfGeneNonCoding$median))) + xlim(c(0,max(DfGeneCoding$cv)))+ 
      geom_vline(xintercept = quantile(CalcCV(microgliaGene$noncoding),percentile),colour = 'red')+ 
      labs(colour='Density')
    
    p3 <- ggplot(DfTranscriptCoding) + aes(x=cv,y = median,colour=density) + geom_point(size=0.1) + 
      labs(y='Median Batch Corrected Abundance',x='Coefficient of Variation') + ggtitle('') + 
      theme_bw() + ylim(c(0,max(DfTranscriptNonCoding$median))) + 
      xlim(c(0,max(DfTranscriptCoding$cv)))+ 
      geom_vline(xintercept = quantile(CalcCV(microgliaTranscript$coding),percentile),colour = 'red')+ 
      labs(colour='Density')
    
    p4 <- ggplot(DfTranscriptNonCoding) + aes(x=cv,y=median,colour=density) + geom_point(size=0.1) +
      labs(y='Median Batch Corrected Abundance',x='Coefficient of Variation') + ggtitle('') + 
      theme_bw() + ylim(c(0,max(DfTranscriptNonCoding$median))) +
      xlim(c(0,max(DfTranscriptCoding$cv))) +  
      geom_vline(xintercept = quantile(CalcCV(microgliaTranscript$noncoding),percentile),colour = 'red')+ 
      labs(colour='Density')
    library(ggpubr)
    ggpubr::ggexport(ggpubr::ggarrange(plotlist = list(p1,p2,p3,p4),nrow=2,ncol=2),
                     filename = '../../FinalFigures/TPMVsCV.pdf',height = 8,width = 10)
    tiff('../../FinalFigures/TPMVsCV.tiff',height = 600,width = 800)
    egg::ggarrange(plots =list(p1,p2,p3,p4),nrow=2,ncol=2)
    dev.off()
}
PlotTPMCV(list('../../Count_Data/Batch_Corrected/','../../Count_Data/TPM_Filtered/'),0.5)
