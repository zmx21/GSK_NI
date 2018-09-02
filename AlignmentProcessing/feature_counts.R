######################################################################################################
#Aggegrate STAR output to Transcript and Exon level counts using FeatureCounts.
######################################################################################################
library(Rsubread)
ImportSTARCounts <- function(path){
  allDir <- dir(path)
  #All ouput samples
  runDir <- allDir[union(grep('SRR',allDir),grep('GSM',allDir))]
  #Get paths for all samples
  allPaths <- sapply(runDir,function(x) (paste0(path,'/',x,'/Aligned.out.sam')))
  
  GTFPath <- '/local/data/public/zmx21/zmx21_private/GSK/GRCh37_Ensembl75/Homo_sapiens.GRCh37.75.gtf'
  
  #Summarize exon level into gene level and transcript level. 
  #useMetaFeatures is thus TRUE. Disallow fractional counts. 
  geneLevelCounts <- Rsubread::featureCounts(files=allPaths,isPairedEnd = T,annot.ext=GTFPath,isGTFAnnotationFile=TRUE,
                                   GTF.featureType = 'exon',GTF.attrType = 'gene_id',requireBothEndsMapped=T,fraction = F,
                                   allowMultiOverlap=T,useMetaFeatures=T,nthreads=10)
  transcriptLevelCounts <- Rsubread::featureCounts(files=allPaths,isPairedEnd = T,annot.ext=GTFPath,isGTFAnnotationFile=TRUE,
                                             GTF.featureType = 'exon',GTF.attrType = 'transcript_id',requireBothEndsMapped=T,fraction = F,
                                             allowMultiOverlap=T,useMetaFeatures=T,nthreads=10)
  
  #Need to allowMultiOverlap, where all exons are assigned a count if a read overlaps multiple exons. 
  #Feature type is exon, don't do summarizations, and output the exon_id from the GTF file. Disallow fractional counts
  exonLevelCounts <- Rsubread::featureCounts(files=allPaths,isPairedEnd = T,annot.ext=GTFPath,isGTFAnnotationFile=TRUE,
                                             GTF.featureType = 'exon',GTF.attrType = 'exon_id',useMetaFeatures=F,nthreads = 10,
                                             allowMultiOverlap=T,fraction = F)
  #Name count matrix
  colnames(geneLevelCounts$counts) <- runDir; colnames(geneLevelCounts$stat) <- c('stat_type',runDir); geneLevelCounts$targets <- runDir
  colnames(exonLevelCounts$counts) <- runDir; colnames(exonLevelCounts$stat) <- c('stat_type',runDir); exonLevelCounts$targets <- runDir
  colnames(transcriptLevelCounts$counts) <- runDir; colnames(transcriptLevelCounts$stat) <- c('stat_type',runDir); transcriptLevelCounts$targets <- runDir
  
  #Save count matrix
  save(transcriptLevelCounts,file='/local/data/public/zmx21/zmx21_private/GSK/Count_Data/STARCounts_TranscriptLevel_Microglia.rda')
  save(geneLevelCounts,file='/local/data/public/zmx21/zmx21_private/GSK/Count_Data/STARCounts_GeneLevel_Microglia.rda')
  save(exonLevelCounts,file='/local/data/public/zmx21/zmx21_private/GSK/Count_Data/STARCounts_ExonLevel_Microglia.rda')
  
}
#Plot statistics, by parsing output of featureCounts. 
PlotStats <- function(stats){
  #Remove description col
  stats <- stats[,-1]
  
  sampleNames <- colnames(stats)
  assignedRate <- stats[1,] / colSums(stats) * 100
  unmappedRate <- stats[2,] / colSums(stats) * 100
  multMapRate <- stats[7,] / colSums(stats) * 100
  noFeatureRate <- stats[10,] / colSums(stats) * 100
  mappingRates <- rbind(noFeatureRate,multMapRate,unmappedRate,assignedRate)
  dat <- as.data.frame(mappingRates)
  library(reshape2)
  dat$RateType <- factor(rev(c('assigned rate','unmapped rate','multi mapping rate','no feature rate')),
                         levels = rev(c('assigned rate','unmapped rate','multi mapping rate','no feature rate')))
  dat2 <- melt(dat, id.vars = "RateType")
  
  library(ggplot2)
  
  ggplot(dat2, aes(x=variable, y=value, fill=RateType)) + 
    geom_bar(stat="identity") +
    coord_flip()+
    xlab("Sample") +
    ylab("Percentage") +
    theme(legend.position = "top")
  
}

path <- '/local/data/public/zmx21/zmx21_private/GSK/Galatro/STAR_aligned_merged/'
# path <- '/local/data/public/zmx21/zmx21_private/GSK/Galatro_Brain/STAR_aligned_whole_brain'
ImportSTARCounts(path)
# PlotStats(geneLevelCounts$stat)
