library(Rsubread)
ImportSTARCounts <- function(path){
  allDir <- dir(path)
  runDir <- allDir[union(grep('SRR',allDir),grep('GSM',allDir))]
  #Get paths for all samples
  allPaths <- sapply(runDir,function(x) (paste0(path,'/',x,'/Aligned.out.sam')))
  
  GTFPath <- '/local/data/public/zmx21/zmx21_private/GSK/GRCh38/Homo_sapiens.GRCh38.92.gtf'
  #GTFPath <- '/local/data/public/zmx21/zmx21_private/GSK/GRCh37_Ensembl75/Homo_sapiens.GRCh37.75.gtf'
  
  #Summarize exon level into gene level. useMetaFeatures is thus TRUE. Same parameters as Galatro et al.
  geneLevelCounts <- Rsubread::featureCounts(files=allPaths,isPairedEnd = T,annot.ext=GTFPath,isGTFAnnotationFile=TRUE,
                                   GTF.featureType = 'exon',GTF.attrType = 'gene_id',requireBothEndsMapped=T,allowMultiOverlap=T,fraction=T,useMetaFeatures=T,nthreads=30)
  
  #Need to allowMultiOverlap, where all exons are assigned a count if a read overlaps multiple exons. 
  #Feature type is exon, don't do summarizations, and output the exon_id from the GTF file
  #exonLevelCounts <- Rsubread::featureCounts(files=allPaths,isPairedEnd = T,annot.ext=GTFPath,isGTFAnnotationFile=TRUE,
                                   #GTF.featureType = 'exon',GTF.attrType = 'exon_id',useMetaFeatures=F,nthreads = 10,allowMultiOverlap=T)
  
  colnames(geneLevelCounts$counts) <- runDir; colnames(geneLevelCounts$stat) <- c('stat_type',runDir); geneLevelCounts$targets <- runDir
  #colnames(exonLevelCounts$counts) <- runDir; colnames(exonLevelCounts$stat) <- c('stat_type',runDir); exonLevelCounts$targets <- runDir
  
  save(geneLevelCounts,file='../STARCounts_GeneLevel_GRCh38.rda')
  #save(exonLevelCounts,file='../STARCounts_ExonLevel_GRCh38.rda')
  
}
PlotStats <- function(stats){
  #Remove descriptio col
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

# path <- '/local/data/public/zmx21/zmx21_private/GSK/Galatro/STAR_aligned_merged_GRCh38/'
# ImportSTARCounts(path)
PlotStats(geneLevelCounts$stat)