path <- '/local/data/public/zmx21/zmx21_private/GSK/Galatro_Brain/fastq_trimmed/'
# path <- '/local/data/public/zmx21/zmx21_private/GSK/Galatro/fastq_untrimmed_merged'
setwd(path)
allDir <- dir(path)
allFastq <- allDir[grep('.fq',allDir)]

metaData <- read.table('/local/data/public/zmx21/zmx21_private/GSK/Galatro_Brain/SraRunTable.txt',header = T,sep='\t')
# metaData <- read.table('/local/data/public/zmx21/zmx21_private/GSK/Galatro/SraRunTable.txt',header = T,sep='\t')
allExperiments <- unique(metaData$Sample_Name)
library(parallel)
# for(i in 1:length(allExperiments)){
mclapply(1:length(allExperiments),function(i){
  currentExperiment <- allExperiments[i]
  currentRuns <- sort(metaData$Run[metaData$Sample_Name==currentExperiment])
  if(length(currentRuns) > 1){
    command <- paste0('mv ',currentRuns[1],'_1_val_1.fq ',currentExperiment,'_1.fq')
    system(command = command)
    command <- paste0('cat ',currentRuns[2],'_1_val_1.fq >> ',currentExperiment,'_1.fq')
    system(command = command)
    command <- paste0('mv ',currentRuns[1],'_2_val_2.fq ',currentExperiment,'_2.fq')
    system(command = command)
    command <- paste0('cat ',currentRuns[2],'_2_val_2.fq >> ',currentExperiment,'_2.fq')
    system(command = command)
  }else{
    command <- paste0('mv ',currentRuns[1],'_1_val_1.fq ',currentExperiment,'_1.fq')
    system(command = command)
    command <- paste0('mv ',currentRuns[1],'_2_val_2.fq ',currentExperiment,'_2.fq')
    system(command = command)
  }
},mc.cores=42)
# }