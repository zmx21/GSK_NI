WriteSettingsFile <- function(outName,outDir,type,calcGeneScore){
  if(type=='all'){
    rawSetting <- readLines('DREAM_Settings_All_Transcripts.txt')
  }else{
    rawSetting <- readLines('DREAM_Settings_Coding_Transcripts.txt')
  }
  rawSetting[grep('outputDirectory',rawSetting)] <- paste0('outputDirectory = ',outDir)
  rawSetting[grep('writeUsedSettings',rawSetting)] <- paste0('writeUsedSettings = ',outDir,'settingsOut.txt')
  if(calcGeneScore){
    rawSetting[grep('loadSingleGeneScoresFromFiles',rawSetting)] <- 'loadSingleGeneScoresFromFiles = 0'
  }
  writeLines(rawSetting,con = file(outName))
}
RunPASCAL <- function(geneSetDir,PASCALPath,outPath,type=c('coding','all'),calcGeneScore=F,study=NULL){
  library(parallel)
  setwd(PASCALPath)
  geneSets <- dir(geneSetDir,pattern = 'gmt')
  if(length(type) == 1){
    if(type=='coding'){
      geneSets <- geneSets[sapply(geneSets,function(x) grepl(pattern = 'Coding',x = x))]
    }else if(type=='all'){
      geneSets <- geneSets[sapply(geneSets,function(x) grepl(pattern = 'All',x = x))]
    }
  }
  geneSetPath <- paste0(geneSetDir,geneSets)
  print(geneSetPath)
  if(!is.null(study)){
    mclapply(1:length(geneSets),function(i){
      outDir <- paste0(outPath,gsub(x=geneSets[i],pattern = '.gmt',replacement = '/'))
      system(paste0('mkdir -p ',outDir))
      settingsName <- gsub(x=geneSets[i],pattern = '.gmt',replacement = '_settings.txt')
      WriteSettingsFile(settingsName,outDir,type = ifelse(grepl(pattern = 'All',x = geneSets[i]),'all','coding'),calcGeneScore=T)
      cmd <- paste(paste0('find ../parsed_studies/',study),' | xargs -n 1 -i ./run_PASCAL_genescore_transcripts {}',settingsName,geneSetPath[i],outDir,'> log.out',sep = ' ')
      system(cmd,intern = T)
    },mc.cores = 1)
    return()
  }
  if(calcGeneScore){
    mclapply(1:length(geneSets),function(i){
      outDir <- paste0(outPath,gsub(x=geneSets[i],pattern = '.gmt',replacement = '/'))
      system(paste0('mkdir -p ',outDir))
      settingsName <- gsub(x=geneSets[i],pattern = '.gmt',replacement = '_settings.txt')
      WriteSettingsFile(settingsName,outDir,type = ifelse(grepl(pattern = 'All',x = geneSets[i]),'all','coding'),calcGeneScore=T)
      cmd <- paste('find ../parsed_studies/*txt | parallel -j 13 --memfree 50G ./run_PASCAL_genescore_transcripts {1}',settingsName,geneSetPath[i],outDir,sep = ' ')
      system(cmd,intern = T)
    },mc.cores = 1)
  }else{
    mclapply(1:length(geneSets),function(i){
      outDir <- paste0(outPath,gsub(x=geneSets[i],pattern = '.gmt',replacement = '/'))
      system(paste0('mkdir -p ',outDir))
      settingsName <- gsub(x=geneSets[i],pattern = '.gmt',replacement = '_settings.txt')
      WriteSettingsFile(settingsName,outDir,type = ifelse(grepl(pattern = 'All',x = geneSets[i]),'all','coding'),calcGeneScore=F)
      cmd <- paste('find ../parsed_studies/*txt | parallel -j 13 --memfree 50G ./run_PASCAL_transcripts {1}',settingsName,geneSetPath[i],outDir,sep = ' ')
      system(cmd,intern = T)
    },mc.cores = 1)
  }
}
# RunPASCAL(geneSetDir = './resources/genesets/Pearson_Cor0p2_Transcripts/',
#           PASCALPath = '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_New/',
#           outPath = '../PASCAL_results/Pearson_Cor0p2_Transcripts/',
#           type = c('coding'),calcGeneScore=T)

# RunPASCAL(geneSetDir = './resources/genesets/Pearson_Cor0p2_Transcripts/',
#           PASCALPath = '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_New/',
#           outPath = '../PASCAL_results/Pearson_Cor0p2_Transcripts/',
#           type = c('coding'),calcGeneScore=T,study = 'tau.Deming.txt')

RunPASCAL(geneSetDir = './resources/genesets/WGCNA_size3_transcripts/',
          PASCALPath = '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_New/',
          outPath = '../PASCAL_results/WGCNA_size3_transcripts/',
          type = c('coding'),calcGeneScore=F)

# RunPASCAL(geneSetDir = './resources/genesets/Pearson_Cor0p2_Transcripts/',
#           PASCALPath = '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_New/',
#           outPath = '../PASCAL_results/Pearson_Cor0p2_Transcripts/',
#           type = c('all'),calcGeneScore=T)

