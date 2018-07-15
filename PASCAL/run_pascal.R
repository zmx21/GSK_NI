WriteSettingsFile <- function(outName,outDir,type,calcGeneScore){
  if(type=='all'){
    rawSetting <- readLines('DREAM_Settings_All.txt')
  }else{
    rawSetting <- readLines('DREAM_Settings_Coding.txt')
  }
  rawSetting[grep('outputDirectory',rawSetting)] <- paste0('outputDirectory = ',outDir)
  rawSetting[grep('writeUsedSettings',rawSetting)] <- paste0('writeUsedSettings = ',outDir,'settingsOut.txt')
  if(calcGeneScore){
    rawSetting[grep('loadSingleGeneScoresFromFiles',rawSetting)] <- 'loadSingleGeneScoresFromFiles = 0'
  }
  writeLines(rawSetting,con = file(outName))
}
RunPASCAL <- function(geneSetDir,PASCALPath,outPath,type=c('coding','all'),calcGeneScore=F){
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
  if(calcGeneScore){
    mclapply(1:length(geneSets),function(i){
      outDir <- paste0(outPath,gsub(x=geneSets[i],pattern = '.gmt',replacement = '/'))
      system(paste0('mkdir ',outDir))
      settingsName <- gsub(x=geneSets[i],pattern = '.gmt',replacement = '_settings.txt')
      WriteSettingsFile(settingsName,outDir,type = ifelse(grepl(pattern = 'All',x = geneSets[i]),'all','coding'),calcGeneScore=T)
      cmd <- paste('find ../parsed_studies/*txt | parallel -j 11 ./run_PASCAL_genescore {1}',settingsName,geneSetPath[i],outDir,sep = ' ')
      system(cmd,intern = T)
    },mc.cores = 1)
  }else{
    mclapply(1:length(geneSets),function(i){
      outDir <- paste0(outPath,gsub(x=geneSets[i],pattern = '.gmt',replacement = '/'))
      system(paste0('mkdir ',outDir))
      settingsName <- gsub(x=geneSets[i],pattern = '.gmt',replacement = '_settings.txt')
      WriteSettingsFile(settingsName,outDir,type = ifelse(grepl(pattern = 'All',x = geneSets[i]),'all','coding'),calcGeneScore=F)
      cmd <- paste('find ../parsed_studies/*txt | parallel -j 11 ./run_PASCAL {1}',settingsName,geneSetPath[i],outDir,sep = ' ')
      system(cmd,intern = T)
    },mc.cores = 1)
  }
}
RunPASCAL(geneSetDir = './resources/genesets/Jaccard_Cor0p2/',
          PASCALPath = '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_New/',
          outPath = '../PASCAL_results/Jaccard_Cor0p2/',
          type = c('coding'),calcGeneScore=T)
RunPASCAL(geneSetDir = './resources/genesets/Pearson_Cor0p2/',
          PASCALPath = '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_New/',
          outPath = '../PASCAL_results/Pearson_Cor0p2/',
          type = c('coding'),calcGeneScore=F)
RunPASCAL(geneSetDir = './resources/genesets/WGCNA_size3/',
          PASCALPath = '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_New/',
          outPath = '../PASCAL_results/WGCNA_size3/',
          type = c('coding'),calcGeneScore=F)

RunPASCAL(geneSetDir = './resources/genesets/Jaccard_Cor0p2/',
          PASCALPath = '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_New/',
          outPath = '../PASCAL_results/Jaccard_Cor0p2/',
          type = c('all'),calcGeneScore=T)
RunPASCAL(geneSetDir = './resources/genesets/Pearson_Cor0p2/',
          PASCALPath = '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_New/',
          outPath = '../PASCAL_results/Pearson_Cor0p2/',
          type = c('all'),calcGeneScore=F)
RunPASCAL(geneSetDir = './resources/genesets/WGCNA_size3/',
          PASCALPath = '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_New/',
          outPath = '../PASCAL_results/WGCNA_size3/',
          type = c('all'),calcGeneScore=F)


