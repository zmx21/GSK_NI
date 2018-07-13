WriteSettingsFile <- function(outName,outDir,type){
  if(type=='all'){
    rawSetting <- readLines('DREAM_Settings_All.txt')
  }else{
    rawSetting <- readLines('DREAM_Settings_Coding.txt')
  }
  rawSetting[grep('outputDirectory',rawSetting)] <- paste0('outputDirectory = ',outDir)
  rawSetting[grep('writeUsedSettings',rawSetting)] <- paste0('writeUsedSettings = ',outDir,'settingsOut.txt')
  writeLines(rawSetting,con = file(outName))
}
RunPASCAL <- function(geneSetDir,PASCALPath,outPath,type=c('coding','all')){
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
  mclapply(1:length(geneSets),function(i){
    outDir <- paste0(outPath,gsub(x=geneSets[i],pattern = '.gmt',replacement = '/'))
    system(paste0('mkdir ',outDir))
    settingsName <- gsub(x=geneSets[i],pattern = '.gmt',replacement = '_settings.txt')
    WriteSettingsFile(settingsName,outDir,type = ifelse(grepl(pattern = 'All',x = geneSets[i]),'all','coding'))
    cmd <- paste('find ../parsed_studies/*txt | parallel -j 9 ./run_PASCAL {1}',settingsName,geneSetPath[i],outDir,sep = ' ')
    system(cmd,intern = T)
  },mc.cores = 2)
}
RunPASCAL(geneSetDir = './resources/genesets/',
          PASCALPath = '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_New/',
          outPath = '../PASCAL_results2/',
          type = 'all')
# RunPASCAL(geneSetDir = './resources/genesets/WGCNA_clusters/',
#           PASCALPath = '/local/data/public/zmx21/zmx21_private/GSK/GWAS/PASCAL_New/',
#           outPath = '../PASCAL_results_WGCNA/')
