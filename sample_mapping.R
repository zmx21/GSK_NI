GetSampleMapping <- function(){
  runTable <- read.table(file = '/local/data/public/zmx21/zmx21_private/GSK/Galatro/SraRunTable.txt',header = T,sep = '\t')
  SRRToGSM <- data.frame(GSM = runTable$Sample_Name,SRR = runTable$Run)
  matrixFile <- readLines(con = "/local/data/public/zmx21/zmx21_private/GSK/Galatro/GSE99074-GPL16791_series_matrix.txt")
  sampleTitles <- matrixFile[which(sapply(matrixFile,function(x) grepl('Sample_title',x)))]
  sampleTitles <- unlist(strsplit(sampleTitles,'\t'))
  sampleTitles <- sampleTitles[2:length(sampleTitles)]
  sampleTitles <- sapply(sampleTitles,function(x) paste0(unlist(strsplit(x,""))[2:(length(unlist(strsplit(x,"")))-1)],collapse = ""))
  
  sampleGSM <- matrixFile[which(sapply(matrixFile,function(x) grepl('Sample_geo_accession',x)))]
  sampleGSM <- unlist(strsplit(sampleGSM,'\t'))
  sampleGSM <- sampleGSM[2:length(sampleGSM)]
  sampleGSM <- sapply(sampleGSM,function(x) paste0(unlist(strsplit(x,""))[2:(length(unlist(strsplit(x,"")))-1)],collapse = ""))
  
  titleToGSM <- data.frame(title=sampleTitles,GSM = sampleGSM)
  mapping <- merge(SRRToGSM,titleToGSM)
  return(mapping)
}