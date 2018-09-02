AppendAdjustedPValue <- function(df){
  uniqueLevels <- unique(df$Level)
  uniqueStudy <- unique(df$StudyName)
  uniqueBiotype <- unique(df$Biotype)
  
  subGroups <- expand.grid(Level=uniqueLevels,StudyName=uniqueStudy,Biotype=uniqueBiotype)
  adjustedDf <- vector(mode='list',length = nrow(subGroups))
  for(i in 1:nrow(subGroups)){
    curDf <- df %>% dplyr::filter(Level==subGroups$Level[i],StudyName==subGroups$StudyName[i],Biotype==subGroups$Biotype[i])
    curDf$adjPval <- p.adjust(curDf$empPvalue,method = 'BH')
    adjustedDf[[i]] <- curDf
  }
  return(do.call(rbind,adjustedDf) %>% dplyr::arrange(adjPval))
}