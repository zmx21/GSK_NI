#Pass in data frame from PASCAL result, for each cluster, append a corrected by value
#Do this for each leve, each study, and each biotype
AppendCorrectedPVal <- function(resultDf){
  library(dplyr)
  uniqueLevels <- unique(resultDf$Level)
  uniqueBiotypes <- unique(resultDf$Biotype)
  uniqueStudies <- unique(resultDf$StudyName)
  allGroups <- expand.grid(Level = uniqueLevels,Biotype=uniqueBiotypes,StudyName = uniqueStudies)
  subGroupDf <- vector(mode='list',length = nrow(allGroups))
  for(i in 1:nrow(allGroups)){
    #Get current study at the same level and same biotype
    currentDf <- resultDf %>% dplyr::filter(Level==allGroups$Level[i],
                                            Biotype==allGroups$Biotype[i],
                                            StudyName==allGroups$StudyName[i])
    #Do FDR p-val adj using BH procedure
    currentDf$adjPvalue <- p.adjust(currentDf$chi2Pvalue,method = 'BH')
    subGroupDf[[i]] <- currentDf
  }
  return(do.call(rbind,subGroupDf) %>% dplyr::arrange(adjPvalue))
}
