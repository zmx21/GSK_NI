library(Hmisc)
CreateEdgeFilesSigned <- function(){
  if(!'MicrogliaGeneCodingCorr' %in% ls()){
    load('../../Count_Data/CV_Filtered/MicrogliaGeneCVFiltered.rda')
    MicrogliaGeneCodingCorr <- rcorr(t(MicrogliaGeneCVFiltered$coding),type = 'pearson')
  }
  if(!'MicrogliaGeneAllCorr' %in% ls()){

    MicrogliaGeneAll <- rbind(MicrogliaGeneCVFiltered$coding,MicrogliaGeneCVFiltered$noncoding)
    MicrogliaGeneAllCorr <- rcorr(t(MicrogliaGeneAll),type = 'pearson')
  }
  WriteFile <- function(corrList,pValCutOff,path){
    #Convert correlation matrix to a 3 columns. 1st and 2nd columns are nodes, and 3rd column is the weight.
    corrList$P[lower.tri(corrList$P,diag = T)] <- Inf #only consider upper diagonal since symmetric
    significantIndex <- which(corrList$P < pValCutOff,arr.ind = T,useNames = F)
    geneNames<- rownames(corrList$r)
    print(paste('writing: ',path))
    outFile <- file(path,'w')
    for(i in 1:nrow(significantIndex)){
      x <- significantIndex[i,]
      currentRow <- c(geneNames[x[1]],
                      geneNames[x[2]],
                      as.character(round((corrList$r[x[1],x[2]] + 1)/2,3)))   #unsigned network as (cor + 1) / 2
      write(currentRow,outFile,append = T,ncolumns = 4)
    }
    close.connection(outFile)
  }
  
  #microglia coding genes
  WriteFile(MicrogliaGeneCodingCorr,pValCutOff=0.05,'../../Louvain_Edge_List/CodingGenesEdgeListMicrogliaSigned.txt')

  #microglia all genes
  WriteFile(MicrogliaGeneAllCorr,pValCutOff=0.05,'../../Louvain_Edge_List/AllGenesEdgeListMicrogliaSigned.txt')
}
