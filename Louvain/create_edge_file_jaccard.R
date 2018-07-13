library(Hmisc)
library(hashmap)
library(parallel)
CalcJaccard <- function(i,j,neighbourMap){
  #Find neighbours, but remove self connections
  iNeighbours <- neighbourMap[[as.character(i)]]
  jNeighbours <- neighbourMap[[as.character(j)]]
  return(length(intersect(iNeighbours,jNeighbours))/length(union(iNeighbours,jNeighbours)))
}
CreateEdgeFilesJaccard <- function(){
  if(!'MicrogliaGeneCodingCorr' %in% ls()){
    # load('../../Count_Data/CV_Filtered/MicrogliaGeneCVFiltered.rda')
    # MicrogliaGeneCodingCorr <- rcorr(t(MicrogliaGeneCVFiltered$coding),type = 'pearson')
    print('loading coding corr')
    load('../../Count_Data/Correlation_Matrices/MicrogliaGeneCodingCorr.rda')
  }
  if(!'MicrogliaGeneAllCorr' %in% ls()){
    # MicrogliaGeneAll <- rbind(MicrogliaGeneCVFiltered$coding,MicrogliaGeneCVFiltered$noncoding)
    # MicrogliaGeneAllCorr <- rcorr(t(MicrogliaGeneAll),type = 'pearson')
    print('loading all corr')
    load('../../Count_Data/Correlation_Matrices/MicrogliaGeneAllCorr.rda')
  }
  WriteFile <- function(corrList,corCutOff,pValCutOff,absoluteValue,path){
    #P value multiple correction adjustment, only taking upper diagonal
    adjPValue <- p.adjust(corrList$P[upper.tri(corrList$P,diag = F)],method = 'BH')
    adjMatrixFullP <- matrix(Inf,nrow = nrow(corrList$P),ncol = ncol(corrList$P))
    adjMatrixFullP[upper.tri(adjMatrixFullP,diag = F)] <- adjPValue
    adjMatrixFullP <- t(adjMatrixFullP)
    adjMatrixFullP[upper.tri(adjMatrixFullP,diag = F)] <- adjPValue  #make symmetric matrix, by transposing.
    
    if(absoluteValue){
      adjMatrixFull <- ifelse(abs(corrList$r)>corCutOff & adjMatrixFullP < pValCutOff,1,0)
    }else{
      adjMatrixFull <- ifelse(corrList$r>corCutOff & adjMatrixFullP < pValCutOff,1,0)
    }
    #For fast lookup, create hashmap where keys are nodes, and values are all it's neighbours, but no self connection. 
    neighbourMap <- new.env(hash = TRUE)  # equivalent to new.env()
    list2env(
      setNames(
        mclapply(1:nrow(adjMatrixFull),function(i) setdiff(which(adjMatrixFull[i,] == 1),i),mc.cores=50), 
        1:nrow(adjMatrixFull)
      ),
      envir = neighbourMap
    )
    #Convert correlation matrix to a 3 columns. 1st and 2nd columns are nodes, and 3rd column is the weight.
    adjMatrix <- adjMatrixFull
    adjMatrix[lower.tri(adjMatrix,diag = T)] <- 0 #only consider upper diagonal since symmetric. Also, no self connections.
    significantIndex <- which(adjMatrix == 1,arr.ind = T,useNames = F)
    geneNames<- rownames(corrList$r)
    print(paste('writing: ',path))
    # outFile <- file(path,'w')
    edgeTbl <- do.call(rbind,mclapply(1:nrow(significantIndex),function(i) {
      x <- significantIndex[i,];
      c(geneNames[x[1]],
        geneNames[x[2]],
        as.character(round(CalcJaccard(x[1],x[2],neighbourMap),3)))   #Use jaccard index as weight, where J = |intersect(A,B)| / |Union(A,B)|
    },mc.cores=20))
    data.table::fwrite(as.data.frame(edgeTbl),path,sep = " ",quote = F,col.names = F,row.names = F)
    # write(edgeTbl,outFile,append = F)
    # close.connection(outFile)
  }
  #microglia coding genes
  WriteFile(MicrogliaGeneCodingCorr,corCutOff=0.25,pValCutOff = 0.05,absoluteValue=T,'../../Louvain_Edge_List/Jaccard/CodingGenesEdgeListMicroglia_Jaccard_pval0p05_cor0p25_abs.txt')

  #microglia all genes
  WriteFile(MicrogliaGeneAllCorr,corCutOff=0.25,pValCutOff = 0.05,absoluteValue=T,'../../Louvain_Edge_List/Jaccard/AllGenesEdgeListMicroglia_Jaccard_pval0p05_cor0p25_abs.txt')
  
  #microglia coding genes
  WriteFile(MicrogliaGeneCodingCorr,corCutOff=0.25,pValCutOff = 0.05,absoluteValue=F,'../../Louvain_Edge_List/Jaccard/CodingGenesEdgeListMicroglia_Jaccard_pval0p05_cor0p25_noabs.txt')
  
  #microglia all genes
  WriteFile(MicrogliaGeneAllCorr,corCutOff=0.25,pValCutOff = 0.05,absoluteValue=F,'../../Louvain_Edge_List/Jaccard/AllGenesEdgeListMicroglia_Jaccard_pval0p05_cor0p25_noabs.txt')
  
  #microglia coding genes
  sortedCor <- sort(MicrogliaGeneCodingCorr$r[upper.tri(MicrogliaGeneCodingCorr$r,diag = F)],decreasing=F)
  top1milCutOff <- sortedCor[length(sortedCor) - 1e6]
  WriteFile(MicrogliaGeneCodingCorr,corCutOff=top1milCutOff,pValCutOff = Inf,absoluteValue=F,'../../Louvain_Edge_List/Jaccard/CodingGenesEdgeListMicroglia_Jaccard_top1milpos.txt')
  
  #microglia all genes
  sortedCor <- sort(MicrogliaGeneAllCorr$r[upper.tri(MicrogliaGeneAllCorr$r,diag = F)],decreasing=F)
  top1milCutOff <- sortedCor[length(sortedCor) - 1e6]
  WriteFile(MicrogliaGeneAllCorr,corCutOff=top1milCutOff,pValCutOff = Inf,absoluteValue=F,'../../Louvain_Edge_List/Jaccard/AllGenesEdgeListMicroglia_Jaccard_top1milpos.txt')
  
}
