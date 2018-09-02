##########################################################################################################
# Fig 3.4.1.1 and Fig 3.4.2.2
##########################################################################################################
library(dplyr)
library(geomnet)
#Write edge list and node table (for weights)
#Input are the list of modules, and count matrix
WriteNetworkFile <- function(cluster,expMatrix){
  allGenes <- unlist(cluster$Genes)
  #Gene Name to Gene ID conversion
  load(file='../../Count_Data/geneGtfTableFull.rda')
  geneIdToName <- hashmap::hashmap(keys = geneGtfTableFull$gene_id,values = geneGtfTableFull$gene_name)
  allGeneNames <- data.frame(ensembl_gene_id=allGenes,gene_name=geneIdToName[[allGenes]])
  
  #GWAS P-value
  studyName <- cluster$StudyName
  geneScoreTbl <- data.table::fread(file=paste0('../../GWAS/PASCAL_results/KEGG_Ensembl_All/',
                                                gsub('_','.',studyName),'.sum.genescores.txt'))
  if(cluster$Biotype=='coding'){
   expMatrix <- expMatrix$coding 
  }else{
   expMatrix <- rbind(expMatrix$noncoding,expMatrix$coding)
  }
  
  #Construct node table, with GWAS p value and gene name
  nodeTbl <- left_join(data.frame(gene_id=allGenes,stringsAsFactors = F),geneScoreTbl,by=c('gene_id'='gene_id')) %>% dplyr::select(gene_id,pvalue) %>% 
    dplyr::left_join(allGeneNames,by=c('gene_id'='ensembl_gene_id'))
  nodeTbl$pvalue[is.na(nodeTbl$pvalue)] <- 1
  nodeTbl$pvalue <- -1*log10(nodeTbl$pvalue)
  
  #Construct table of edge weights
  edgeInteractions <- cor(t(expMatrix[allGenes,]))
  rownames(edgeInteractions) <- geneIdToName[[allGenes]]; colnames(edgeInteractions) <- geneIdToName[[allGenes]];
  upperTri <- which(upper.tri(edgeInteractions,diag = F),arr.ind = T)
  edgeTbl <- data.frame(Node1=character(),Node2=character(),W=numeric())
  
  for(i in 1:nrow(upperTri)){
    x <- upperTri[i,]
    edgeTbl <- rbind(edgeTbl,data.frame(Node1=geneIdToName[[allGenes[x[1]]]],
                                        Node2=geneIdToName[[allGenes[x[2]]]],
                                        W=abs(edgeInteractions[x[1],x[2]]),stringsAsFactors = F))
  }
  #Save table as txt files
  write.table(nodeTbl,file=paste0('../../Cytoscape_Networks/nodeTbl_',cluster$clusterNames,'.txt'),sep = '\t',quote = F,row.names = F)
  write.table(edgeTbl,file=paste0('../../Cytoscape_Networks/edgeTbl_',cluster$clusterNames,'.txt'),sep= '\t',quote = F,row.names = F)
  nodeTbl$gene_name <- as.character(nodeTbl$gene_name)
  microgliaGenes <- unlist(cluster$MicrogliaGenes)
  nodeTbl$Microglia <- sapply(nodeTbl$gene_name,function(x) x %in% microgliaGenes)
  return(list(nodeTbl=nodeTbl %>% dplyr::select(gene_name,pvalue,Microglia),edgeTbl=edgeTbl,adjMat=edgeInteractions))
}
#Plot module as interaction network using geom_net package
PlotCluster <- function(curCluster,save=T){
  #Write node and edge file
  load('../../Count_Data/CV_Filtered/MicrogliaGeneCVFiltered.rda')
  MSSig1 <- WriteNetworkFile(curCluster,MicrogliaGeneCVFiltered)
  library(geomnet)
  #Remove edges with small weight
  MSSig1$edgeTbl <- MSSig1$edgeTbl %>% dplyr::filter(W>0.2)
  if(length(unlist(curCluster$Genes)) > 50){
    cutoff <- quantile(MSSig1$edgeTbl$W,0.9)
    MSSig1$edgeTbl <- MSSig1$edgeTbl %>% dplyr::filter(W>cutoff)
  }
  #Construct network object
  net <- as.edgedf(MSSig1$edgeTbl) %>% dplyr::left_join(MSSig1$nodeTbl,c("from_id"="gene_name"))
  net <- fortify(as.edgedf(MSSig1$edgeTbl),MSSig1$nodeTbl)
  net$scaledW <- ((net$W)^5)*5
  #plot geomnet
  networkPlot <- ggplot(data = net,
                        aes(from_id = from_id, to_id = to_id)) +
    geom_net(aes(colour = pvalue,shape=Microglia,width=net$W),labelon = TRUE,
             ealpha = 0.25,show.legend = T,size=6,
             curvature = 0,fontsize=2.5,
             directed = F) + scale_colour_gradientn('-log(P)',colours = c('black','blue','green','orange','red'),
                                                    values = c(0,0.05,0.1,0.2,1)) + 
    theme_net() + scale_shape_discrete()+guides(shape = guide_legend(override.aes = list(size = 10))) + 
    ggtitle(curCluster$clusterNames) + theme(plot.title = element_text(size = 20, face = "bold",vjust = -3))
  if(save){
    ggpubr::ggexport(plotlist = list(ggpubr::ggarrange(networkPlot)),width = 12,height = 6,filename = paste0('../../FinalFigures/NetworkVisualizations/',curCluster$clusterNames,'.pdf'))
    ggpubr::ggexport(plotlist = list(ggpubr::ggarrange(networkPlot)),width = 1000,height = 800,filename = paste0('../../FinalFigures/NetworkVisualizations/',curCluster$clusterNames,'.tiff'))
  }else{
    networkPlot
  }
}
#Plot dendrogram of module, using average linkage
PlotClusterDendro <- function(curCluster,save=T){
  library(dplyr)
  library(dendextend)
  load('../../Count_Data/CV_Filtered/MicrogliaGeneCVFiltered.rda')
  #Get node and edge table
  clusterData <- WriteNetworkFile(curCluster,expMatrix = MicrogliaGeneCVFiltered)
  corrMat <- clusterData$adjMat
  #dissimilarity is 1-|cor|
  dissimilarityMat <- as.dist(1-abs(corrMat))
  #Average linkage method
  hc <- hclust(dissimilarityMat,method = 'average')
  #Plot dendrogram
  hcDendro <- as.dendrogram(hc,hang = 0.1)
  #join labels of dendrogram with gene names.
  labelInfo <- data.frame(gene_name=labels(hcDendro),stringsAsFactors = F) %>% 
    dplyr::left_join(clusterData$nodeTbl,by=c('gene_name'='gene_name'))
  #set shape based on whether there is microglia enrichment
  hcDendro <- hcDendro %>% set("leaves_pch", ifelse(labelInfo$Microglia, 15, 17))
  #Set color based on GWAS p value
  GetColour <- function(value){
    if(value<1){
      return('grey')
    }else if(value >= 1 & value < 2){
      return('lightgreen')
    }else if(value >= 2 & value < 3){
      return('skyblue')
    }else if(value >= 3 & value < 5){
      return('orange')
    }else if(value >= 5 & value < -log10(5e-8)){
      return('red') 
    }else{
      return('darkred')
    }
  }
  print(paste0('../../FinalFigures/ClusterDendro/',curCluster$clusterNames,'.pdf'))
  if(save){
    # pdf(file = paste0('../../FinalFigures/ClusterDendro/',curCluster$clusterNames,'.pdf'),width = 25,height = 6)
    pdf(file = paste0('../../FinalFigures/ClusterDendro/',curCluster$clusterNames,'.pdf'),width = 10,height = 6)
    
  }
  #Set plot properties
  hcDendro <- hcDendro %>% set("labels_col",sapply(labelInfo$pvalue,GetColour))
  # hcDendro <- hcDendro %>% set("labels_cex",rep(0.8,nrow(labelInfo)))
  hcDendro <- hcDendro %>% set("labels_cex",rep(1.3,nrow(labelInfo)))
  hcDendro <- hcDendro %>% set("leaves_cex",rep(1.5,nrow(labelInfo)))
  
  layout(mat = matrix(c(2,1), 
                      nrow = 2, 
                      ncol = 1),
         heights = c(3,1))
  # par(mar = c(3, 57, 4, 55))
  par(mar = c(3, 18, 4, 16))
  #Conctruct plot
  barplot(as.matrix(diff(c(0,1,2,3,5,8,11))),horiz = T,axes = F,col = c('grey','lightgreen','skyblue','orange','red','darkred'))
  axis(1, labels=c('0','1','2','3','5','>8'), at=c(0,1,2,3,5,8))
  title(main='-log10(GWAS P-value)',line = 0.6)
  par(mar = c(4, 4, 2, 0), xpd=TRUE)
  plot(hcDendro,main=paste0(curCluster$clusterNames,'-',curCluster$`KEGG Term`),ylab='1-|cor|',cex.lab=1.5)
  legend("bottom",inset=c(0,-0.24),legend = c("Gene","Microglia Enriched Gene"),ncol=2,
         col=c('black','black'), pch= c(16,17))
  if(save){
    dev.off()
  }
}
#load all modules
library(ggpubr)
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingPearson.rda')
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersNonNeurologicalPearsonCoding.rda')
sigClustersCodingPearson <- rbind(sigClustersCodingPearson,sigClustersNonNeurologicalPearsonCoding)
sigClustersCodingPearson$Method <- rep('CPG',nrow(sigClustersCodingPearson))
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersAllPearson.rda')
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersNonNeurologicalPearsonAll.rda')
sigClustersAllPearson <- rbind(sigClustersAllPearson,sigClustersNonNeurologicalPearsonAll)
sigClustersAllPearson$Method <- rep('APG',nrow(sigClustersAllPearson))
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingWGCNAUnsigned.rda')
sigClustersCodingWGCNAUnsigned$Method <- rep('CWG',nrow(sigClustersCodingWGCNAUnsigned))
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersAllWGCNAUnsigned.rda')
sigClustersAllWGCNAUnsigned$Method <- rep('AWG',nrow(sigClustersAllWGCNAUnsigned))
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingPearsonTranscripts.rda')
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingNonNeurologicalPearsonTranscripts.rda')
sigClustersCodingPearsonTranscripts <- rbind(sigClustersCodingPearsonTranscripts,
                                             sigClustersCodingNonNeurologicalPearsonTranscripts)
sigClustersCodingPearsonTranscripts$Method <- rep('CPT',nrow(sigClustersCodingPearsonTranscripts))
load('../../Count_Data/PASCAL_Results/sigClusters/sigClustersCodingWGCNATranscripts.rda')
sigClustersCodingWGCNATranscripts$Method <- rep('CWT',nrow(sigClustersCodingWGCNATranscripts))
clusterNames <-   c(paste0('CPG',1:nrow(sigClustersCodingPearson)),
                    paste0('APG',1:nrow(sigClustersAllPearson)),
                    paste0('CWG',1:nrow(sigClustersCodingWGCNAUnsigned)),
                    paste0('AWG',1:nrow(sigClustersAllWGCNAUnsigned)),
                    paste0('CPT',1:nrow(sigClustersCodingPearsonTranscripts)),
                    paste0('CWT',1:nrow(sigClustersCodingWGCNATranscripts)))
allMethodsDf <- rbind(sigClustersCodingPearson,
                      sigClustersAllPearson,
                      sigClustersCodingWGCNAUnsigned,
                      sigClustersAllWGCNAUnsigned,
                      sigClustersCodingPearsonTranscripts,
                      sigClustersCodingWGCNATranscripts)
allMethodsDf$clusterNames <- clusterNames
#Plot dendrograms and 
sapply(1:nrow(allMethodsDf),function(i) PlotCluster(allMethodsDf[i,]))
sapply(1:nrow(allMethodsDf),function(i) PlotClusterDendro(allMethodsDf[i,]))
