load('../../Count_Data/Batch_Corrected/SalmonTPM_Combat_ExpCorrected.rda')
#load('../../Count_Data/SalmonTPM_Gene_Merged_GeneFilt.rda')
load('../../Count_Data/geneGtfTableFull.rda')
proteinCodingGenes <- geneGtfTableFull %>% dplyr::filter(gene_biotype=='protein_coding') %>% dplyr::select(gene_id) %>% t() %>% as.vector()

batchCorrected <- SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged
proteinCodingBatchCorrected <- batchCorrected[rownames(batchCorrected)%in%proteinCodingGenes,]
cvBatchCorrected <- apply(proteinCodingBatchCorrected,1,function(x) sd(x)/mean(x))
#hist(cvBatchCorrected,breaks=50,main='Batch Corrected TPM - CV Distribution',xlab='CV')

batchCorrectedLowCV <- proteinCodingBatchCorrected[cvBatchCorrected < quantile(cvBatchCorrected,0.1),]
batchCorrectedHighCV <- proteinCodingBatchCorrected[cvBatchCorrected > quantile(cvBatchCorrected,0.9),]

batchCorrectedLowCVGenes <- data.frame(gene_id=rownames(batchCorrectedLowCV),stringsAsFactors = F) %>% 
  left_join(dplyr::select(geneGtfTableFull,gene_id,gene_name)) %>%
  dplyr::select(gene_name) %>% t() %>% as.vector()
batchCorrectedHighCVGene <- data.frame(gene_id=rownames(batchCorrectedHighCV),stringsAsFactors = F) %>% 
  left_join(dplyr::select(geneGtfTableFull,gene_id,gene_name)) %>%
  dplyr::select(gene_name) %>% t() %>% as.vector()
write(batchCorrectedLowCVGenes,file='../../CV_GOAnno/lowCVGenes.txt')
write(batchCorrectedHighCVGene,file='../../CV_GOAnno/highCVGenes.txt')

  
library(enrichR)
batchCorrectedLowCVGO <- enrichr(batchCorrectedLowCVGenes,databases = c('GO_Molecular_Function_2017','GO_Biological_Process_2017','GO_Cellular_Component_2017'))
batchCorrectedHighCVGO <- enrichr(batchCorrectedHighCVGene,databases = c('GO_Molecular_Function_2017','GO_Biological_Process_2017','GO_Cellular_Component_2017'))

batchCorrectedHighCVGO$GO_Molecular_Function_2017$Overlap <- sapply(batchCorrectedHighCVGO$GO_Molecular_Function_2017$Overlap,function(x) gsub('/',' out of ',x))
batchCorrectedLowCVGO$GO_Molecular_Function_2017$Overlap <- sapply(batchCorrectedLowCVGO$GO_Molecular_Function_2017$Overlap,function(x) gsub('/',' out of ',x))

batchCorrectedHighCVGO$GO_Biological_Process_2017$Overlap <- sapply(batchCorrectedHighCVGO$GO_Biological_Process_2017$Overlap,function(x) gsub('/',' out of ',x))
batchCorrectedLowCVGO$GO_Biological_Process_2017$Overlap <- sapply(batchCorrectedLowCVGO$GO_Biological_Process_2017$Overlap,function(x) gsub('/',' out of ',x))

batchCorrectedHighCVGO$GO_Cellular_Component_2017$Overlap <- sapply(batchCorrectedHighCVGO$GO_Cellular_Component_2017$Overlap,function(x) gsub('/',' out of ',x))
batchCorrectedLowCVGO$GO_Cellular_Component_2017$Overlap <- sapply(batchCorrectedLowCVGO$GO_Cellular_Component_2017$Overlap,function(x) gsub('/',' out of ',x))


library(openxlsx)
write.xlsx(list(batchCorrectedLowCVGO$GO_Molecular_Function_2017,batchCorrectedLowCVGO$GO_Biological_Process_2017,batchCorrectedLowCVGO$GO_Cellular_Component_2017),file='../../batchCorrectedLowCVGO.xlsx',sheetName=c('Molecular Process','Biological Process','Cellular Component'))
write.xlsx(list(batchCorrectedHighCVGO$GO_Molecular_Function_2017,batchCorrectedHighCVGO$GO_Biological_Process_2017,batchCorrectedHighCVGO$GO_Cellular_Component_2017),file='../../batchCorrectedHighCVGO.xlsx',sheetName=c('Molecular Process','Biological Process','Cellular Component'))



batchCorrectedLowCVPanther <- read.table(file='../../CV_GOAnno/lowCVPanther.txt',sep='\t',header=F)
batchCorrectedHighCVPanther <- read.table(file='../../CV_GOAnno/highCVPanther.txt',sep='\t',header=F)
PantherDf <- data.frame(CV=rep(c('LowCV','HighCV'),each=nrow(batchCorrectedHighCVPanther)),
                        GOAnnotation=c(batchCorrectedLowCVPanther$V2,batchCorrectedHighCVPanther$V2),
                        Fraction=c(as.numeric(gsub('%','',batchCorrectedLowCVPanther$V5)),as.numeric(gsub('%','',batchCorrectedHighCVPanther$V5))))
ggplot(PantherDf) + aes(x=GOAnnotation,y=Fraction,fill=CV) + geom_bar(stat="identity", position=position_dodge()) + coord_flip() + 
  labs(y='Fraction of GO Annotation Hits',x='GO Annotation Category') + ggtitle('Categories of Genes with Low CV vs High CV')

#PercentileAbundance <- apply(SalmonTPM_Gene_Merged[rownames(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged),],2,function(x) rank(x)/length(x))
#cvPercentile <- apply(PercentileAbundance,1,function(x) sd(x)/mean(x))

# percentileCorrectedLowCV <- PercentileAbundance[cvPercentile < quantile(cvPercentile,0.25),]
# percentileCorrectedHighCV <- PercentileAbundance[cvPercentile > quantile(cvPercentile,0.75),]

# percentileCorrectedLowCVGenes <- data.frame(gene_id=rownames(percentileCorrectedLowCV),stringsAsFactors = F) %>% 
#   left_join(dplyr::select(geneGtfTableFull,gene_id,gene_name)) %>% 
#   dplyr::select(gene_name) %>% t() %>% as.vector()
# percentileCorrectedHighCVGene <- data.frame(gene_id=rownames(percentileCorrectedHighCV),stringsAsFactors = F) %>% 
#   left_join(dplyr::select(geneGtfTableFull,gene_id,gene_name)) %>% 
#   dplyr::select(gene_name) %>% t() %>% as.vector()
# percentileCorrectedLowCVGO <- enrichr(percentileCorrectedLowCVGenes,databases = c('GO_Molecular_Function_2017','GO_Biological_Process_2017','GO_Cellular_Component_2017'))
# percentileCorrectedHighCVGO <- enrichr(percentileCorrectedHighCVGene,databases = c('GO_Molecular_Function_2017','GO_Biological_Process_2017','GO_Cellular_Component_2017'))
# percentileCorrectedLowCVGO$GO_Molecular_Function_2017$Overlap <- sapply(percentileCorrectedLowCVGO$GO_Molecular_Function_2017$Overlap,function(x) gsub('/',' out of ',x))
# percentileCorrectedHighCVGO$GO_Molecular_Function_2017$Overlap <- sapply(percentileCorrectedHighCVGO$GO_Molecular_Function_2017$Overlap,function(x) gsub('/',' out of ',x))
# percentileCorrectedLowCVGO$GO_Biological_Process_2017$Overlap <- sapply(percentileCorrectedLowCVGO$GO_Biological_Process_2017$Overlap,function(x) gsub('/',' out of ',x))
# percentileCorrectedHighCVGO$GO_Biological_Process_2017$Overlap <- sapply(percentileCorrectedHighCVGO$GO_Biological_Process_2017$Overlap,function(x) gsub('/',' out of ',x))
# percentileCorrectedLowCVGO$GO_Cellular_Component_2017$Overlap <- sapply(percentileCorrectedLowCVGO$GO_Cellular_Component_2017$Overlap,function(x) gsub('/',' out of ',x))
# percentileCorrectedHighCVGO$GO_Cellular_Component_2017$Overlap <- sapply(percentileCorrectedHighCVGO$GO_Cellular_Component_2017$Overlap,function(x) gsub('/',' out of ',x))
# write.xlsx(list(percentileCorrectedLowCVGO$GO_Molecular_Function_2017,percentileCorrectedLowCVGO$GO_Biological_Process_2017,percentileCorrectedLowCVGO$GO_Cellular_Component_2017),file='../../percentileCorrectedLowCVGO.xlsx',sheetName=c('Molecular Process','Biological Process','Cellular Component'))
# write.xlsx(list(percentileCorrectedHighCVGO$GO_Molecular_Function_2017,percentileCorrectedHighCVGO$GO_Biological_Process_2017,percentileCorrectedHighCVGO$GO_Cellular_Component_2017),file='../../percentileCorrectedHighCVGO.xlsx',sheetName=c('Molecular Process','Biological Process','Cellular Component'))

# df <- data.frame(cvBatchCorrected=cvBatchCorrected,cvPercentile=cvPercentile)
# get_density <- function(x, y, n = 100) {
#   dens <- MASS::kde2d(x = x, y = y, n = n)
#   ix <- findInterval(x, dens$x)
#   iy <- findInterval(y, dens$y)
#   ii <- cbind(ix, iy)
#   return(dens$z[ii])
# }
# 
# library(ggplot2)
# library(viridis)
# df$density <- get_density(df$cvBatchCorrected,df$cvPercentile)
# ggplot(df,aes(x=cvBatchCorrected,y=cvPercentile,colour=density)) +
#   geom_point(size=0.5) + scale_y_continuous(limits = c(0,5))+
#   theme_bw() +  scale_color_viridis()
# 
