##############################################################################################
#Fig. 3.1.0.2
##############################################################################################

load('../../Count_Data/Batch_Corrected/SalmonTPM_Combat_ExpCorrected.rda')
GeneLevelCountMatrix <- SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged
#import data
microgliaGenes <- data.table::fread('../../Count_Data/Galatro_Microglia_Core_Genes.txt',header = F)
microgliaGenes <- unique(microgliaGenes$V2)

housekeepingGenes <- data.table::fread('../../Count_Data/HK_genes.txt',header = F)
housekeepingGenes <- unique(housekeepingGenes$V1)
#Remove microglia enriched genes from housekeeping gene set
housekeepingGenes <- setdiff(housekeepingGenes,microgliaGenes)

#Map ENSG to Gene Names
load(file='../../Count_Data/geneGtfTableFull.rda')
#only consider protein coding genes
geneIdToName <- hashmap::hashmap(keys = geneGtfTableFull$gene_id,values = geneGtfTableFull$gene_name)

#Get list of genes
rownames(GeneLevelCountMatrix) <- geneIdToName[[rownames(GeneLevelCountMatrix)]]
GeneLevelCV <- apply(GeneLevelCountMatrix,1,function(x) sd(x)/mean(x))
microgliaGenes <- intersect(names(GeneLevelCV),microgliaGenes)
otherGenes <- names(GeneLevelCV)#setdiff(names(GeneLevelCV),c(microgliaGenes,housekeepingGenes))
housekeepingGenes <- intersect(names(GeneLevelCV),housekeepingGenes)


NameToBiotype <- hashmap::hashmap(keys = geneGtfTableFull$gene_name,values = geneGtfTableFull$gene_biotype)
#Get cv for protein coding genes, which are Houskeeping, microglia enriched, and all
proteinCodingMicrogliaCV <- GeneLevelCV[microgliaGenes[which(NameToBiotype[[microgliaGenes]] == 'protein_coding')]]
proteinCodingHKCV <- GeneLevelCV[housekeepingGenes[which(NameToBiotype[[housekeepingGenes]] == 'protein_coding')]]
allGenesCV <- GeneLevelCV[otherGenes[which(NameToBiotype[[otherGenes]]=='protein_coding')]]

#Merge data for each type of genes into data frame
proteinCodingDf <- rbind(data.frame(Type=rep('All\nGenes',length(allGenesCV)),CV=allGenesCV),
                         data.frame(Type=rep('Microglia\nEnriched\nGenes',length(proteinCodingMicrogliaCV)),CV=proteinCodingMicrogliaCV),
                         data.frame(Type=rep('Housekeeping\nGenes',length(proteinCodingHKCV)),CV=proteinCodingHKCV))
library(ggpubr)
library(ggsignif)
wcTest <- compare_means(CV ~ Type, data = proteinCodingDf)
p <- ggplot(proteinCodingDf, aes(x=Type, y=CV)) + 
  geom_boxplot(outlier.shape = NA,fill="#999999") + ylim(0,0.65) + guides(fill=F)  + 
  labs(x='Gene Type',y='Coefficient of Variation') + theme(text = element_text(size = 16),axis.text = element_text(size=16)) + 
  geom_signif(comparisons = list(c("Housekeeping\nGenes", "Microglia\nEnriched\nGenes")), 
              map_signif_level=F,y_position=0.65,test = 'wilcox.test') + geom_hline(yintercept = median(allGenesCV),colour='red',size=3)
ggpubr::ggexport(plotlist = list(ggpubr::ggarrange(p)),
                 filename = '../../FinalFigures/MicrogliaVsHousekeepingCV.pdf',
                 width = 5,height = 4.5)