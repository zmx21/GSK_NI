##############################################################################################
#Fig. 3.1.0.3
##############################################################################################

load('../../Count_Data/Batch_Corrected/SalmonTPM_Combat_ExpCorrected.rda')
load('../../Count_Data/geneGtfTableFull.rda')

#Only consider protein coding genes
proteinCodingGenes <- geneGtfTableFull %>% dplyr::filter(gene_biotype=='protein_coding') %>% dplyr::select(gene_id) %>% t() %>% as.vector()
batchCorrected <- SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged
proteinCodingBatchCorrected <- batchCorrected[rownames(batchCorrected)%in%proteinCodingGenes,]
cvBatchCorrected <- apply(proteinCodingBatchCorrected,1,function(x) sd(x)/mean(x))

#Low CV and high cv genes
batchCorrectedLowCV <- proteinCodingBatchCorrected[cvBatchCorrected < quantile(cvBatchCorrected,0.1),]
batchCorrectedHighCV <- proteinCodingBatchCorrected[cvBatchCorrected > quantile(cvBatchCorrected,0.9),]

#Extract genes
batchCorrectedLowCVGenes <- data.frame(gene_id=rownames(batchCorrectedLowCV),stringsAsFactors = F) %>% 
  left_join(dplyr::select(geneGtfTableFull,gene_id,gene_name)) %>%
  dplyr::select(gene_name) %>% t() %>% as.vector()
batchCorrectedHighCVGene <- data.frame(gene_id=rownames(batchCorrectedHighCV),stringsAsFactors = F) %>% 
  left_join(dplyr::select(geneGtfTableFull,gene_id,gene_name)) %>%
  dplyr::select(gene_name) %>% t() %>% as.vector()

library(enrichR)
source('../BatchCorrection_and_QC/load_GTF.R')
#convert gene ID to gene names
geneIdToName <- hashmap::hashmap(keys = geneGtfTableFull$gene_id,values = geneGtfTableFull$gene_name)
codingGenes <- hashmap::hashmap(keys = geneGtfTableFull$gene_name,values = ifelse(geneGtfTableFull$gene_biotype=='protein_coding',T,F))
librariesToRun <- c('KEGG_2016','TRANSFAC_and_JASPAR_PWMs','Jensen_DISEASES')

#Run EnrichR
batchCorrectedLowCVKEGG <- enrichr(batchCorrectedLowCVGenes,databases = librariesToRun)
batchCorrectedHighCVKEGG <- enrichr(batchCorrectedHighCVGene,databases = librariesToRun)

library(grid)
library(gtable)
library(gridExtra)

#Low CV Terms
table1 <- tableGrob(batchCorrectedLowCVKEGG$KEGG_2016[1:5,] %>% 
            dplyr::select(KEGG_Term=Term,Overlap=Overlap,Adjusted.P.value) %>% 
              dplyr::mutate(KEGG_Term=sapply(KEGG_Term,function(x) gsub(x=x,replacement = "",pattern = "_Homo.*"))))
title <- textGrob("Bottom 10% of Variable Genes", gp = gpar(fontsize = 16))
padding <- unit(0.5,"line")
table1 <- gtable_add_rows(
  table1, heights = grobHeight(title) + padding, pos = 0
)
table1 <- gtable_add_grob(
  table1, list(title),
  t = 1, l = 1, r = ncol(table1)
)

#High CV Terms
table2 <- tableGrob(batchCorrectedHighCVKEGG$KEGG_2016[1:5,] %>% 
                      dplyr::select(KEGG_Term=Term,Overlap=Overlap,Adjusted.P.value) %>% 
                      dplyr::mutate(KEGG_Term=sapply(KEGG_Term,function(x) gsub(x=x,replacement = "",pattern = "_Homo.*"))))
title <- textGrob("Top 10% of Variable Genes", gp = gpar(fontsize = 16))
padding <- unit(0.5,"line")
table2 <- gtable_add_rows(
  table2, heights = grobHeight(title) + padding, pos = 0
)
table2 <- gtable_add_grob(
  table2, list(title),
  t = 1, l = 1, r = ncol(table2)
)

library(ggpubr)
ggpubr::ggexport(ggpubr::ggarrange(table1,table2,nrow = 2),filename = '../../FinalFigures/Supplementary/KEGGCV.pdf',height = 4.5,width = 5.5)
