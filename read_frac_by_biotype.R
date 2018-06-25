source('load_GTF.R')
source('import_salmon.R')
fullGTF <- LoadGTF(full=T)
# allPaths <- c('/local/data/public/zmx21/zmx21_private/GSK/Galatro/Salmon_aligned_merged/',
#               '/local/data/public/zmx21/zmx21_private/GSK/Gosselin/salmon/',
#               '/local/data/public/zmx21/zmx21_private/GSK/Olah/salmon/')
# tx2gene <- ImportTx2gene()
# allCounts <- lapply(allPaths,function(x)ImportSalmonCounts(x,tx2gene)$geneLevel$counts)
# allCounts <- do.call(cbind,allCounts)

#Remove Bad Samples
#Remove samples from Galatro which have bad mapping. 
# BadMappingGalatro <- read.table(file='/local/data/public/zmx21/zmx21_private/GSK/Galatro/STAR_aligned_merged/multiqc_data/multiqc_rseqc_read_distribution.txt',
#                                 header = T,stringsAsFactors = F) %>% dplyr::select('Sample','cds_exons_tag_count','total_tags','other_intergenic_tag_count') %>% 
#   dplyr::mutate(Sample_Name = sapply(.$Sample,function(x) unlist(strsplit(x,'[.]'))[1]),
#                 PercExonic = round(as.numeric(cds_exons_tag_count/total_tags),2),
#                 ExonicTags = round(as.numeric(cds_exons_tag_count),2),
#                 PercIntergenic = round(as.numeric((other_intergenic_tag_count)/total_tags),2)) %>% 
#   dplyr::filter(ExonicTags < 4e6 | PercIntergenic > 0.25) %>% dplyr::select(Sample_Name) %>% t() %>% as.vector() #Filter by mapping rate 
# #Remove samples from Olah which have bad mapping.
# BadMappingOlah <- read.table(file='/local/data/public/zmx21/zmx21_private/GSK/Olah/QC_Reports/star_multiqc_data/multiqc_rseqc_read_distribution.txt',
#                              header = T,stringsAsFactors = F) %>% dplyr::select('Sample_Name'='Sample','cds_exons_tag_count','total_tags') %>% 
#   dplyr::mutate(PercExonic = round(as.numeric(cds_exons_tag_count/total_tags),2),
#                 ExonicTags = round(as.numeric(cds_exons_tag_count),2)) %>% 
#   dplyr::filter(ExonicTags < 4e6) %>% dplyr::select(Sample_Name) %>% t() %>% as.vector() #Filter by mapping rate 
# 
# 
# #Remove bad mapping samples
# allCounts <- allCounts[,!colnames(allCounts)%in% c(BadMappingGalatro,BadMappingOlah)]

bioTypeTable <- as.data.frame(fullGTF %>% 
                        dplyr::group_by(gene_id, gene_biotype) %>% 
                        dplyr::filter(row_number() == 1)) %>% dplyr::select(gene_id,gene_biotype)

proteinCodingTypes <- c("protein_coding",
                   unique(bioTypeTable$gene_biotype)[grepl("_gene",unique(bioTypeTable$gene_biotype))],
                   'polymorphic_pseudogene')
lncRNATypes <- c('lincRNA','3prime_overlapping_ncrna','antisense','processed_transcript','sense_overlapping','sense_intronic')
OthersncRNATypes <- c('misc_RNA','Mt_tRNA','snRNA','miRNA','snoRNA')
pseudoGeneTypes <- c('pseudogene',
                unique(bioTypeTable$gene_biotype)[grepl("_pseudogene",unique(bioTypeTable$gene_biotype))])

rRNAGenes <- bioTypeTable %>% dplyr::filter(gene_biotype%in%c('rRNA','Mt_rRNA')) %>% select(gene_id) %>% t() %>% as.vector()
proteinCodingGenes <- bioTypeTable %>% dplyr::filter(gene_biotype%in%proteinCodingTypes)%>% select(gene_id) %>% t() %>% as.vector()
pseudoGene <- bioTypeTable %>% dplyr::filter(gene_biotype%in%pseudoGeneTypes)%>% select(gene_id) %>% t() %>% as.vector()
lncRNA <- bioTypeTable%>% dplyr::filter(gene_biotype%in%lncRNATypes)%>% select(gene_id) %>% t() %>% as.vector()
OthersncRNA <- bioTypeTable%>% dplyr::filter(gene_biotype%in%OthersncRNATypes)%>% select(gene_id) %>% t() %>% as.vector()


CalcCountFraction <- function(geneNames,countMatrix){
  indices <- which(rownames(countMatrix) %in% geneNames)
  frac <- apply(countMatrix,2,function(x) sum(x[indices])/sum(x))
  return(frac)
}

BiotypeStats <- data.frame(Sample = rep(colnames(SalmonTPM_Gene_Merged),5),
                           Biotype = c(rep('rRNA',ncol(SalmonTPM_Gene_Merged)),
                                       rep('Protein_Coding',ncol(SalmonTPM_Gene_Merged)),
                                       rep('PseudoGene',ncol(SalmonTPM_Gene_Merged)),
                                       rep('lncRNA',ncol(SalmonTPM_Gene_Merged)),
                                       rep('OthersncRNA',ncol(SalmonTPM_Gene_Merged))),
                           Fraction = do.call(c,lapply(list(rRNAGenes,proteinCodingGenes,pseudoGene,lncRNA,OthersncRNA),function(x) CalcCountFraction(x,SalmonTPM_Gene_Merged))))
BiotypeStats$Dataset = ifelse(substr(BiotypeStats$Sample,1,3) == 'GSM','Galatro',ifelse(substr(BiotypeStats$Sample,1,3) == 'SRR','Gosselin','Olah'))
BiotypeStats <- BiotypeStats %>% left_join(CollectMetadata(SalmonTPM_Gene_Merged),by=c('Sample'='Sample_Name'))


ggplot(data=BiotypeStats, aes(x=Sample, y=Fraction, fill=Biotype)) +
  geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.6,size =15))

nonProteinCoding <- sapply(unique(BiotypeStats$Sample),function(x) sum(dplyr::filter(BiotypeStats,Sample==x,Biotype=='rRNA')$Fraction))
nonProteinCoding <- cbind(data.frame(Sample=unique(BiotypeStats$Sample)),data.frame(Fraction = nonProteinCoding)) %>% 
  left_join(CollectMetadata(SalmonTPM_Gene_Merged),by=c('Sample'='Sample_Name')) 
nonProteinCoding$sdTPM <- sapply(nonProteinCoding$Sample,function(x) sd(SalmonTPM_Gene_Merged[,x]))
nonProteinCoding$medianTPM <- sapply(nonProteinCoding$Sample,function(x) median(SalmonTPM_Gene_Merged[,x]))

nonProteinCoding$Dataset = ifelse(substr(nonProteinCoding$Sample,1,3) == 'GSM','Galatro',ifelse(substr(nonProteinCoding$Sample,1,3) == 'SRR','Gosselin','Olah'))


ggarrange(ggplot(nonProteinCoding) + aes(x=Fraction,y=medianTPM,color=Dataset) + geom_point() +
  labs(x='Fraction of TPM Abundance which are Non-Protein Coding',y='Median TPM') + ggtitle('Median TPM vs Fraction of Non-Protein Coding'),
ggplot(nonProteinCoding) + aes(x=mappingRate,y=medianTPM,color=Dataset) + geom_point() +
  labs(x='Salmon Mapping Rate',y='Median TPM') + ggtitle('Median TPM vs Mapping Rate'),ncol=3)
# 
# ggplot(data = nonProteinCoding,aes(x=Fraction,y=mappingRate)) + geom_point(aes(colour=medianTPM))

ggplot(data = nonProteinCoding,aes(x=Fraction,y=mappingRate)) + geom_point(aes(shape=Dataset,colour=medianTPM),size=3) +
  labs(x='Fraction of rRNA and sncRNA Genes',y='Mapping Rate',colour='Median TPM')
ggplot(data=nonProteinCoding,aes(x=sdTPM,y=medianTPM)) + geom_point(aes(shape=Dataset),size=3) + labs(x='TPM Interquartile Range',y='TPM Median') + ggtitle('Median vs IQR of TPM')
