source('load_GTF.R')
# fullGTF <- LoadGTF(full=T)
library(seqinr)

fasta <- read.fasta('/local/data/public/zmx21/zmx21_private/GSK/GRCh37_Ensembl75/Homo_sapiens.GRCh37.75.cdna.ncrna.fa')
proteinCodingTypes <- c("protein_coding",
                        unique(bioTypeTable$gene_biotype)[grepl("_gene",unique(bioTypeTable$gene_biotype))],
                        'polymorphic_pseudogene')
lncRNATypes <- c('lincRNA')
OthersncRNATypes <- c('miRNA')
pseudoGeneTypes <- c('pseudogene',
                     unique(bioTypeTable$gene_biotype)[grepl("_pseudogene",unique(bioTypeTable$gene_biotype))])

inclTypes <- c(proteinCodingTypes,lncRNATypes,OthersncRNATypes,pseudoGeneTypes)

proteinCodingTranscripts <- dplyr::filter(bioTypeTable$transcriptTable,source %in% inclTypes) %>% dplyr::select('transcript_id') %>% t() %>% as.vector()

filtFASTA <- fasta[names(fasta) %in% proteinCodingTranscripts]
write.fasta(sequences = filtFASTA,names = names(filtFASTA),file.out = '/local/data/public/zmx21/zmx21_private/GSK/GRCh37_Ensembl75/Homo_sapiens.GRCh37.75.cdna.ncrna.filt.fa')



ggplot(data.frame(Samples=rep(colnames(GalatroSamples_Genes),2),MedianTPM <- c(apply(GalatroSamples_Genes,2,median),apply(GalatroFilteredGenes,2,median)),Type=c(rep('all non-coding',ncol(GalatroSamples_Genes)),rep('only miRNA and lincRNA',ncol(GalatroFilteredGenes)))),
       aes(x=Samples,y=MedianTPM,fill=as.factor(Type)))+
  geom_bar(stat="identity",position="dodge",width = 0.7)+
  scale_fill_discrete(name="Type")+
  xlab("Sample")+ylab("Median TPM") + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.6,size =15)) + ggtitle('Galatro Pre and Pos Filtering for non coding RNA')
