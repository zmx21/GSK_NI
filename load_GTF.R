LoadBiotypeMapping <- function(full=F){
  library(refGenome)
  library(dplyr)
  GTFPath <- '/local/data/public/zmx21/zmx21_private/GSK/GRCh37_Ensembl75/Homo_sapiens.GRCh37.75.gtf'
  ens <- ensemblGenome()
  read.gtf(ens, GTFPath,useBasedir = F)
  gtfDf <- getGtf(ens)
  if(full){
    geneTable <- as.data.frame(gtfDf %>%
                                 dplyr::group_by(gene_id, gene_biotype) %>% 
                                 dplyr::filter(row_number() == 1))
    allBiotypes <- unique(geneTable$gene_biotype)
    transcriptTable <- dplyr::filter(gtfDf,feature=='transcript')
  }else{
    geneTable <- as.data.frame(dplyr::select(gtfDf,"gene_id","gene_biotype") %>%
                                 dplyr::group_by(gene_id, gene_biotype) %>% 
                                 dplyr::filter(row_number() == 1))
    allBiotypes <- unique(geneTable$gene_biotype)
    transcriptTable <- dplyr::filter(gtfDf,feature=='transcript') %>% dplyr::select("transcript_id","source")
  }
  return(list(geneTable = geneTable,transcriptTable = transcriptTable))
}
ParseBiotypeTable <- function(categoricalCounts){
  type <- ifelse(colnames(categoricalCounts)[1] == 'gene_id','gene','transcript')
  if(type=='gene'){
    proteinCoding <- c("protein_coding",
                       unique(categoricalCounts$gene_biotype)[grepl("_gene",unique(categoricalCounts$gene_biotype))],
                       'polymorphic_pseudogene')
    pseudoGene <- c('pseudogene',
                    unique(categoricalCounts$gene_biotype)[grepl("_pseudogene",unique(categoricalCounts$gene_biotype))])
    lncRNA <- c('lincRNA','3prime_overlapping_ncrna','antisense','processed_transcript','sense_overlapping','sense_intronic')
    sncRNA <- c('rRNA','misc_RNA','Mt_tRNA','snRNA','Mt_rRNA','miRNA','snoRNA')
    
    categoricalCounts$gene_biotype[categoricalCounts$gene_biotype %in% proteinCoding] <- 'protein_coding'
    categoricalCounts$gene_biotype[categoricalCounts$gene_biotype %in% pseudoGene] <- 'pseudogene'
    categoricalCounts$gene_biotype[categoricalCounts$gene_biotype %in% lncRNA] <- 'lncRNA'
    categoricalCounts$gene_biotype[categoricalCounts$gene_biotype %in% sncRNA] <- 'sncRNA'
    
    return(categoricalCounts)
  }else{
    proteinCoding <- c("protein_coding","retained_intron",
                       unique(categoricalCounts$source)[grepl("_gene",unique(categoricalCounts$source))],
                       'polymorphic_pseudogene',unique(categoricalCounts$source)[grepl("_decay",unique(categoricalCounts$source))])
    
    pseudoGene <- c('pseudogene',
                    unique(categoricalCounts$source)[grepl("_pseudogene",unique(categoricalCounts$source))])
    lncRNA <- c('lincRNA','3prime_overlapping_ncrna','antisense','processed_transcript','sense_overlapping','sense_intronic')
    sncRNA <- c('rRNA','misc_RNA','Mt_tRNA','snRNA','Mt_rRNA','miRNA','snoRNA')
    categoricalCounts$source[categoricalCounts$source %in% proteinCoding] <- 'protein_coding'
    categoricalCounts$source[categoricalCounts$source %in% pseudoGene] <- 'pseudogene'
    categoricalCounts$source[categoricalCounts$source %in% lncRNA] <- 'lncRNA'
    categoricalCounts$source[categoricalCounts$source %in% sncRNA] <- 'sncRNA'
    return(categoricalCounts)
  }
}
