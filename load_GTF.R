##########################################################################################################
#LoadBiotypeMapping Loads a GTF file, and creates mapping between gene/transcript and their biotype
#ParseBiotypeTable aggregates indivdual biotype categories into 4 broader categories
##########################################################################################################
#Read GTF, create mapping between gene and transcript and exon ID
LoadGTF <- function(full=F){
  library(refGenome)
  library(dplyr)
  GTFPath <- '/local/data/public/zmx21/zmx21_private/GSK/GRCh37_Ensembl75/Homo_sapiens.GRCh37.75.gtf'
  ens <- ensemblGenome()
  read.gtf(ens, GTFPath,useBasedir = F)
  gtfDf <- getGtf(ens)
  if(full){
    return(gtfDf)
  }else{
    return(dplyr::select(gtfDf,transcript_name,gene_name,gene_id,transcript_id,exon_id))
  }
  return(gtfDf)
}

#Reads GTF file, and creates biotype mapping between gene/transcript ID and biotype
#If full=T, the full GTF file columns are returns. Otherwise, only ID and biotype.
LoadBiotypeMapping <- function(full=F){
  gtfDf <- LoadGTF(full=T)
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


#Input should be a data.frame, containg both ID and biotype. 
#Output will be a data.frame, where biotypes are converted to a broader category
ParseBiotypeTable <- function(categoricalCounts){
  type <- ifelse(colnames(categoricalCounts)[1] == 'gene_id','gene','transcript')
  if(type=='gene'){
    #Define which sub-category belong to which full broad category.
    proteinCoding <- c("protein_coding",
                       unique(categoricalCounts$gene_biotype)[grepl("_gene",unique(categoricalCounts$gene_biotype))],
                       'polymorphic_pseudogene')
    pseudoGene <- c('pseudogene',
                    unique(categoricalCounts$gene_biotype)[grepl("_pseudogene",unique(categoricalCounts$gene_biotype))])
    lncRNA <- c('lincRNA','3prime_overlapping_ncrna','antisense','processed_transcript','sense_overlapping','sense_intronic')
    sncRNA <- c('rRNA','misc_RNA','Mt_tRNA','snRNA','Mt_rRNA','miRNA','snoRNA')
    
    #Convert detailed category into broad category, by replacing biotype to the broad category.
    categoricalCounts$gene_biotype[categoricalCounts$gene_biotype %in% proteinCoding] <- 'protein_coding'
    categoricalCounts$gene_biotype[categoricalCounts$gene_biotype %in% pseudoGene] <- 'pseudogene'
    categoricalCounts$gene_biotype[categoricalCounts$gene_biotype %in% lncRNA] <- 'lncRNA'
    categoricalCounts$gene_biotype[categoricalCounts$gene_biotype %in% sncRNA] <- 'sncRNA'
    
    return(categoricalCounts)
  }else{
    #Define which sub-category belong to which full broad category.
    proteinCoding <- c("protein_coding","retained_intron",
                       unique(categoricalCounts$source)[grepl("_gene",unique(categoricalCounts$source))],
                       'polymorphic_pseudogene',unique(categoricalCounts$source)[grepl("_decay",unique(categoricalCounts$source))])
    
    pseudoGene <- c('pseudogene',
                    unique(categoricalCounts$source)[grepl("_pseudogene",unique(categoricalCounts$source))])
    lncRNA <- c('lincRNA','3prime_overlapping_ncrna','antisense','processed_transcript','sense_overlapping','sense_intronic')
    sncRNA <- c('rRNA','misc_RNA','Mt_tRNA','snRNA','Mt_rRNA','miRNA','snoRNA')
    
    #Convert detailed category into broad category, by replacing biotype to the broad category.
    categoricalCounts$source[categoricalCounts$source %in% proteinCoding] <- 'protein_coding'
    categoricalCounts$source[categoricalCounts$source %in% pseudoGene] <- 'pseudogene'
    categoricalCounts$source[categoricalCounts$source %in% lncRNA] <- 'lncRNA'
    categoricalCounts$source[categoricalCounts$source %in% sncRNA] <- 'sncRNA'
    return(categoricalCounts)
  }
}

#Return Categorical counts from a list of genes/transcripts
ExtractBioType <- function(countMatrixInput,table=T){
  load('../Count_Data/gtfTables.rda')
  type <- ifelse(grepl('ENST',rownames(countMatrixInput)[1]),'transcript_id','gene_id')
  allNames <- data.frame(ids = rownames(countMatrixInput),stringsAsFactors = F)
  colnames(allNames) <- type
  categoricalCounts <- as.data.frame(dplyr::left_join(allNames,
                                                      gtfTables[[ifelse(type=='transcript_id','transcriptTable','geneTable')]],by=c(type,type)))
  broadCategoricalCounts <- ParseBiotypeTable(categoricalCounts)
  if(table){
    return(table(broadCategoricalCounts[,ifelse(type=='transcript_id','source','gene_biotype')]))
  }else{
    return(broadCategoricalCounts[,ifelse(type=='transcript_id','source','gene_biotype')])
  }
}
#Return detailed Categorical counts from a list of genes/transcripts
ExtactDetailedBioType <- function(countMatrixInput,table=T){
  load('../Count_Data/gtfTables.rda')
  type <- ifelse(grepl('ENST',rownames(countMatrixInput)[1]),'transcript_id','gene_id')
  allNames <- data.frame(ids = rownames(countMatrixInput),stringsAsFactors = F)
  colnames(allNames) <- type
  categoricalCounts <- as.data.frame(dplyr::left_join(allNames,
                                                      gtfTables[[ifelse(type=='transcript_id','transcriptTable','geneTable')]],by=c(type,type)))
  if(table){
    return(table(categoricalCounts[,ifelse(type=='transcript_id','source','gene_biotype')]))
  }else{
    return(categoricalCounts[,ifelse(type=='transcript_id','source','gene_biotype')])
  }
}
