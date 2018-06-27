#Filter all samples, according to mapping statistics. 
FilterAllSamplesByMappingRate <- function(plots=T)
{
  library(dplyr)
  library(gridExtra)
  source('pca_analysis.R')
  load(file = '../../Count_Data/gtfTables.rda')
  source('load_GTF.R')
  load('../../Count_Data/Galatro/SalmonTPM_Microglia.rda')
  load('../../Count_Data/Gosselin/SalmonTPM_Microglia.rda')
  load('../../Count_Data/Olah/SalmonTPM_Microglia.rda')
  load('../../Count_Data/Galatro/SalmonTPM_Gene_WholeBrain.rda')
  load('../../Count_Data/Galatro/SalmonTPM_Transcript_WholeBrain.rda')
  
  TPM_Microglia_Gene_Merged <- cbind(GalatroSalmon$geneLevel$abundance,
                                            GosselinSalmon$geneLevel$abundance,
                                            OlahSalmon$geneLevel$abundance)
  TPM_Microglia_Transcript_Merged<- cbind(GalatroSalmon$transcriptLevel$abundance,
                                                 GosselinSalmon$transcriptLevel$abundance,
                                                 OlahSalmon$transcriptLevel$abundance)
  TPM_WholeBrain_Gene <- SalmonTPM_Gene_WholeBrain
  TPM_WholeBrain_Transcript <- SalmonTPM_Transcript_WholeBrain
  
  save(TPM_Microglia_Gene_Merged,file='../../Count_Data/Raw/TPM_Microglia_Gene_Merged.rda')
  save(TPM_Microglia_Transcript_Merged,file='../../Count_Data/Raw/TPM_Microglia_Transcript_Merged.rda')
  save(TPM_WholeBrain_Gene,file='../../Count_Data/Raw/TPM_WholeBrain_Gene.rda')
  save(TPM_WholeBrain_Transcript,file='../../Count_Data/Raw/TPM_WholeBrain_Transcript.rda')
  
  #BoxPlots of All Samples, pre mapping metric filtering. 
  #Samples which should be filtered are highlighted in red
  if(plots){
    library(scales)
    
    #Merge metadata, take intergenic percentage and exonic tags
    microglia_metadata <- CollectMetadata(TPM_Microglia_Gene_Merged) %>% 
      dplyr::mutate(PercIntergenic=ifelse(is.na(intergenicTags) | intergenicTags==0,0,intergenicTags/totalTags)) %>%
      dplyr::select(Sample_Name,ExonicTags=exonTags,PercIntergenic,dataset) %>% dplyr::mutate(ExonicTags=ifelse(is.na(ExonicTags),Inf,ExonicTags)) 
    microglia_gene_df <- stack(as.data.frame(TPM_Microglia_Gene_Merged))
    microglia_gene_df <- dplyr::mutate(microglia_metadata,indShort=ifelse(dataset=="Olah",
                                                                         paste0('OLA_',substr(Sample_Name,11,12),'_',
                                                                                substr(Sample_Name,20,22)),
                                                        ifelse(dataset=='Galatro',
                                                               paste0('GAL',substr(Sample_Name,9,10)),
                                                               paste0('GOS',substr(Sample_Name,9,10))))) %>% 
      dplyr::right_join(microglia_gene_df,by=c("Sample_Name" = "ind"))
    
    tiff(filename = "../../Figures/Mapping_Metric_Filtering/ReadFiltering_TPM_Microglia.tiff", width = 1400, height = 500)
    p1 <- ggplot(microglia_gene_df) + 
      geom_boxplot(aes(x = indShort, y = values,fill=ExonicTags >= 4e+06 & PercIntergenic <= 0.25),outlier.shape = NA) + 
      scale_x_discrete(name = "Sample",labels = wrap_format(5)) +
      scale_y_continuous(name = "TPM",limits = quantile(microglia_gene_df$values, c(0.1, 0.9))) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.6,size =15)) + 
      ggtitle('Micrgolia Gene Level TPM') + 
      labs(fill='Passed Filter')
    grid.arrange(p1)
    dev.off()
    
    brain_metadata <- CollectMetadata(TPM_WholeBrain_Gene) %>% 
      mutate(PercIntergenic=ifelse(is.na(intergenicTags) | intergenicTags==0,0,intergenicTags/totalTags)) %>%
      dplyr::select(Sample_Name,ExonicTags=exonTags,PercIntergenic,dataset) %>% mutate(ExonicTags=ifelse(is.na(ExonicTags),Inf,ExonicTags)) 
    brain_gene_df <- stack(as.data.frame(TPM_WholeBrain_Gene))
    brain_gene_df <- dplyr::mutate(brain_metadata,indShort=ifelse(dataset=="Olah",
                                                                      paste0('OLA_',substr(Sample_Name,11,12),'_',
                                                                             substr(Sample_Name,20,22)),
                                                                      ifelse(dataset=='Galatro',
                                                                             paste0('GAL',substr(Sample_Name,9,10)),
                                                                             paste0('GOS',substr(Sample_Name,9,10))))) %>% 
      dplyr::right_join(brain_gene_df,by=c("Sample_Name" = "ind"))
    tiff(height = 800,width=1000,filename = "../../Figures/Mapping_Metric_Filtering/ReadFiltering_TPM_Brain.tiff")
    p2 <- ggplot(brain_gene_df) + 
      geom_boxplot(aes(x = indShort, y = values,fill=ExonicTags >= 4e+06 & PercIntergenic <= 0.25),outlier.shape = NA) + 
      scale_x_discrete(name = "Sample",labels = wrap_format(5)) +
      scale_y_continuous(name = "TPM",limits = c(0,4)) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.6,size =15)) + 
      ggtitle('Whole Brain Gene Level TPM') + 
      labs(fill='Passed Filter')
    grid.arrange(p2)
    dev.off()
    MappingFilteringEnv <- environment()
  }
  #Remove samples from Galatro which have bad mapping. 
  BadMappingGalatro <- read.table(file='/local/data/public/zmx21/zmx21_private/GSK/Galatro/STAR_aligned_merged/multiqc_data/multiqc_rseqc_read_distribution.txt',
                                  header = T,stringsAsFactors = F) %>% dplyr::select('Sample','cds_exons_tag_count','total_tags','other_intergenic_tag_count') %>% 
    dplyr::mutate(Sample_Name = sapply(.$Sample,function(x) unlist(strsplit(x,'[.]'))[1]),
                  PercExonic = round(as.numeric(cds_exons_tag_count/total_tags),2),
                  ExonicTags = round(as.numeric(cds_exons_tag_count),2),
                  PercIntergenic = round(as.numeric((other_intergenic_tag_count)/total_tags),2)) %>% 
    dplyr::filter(ExonicTags < 4e6 | PercIntergenic > 0.25) %>% dplyr::select(Sample_Name) %>% t() %>% as.vector() #Filter by mapping rate 
  BadMappingGalatroBrain <- read.table(file='/local/data/public/zmx21/zmx21_private/GSK/Galatro_Brain/STAR_aligned_whole_brain/multiqc_data/multiqc_rseqc_read_distribution.txt',
                                       header = T,stringsAsFactors = F) %>% dplyr::select('Sample','cds_exons_tag_count','total_tags','other_intergenic_tag_count') %>% 
    dplyr::mutate(Sample_Name = sapply(.$Sample,function(x) unlist(strsplit(x,'[.]'))[1]),
                  PercExonic = round(as.numeric(cds_exons_tag_count/total_tags),2),
                  ExonicTags = round(as.numeric(cds_exons_tag_count),2),
                  PercIntergenic = round(as.numeric((other_intergenic_tag_count)/total_tags),2)) %>% 
    dplyr::filter(ExonicTags < 4e6 | PercIntergenic > 0.25) %>% dplyr::select(Sample_Name) %>% t() %>% as.vector() #Filter by mapping rate 
  
  #Remove samples from Olah which have bad mapping.
  BadMappingOlah <- read.table(file='/local/data/public/zmx21/zmx21_private/GSK/Olah/QC_Reports/star_multiqc_data/multiqc_rseqc_read_distribution.txt',
                               header = T,stringsAsFactors = F) %>% dplyr::select('Sample_Name'='Sample','cds_exons_tag_count','total_tags') %>% 
    dplyr::mutate(PercExonic = round(as.numeric(cds_exons_tag_count/total_tags),2),
                  ExonicTags = round(as.numeric(cds_exons_tag_count),2)) %>% 
    dplyr::filter(ExonicTags < 4e6) %>% dplyr::select(Sample_Name) %>% t() %>% as.vector() #Filter by mapping rate 
  
  
  #Remove bad mapping samples
  TPM_Microglia_Gene_Merged <- TPM_Microglia_Gene_Merged[,!colnames(TPM_Microglia_Gene_Merged)%in% c(BadMappingGalatro,BadMappingOlah)]
  TPM_Microglia_Transcript_Merged <- TPM_Microglia_Transcript_Merged[,!colnames(TPM_Microglia_Transcript_Merged)%in% c(BadMappingGalatro,BadMappingOlah)]
  TPM_WholeBrain_Gene <- TPM_WholeBrain_Gene[,!colnames(TPM_WholeBrain_Gene)%in% c(BadMappingGalatroBrain)]
  TPM_WholeBrain_Transcript <- TPM_WholeBrain_Transcript[,!colnames(TPM_WholeBrain_Transcript)%in% c(BadMappingGalatroBrain)]
  
  #Remove rRNA,tRNA,scRMA
  genesToFilter <- gtfTables$geneTable %>% dplyr::filter(gene_biotype%in%c('rRNA','Mt_rRNA','Mt_tRNA','tRNA','tRNA_pseudogene','scRNA')) %>% dplyr::select(gene_id) %>% t() %>% as.vector()
  transcriptsToFilter <- gtfTables$transcriptTable %>% dplyr::filter(source%in%c('rRNA','Mt_rRNA','Mt_tRNA','tRNA','tRNA_pseudogene','scRNA')) %>% dplyr::select(transcript_id) %>% t() %>% as.vector()
  TPM_Microglia_Gene_Merged <- TPM_Microglia_Gene_Merged[!rownames(TPM_Microglia_Gene_Merged)%in%genesToFilter,]
  TPM_Microglia_Transcript_Merged <- TPM_Microglia_Transcript_Merged[!rownames(TPM_Microglia_Transcript_Merged)%in%transcriptsToFilter,]
  TPM_WholeBrain_Gene <- TPM_WholeBrain_Gene[!rownames(TPM_WholeBrain_Gene)%in%genesToFilter,]
  TPM_WholeBrain_Transcript <- TPM_WholeBrain_Transcript[!rownames(TPM_WholeBrain_Transcript)%in%transcriptsToFilter,]
  
  
  save(TPM_Microglia_Gene_Merged,file='../../Count_Data/Read_Filtered/TPM_Microglia_Gene_Merged_ReadFilt.rda')
  save(TPM_Microglia_Transcript_Merged,file='../../Count_Data/Read_Filtered/TPM_Microglia_Transcript_Merged_ReadFilt.rda')
  save(TPM_WholeBrain_Gene,file='../../Count_Data/Read_Filtered/TPM_WholeBrain_Gene.rda')
  save(TPM_WholeBrain_Transcript,file='../../Count_Data/Read_Filtered/TPM_WholeBrain_Transcript')
  # save(list = ls(environment()),file='../../CodeImages/MappingFiltering.RData')
  
}
