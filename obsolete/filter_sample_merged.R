# pseudoGeneTypes <- c('pseudogene',
#                      unique(gtfTables$geneTable$gene_biotype)[grepl("_pseudogene",unique(gtfTables$geneTable$gene_biotype))])
# pseudoGeneTypesTranscripts <- c('pseudogene',
#                                 unique(gtfTables$transcriptTable$source)[grepl("_pseudogene",unique(gtfTables$transcriptTable$source))])
# pseudoGenes <- gtfTables$geneTable %>% dplyr::filter(gtfTables$geneTable$gene_biotype%in%pseudoGeneTypes)%>% select(gene_id) %>% t() %>% as.vector()
#pseudoGenesTranscript <- gtfTables$transcriptTable %>% dplyr::filter(gtfTables$transcriptTable$source%in%pseudoGeneTypesTranscripts)%>% select(transcript_id) %>% t() %>% as.vector()
#SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged <- SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged[!rownames(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged)%in%pseudoGenes,]
#SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged <- SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged[!rownames(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged)%in%pseudoGenesTranscript,]

#Filter accroding to CV
meanGene_Microglia <- apply(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged,1,function(x) mean(x))
meanTranscript_Microglia <- apply(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged,1,function(x) mean(x))
cvGene_Microglia <- apply(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged,1,function(x) sd(x)/mean(x))
cvTranscript_Microglia <- apply(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged,1,function(x) sd(x)/mean(x))

proteinCodingGenes <- gtfTables$geneTable %>% dplyr::filter(gene_biotype=='protein_coding') %>% dplyr::select(gene_id) %>% t() %>% as.vector()
proteinCodingTranscripts <- gtfTables$transcriptTable %>% dplyr::filter(source=='protein_coding') %>% dplyr::select(transcript_id) %>% t() %>% as.vector()

proteinCodingGeneBatchCorrected <- SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged[rownames(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged)%in%proteinCodingGenes,]
proteinCodingTranscriptBatchCorrected <- SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged[rownames(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged)%in%proteinCodingTranscripts,]

meanGene_proteinCoding_Microglia <- apply(proteinCodingGeneBatchCorrected,1,function(x) mean(x))
meanTranscript_proteinCoding_Microglia <- apply(proteinCodingTranscriptBatchCorrected,1,function(x) mean(x))
cvGene_proteinCoding_Microglia <- apply(proteinCodingGeneBatchCorrected,1,function(x) sd(x)/mean(x))
cvTranscript_proteinCoding_Microglia <- apply(proteinCodingTranscriptBatchCorrected,1,function(x) sd(x)/mean(x))

par(mfrow=c(2,2))
hist(cvGene_Microglia,breaks=100,main='CV Gene Micrgolia - All',xlab='CV')
lines(c(quantile(cvGene_Microglia,0.75),quantile(cvGene_Microglia,0.75)),y=c(0,100000),col='red',lwd=3)
hist(cvTranscript_Microglia,breaks=100,main='CV Transcript Microglia - All',xlab='CV')
lines(c(quantile(cvTranscript_Microglia,0.75),quantile(cvTranscript_Microglia,0.75)),y=c(0,100000),col='red',lwd=3)
hist(cvGene_proteinCoding_Microglia,breaks=100,main='CV Gene Micrgolia - Protein Coding',xlab='CV')
lines(c(quantile(cvGene_proteinCoding_Microglia,0.75),quantile(cvGene_proteinCoding_Microglia,0.75)),y=c(0,100000),col='red',lwd=3)
hist(cvTranscript_Microglia,breaks=100,main='CV Transcript Microglia - Protein Coding',xlab='CV')
lines(c(quantile(cvTranscript_proteinCoding_Microglia,0.75),quantile(cvTranscript_proteinCoding_Microglia,0.75)),y=c(0,100000),col='red',lwd=3)

geneLevelFiltCounts_Microglia <- FilterCountMatrix(countMatrix = SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged,
                                                                                          meanAbundances = meanGene_Microglia,
                                                                                          cv = cvGene_Microglia,
                                                                                          cvCutOffAbsolute = quantile(cvGene_Microglia,0.5),
                                                                                          meanCutOffAbsolute = -Inf,
                                                                                          mappingTable = gtfTables$geneTable) %>%
                                                                                          {ParseBiotypeTable(.$categoricalCounts)} %>%
                                                                                          {table(.$gene_biotype)}
transcriptLevelFiltCounts_Microglia <- FilterCountMatrix(countMatrix = SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged,
                                                                                          meanAbundances = meanTranscript_Microglia,
                                                                                          cv = cvTranscript_Microglia,
                                                                                          cvCutOffAbsolute = quantile(cvTranscript_Microglia,0.5),
                                                                                          meanCutOffAbsolute = -Inf,
                                                                                          mappingTable = gtfTables$transcriptTable) %>%
                                                                                          {ParseBiotypeTable(.$categoricalCounts)} %>%
                                                                                          {table(.$source)} #%>% {do.call(rbind,.)} %>% {cbind(geneLevelCutOff,.)} %>% 
  #dplyr::select(-TPM) %>% dplyr::mutate(Total_Transcripts = lncRNA+protein_coding+pseudogene+sncRNA) %>%dplyr::mutate(Frac_Filtered = 1-Total_Transcripts/nrow(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged)) 
library(grid)
library(gridExtra)
allGenesNetwork <- data.frame(cbind(rbind(geneLevelFiltCounts_Microglia,transcriptLevelFiltCounts_Microglia)),data.frame(Total = c(sum(geneLevelFiltCounts_Microglia),sum(transcriptLevelFiltCounts_Microglia))))
rownames(allGenesNetwork) <- c('Gene','Transcript')
grid.table(allGenesNetwork)

#Brain Data
SalmonTPM_Gene_WholeBrain <- log2(SalmonTPM_Gene_WholeBrain+1)
SalmonTPM_Transcript_WholeBrain <- log2(SalmonTPM_Transcript_WholeBrain+1)
brainProteinCodingGenes <- SalmonTPM_Gene_WholeBrain[rownames(SalmonTPM_Gene_WholeBrain)%in%proteinCodingGenes,]
brainProteinCodingTranscripts <- SalmonTPM_Transcript_WholeBrain[rownames(SalmonTPM_Transcript_WholeBrain)%in%proteinCodingTranscripts,]
meanGene_Brain <- apply(SalmonTPM_Gene_WholeBrain,1,function(x) mean(x))
meanTranscript_Brain <- apply(SalmonTPM_Transcript_WholeBrain,1,function(x) mean(x))
cvGene_Brain <- apply(SalmonTPM_Gene_WholeBrain,1,function(x) sd(x)/mean(x))
cvTranscript_Brain <- apply(SalmonTPM_Transcript_WholeBrain,1,function(x) sd(x)/mean(x))

meanGene_ProteinCoding_Brain <- apply(brainProteinCodingGenes,1,function(x) mean(x))
meanTranscript_ProteinCoding_Brain <- apply(brainProteinCodingTranscripts,1,function(x) mean(x))
cvGene_ProteinCoding_Brain <- apply(brainProteinCodingGenes,1,function(x) sd(x)/mean(x))
cvTranscript_ProteinCoding_Brain <- apply(brainProteinCodingTranscripts,1,function(x) sd(x)/mean(x))


par(mfrow=c(2,2))
hist(cvGene_Brain,breaks=100,main='CV GeneBrain - All',xlab='CV')
lines(rep(quantile(cvGene_Brain,0.75,na.rm = T),2),y=c(0,100000),col='red',lwd=3)
hist(cvTranscript_Brain,breaks=100,main='CV Transcript Brain - All',xlab='CV')
lines(rep(quantile(cvTranscript_Brain,0.75,na.rm = T),2),y=c(0,100000),col='red',lwd=3)
hist(cvGene_ProteinCoding_Brain,breaks=100,main='CV GeneBrain - Protein Coding',xlab='CV')
lines(rep(quantile(cvGene_ProteinCoding_Brain,0.75,na.rm = T),2),y=c(0,100000),col='red',lwd=3)
hist(cvTranscript_ProteinCoding_Brain,breaks=100,main='CV Transcript Brain - Protein Coding',xlab='CV')
lines(rep(quantile(cvTranscript_ProteinCoding_Brain,0.75,na.rm = T),2),y=c(0,100000),col='red',lwd=3)

geneLevelFiltCounts_Brain <- FilterCountMatrix(countMatrix = SalmonTPM_Gene_WholeBrain,
                                                   meanAbundances = meanGene_Brain,
                                                   cv = cvGene_Brain,
                                                   cvCutOffAbsolute = quantile(cvGene_Brain,0.75),
                                                   meanCutOffAbsolute = -Inf,
                                                   mappingTable = gtfTables$geneTable) %>%
                                                   {ParseBiotypeTable(.$categoricalCounts)} %>%
                                                   {table(.$gene_biotype)}
transcriptLevelFiltCounts_Brain <- FilterCountMatrix(countMatrix = SalmonTPM_Transcript_WholeBrain,
                                                         meanAbundances = meanTranscript_Brain,
                                                         cv = cvTranscript_Brain,
                                                         cvCutOffAbsolute = quantile(cvTranscript_Brain,0.75),
                                                         meanCutOffAbsolute = -Inf,
                                                         mappingTable = gtfTables$transcriptTable) %>%
                                                         {ParseBiotypeTable(.$categoricalCounts)} %>%
                                                         {table(.$source)} #%>% {do.call(rbind,.)} %>% {cbind(geneLevelCutOff,.)} %>% 
#dplyr::select(-TPM) %>% dplyr::mutate(Total_Transcripts = lncRNA+protein_coding+pseudogene+sncRNA) %>%dplyr::mutate(Frac_Filtered = 1-Total_Transcripts/nrow(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged)) 
library(grid)
library(gridExtra)
allGenesNetwork <- data.frame(cbind(rbind(geneLevelFiltCounts_Brain,transcriptLevelFiltCounts_Brain)),data.frame(Total = c(sum(geneLevelFiltCounts_Brain),sum(transcriptLevelFiltCounts_Brain))))
rownames(allGenesNetwork) <- c('Gene','Transcript')
grid.table(allGenesNetwork)


# geneLevelCutOffTPM <- quantile(meanAbundancesGene,prob=seq(0.02,0.2,0.02)) %>% round(2)
# geneLevelCutOffCV <- 0.3
# geneLevelCutOff <- expand.grid(TPM=geneLevelCutOffTPM,CV=geneLevelCutOffCV)
# 
# geneLevelFiltCounts <- lapply(1:nrow(geneLevelCutOff), function(i) FilterCountMatrix(countMatrix = SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged,
#                                                                                           meanAbundances = meanAbundancesGene,
#                                                                                           cv = cvGene,
#                                                                                           cvCutOffAbsolute = geneLevelCutOff[i,1],
#                                                                                           meanCutOffAbsolute = geneLevelCutOff[i,2],
#                                                                                           mappingTable = gtfTables$geneTable) %>%
#                                                                                           {ParseBiotypeTable(.$categoricalCounts)} %>%
#                                                                                           {table(.$gene_biotype)}) %>% {do.call(rbind,.)} %>% {cbind(geneLevelCutOff,.)}
# transcriptLevelFiltCounts <- lapply(1:nrow(geneLevelFiltCounts), function(i) FilterCountMatrix(countMatrix = SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged,
#                                                                                           meanAbundances = meanAbundancesTranscript,
#                                                                                           cv = cvTranscript,
#                                                                                           cvCutOffAbsolute = geneLevelCutOff[i,1],
#                                                                                           meanCutOffAbsolute = geneLevelCutOff[i,2],
#                                                                                           mappingTable = gtfTables$transcriptTable) %>%
#                                                                                           {ParseBiotypeTable(.$categoricalCounts)} %>%
#                                                                                           {table(.$source)}) %>% {do.call(rbind,.)} %>% {cbind(geneLevelCutOff,.)}
# geneLevelFiltCounts$TPM_Percentile <- seq(2,20,2)
# geneLevelFiltCounts <- geneLevelFiltCounts[c(7,1:6)]
# 
# transcriptLevelFiltCounts$TPM_Percentile <- seq(2,20,2)
# transcriptLevelFiltCounts <- transcriptLevelFiltCounts[c(7,1:6)]

# library(gridExtra)
# grid.table(round(geneLevelFiltCounts,3))
# grid.table(round(transcriptLevelFiltCounts,3))


#Bias in read Distribution
# MappingRatePlot <- ggplot() + aes(x=CollectMetadata(SalmonTPM_Gene_Merged,full=F)$mappingRate,y=apply(SalmonTPM_Gene_Merged,2,median),
#                                   color=CollectMetadata(SalmonTPM_Gene_Merged,full=F)$dataset) + geom_point() + 
#                               labs(color='Dataset',x = 'Mapping Rate',y='Median TPM') + ggtitle('Median TPM VS Mapping Rate ')
# 
# TotalReadsPlot <- ggplot() + aes(x=CollectMetadata(SalmonTPM_Gene_Merged,full=F)$numReads,y=apply(SalmonTPM_Gene_Merged,2,median),
#                                   color=CollectMetadata(SalmonTPM_Gene_Merged,full=F)$dataset) + geom_point() + 
#   labs(color='Dataset',x = 'Num Reads',y='Median TPM') + ggtitle('Median TPM VS Num Reads ')
# 
# IntronicRatePlot <- ggplot() + aes(x=CollectMetadata(SalmonTPM_Gene_Merged,full=F)$intronTags/CollectMetadata(SalmonTPM_Gene_Merged,full=F)$totalTags,y=apply(SalmonTPM_Gene_Merged,2,median),
#                                   color=CollectMetadata(SalmonTPM_Gene_Merged,full=F)$dataset) + geom_point() + 
#   labs(color='Dataset',x = 'Intronic Tags',y='Median TPM') + ggtitle('Median TPM VS Intronic Reads ')
# ggarrange(plots=list(MappingRatePlot,TotalReadsPlot,IntronicRatePlot),ncol=3)


# #Remove genes which have 0 expression across all Samples
# SalmonTPM_Gene_Merged <- SalmonTPM_Gene_Merged[!apply(SalmonTPM_Gene_Merged,1,function(x) all(x==0)),]
# SalmonTPM_Transcript_Merged <- SalmonTPM_Transcript_Merged[!apply(SalmonTPM_Transcript_Merged,1,function(x) all(x==0)),]
# 
# SalmonTPM_Combat <- RunCombat(SalmonTPM_Gene_Merged,SalmonTPM_Transcript_Merged,Samples = c('Galatro','Gosselin'),expType = F,full=F)
# SalmonTPM_Combat_ExpCorrected <- RunCombat(SalmonTPM_Gene_Merged,SalmonTPM_Transcript_Merged,Samples = c('Galatro','Gosselin'),expType = T,full=F)
# 
# save(SalmonTPM_Combat,file='../Count_Data/Batch_Corrected/Salmon_TPM_Combat_Par.rda')
# save(SalmonTPM_Combat_ExpCorrected,file='../Count_Data/Batch_Corrected/SalmonTPM_Combat_Par_ExpCorrected.rda')
# 
# #Look at Biotypes under different cutoffs
# meanAbundancesGene <- apply(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged,1,function(x) mean(x))
# cvGene <- apply(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged,1,function(x) sd(x) / mean(x))
# #Get Transcript Level Stats
# meanAbundancesTranscript <- apply(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged,1,function(x) mean(x))
# cvTranscript <- apply(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged,1,function(x) sd(x) / mean(x))
# 
# geneLevelCutOffTPM <- quantile(meanAbundancesGene,prob=seq(0.1,0.3,0.1)) %>% round(2)
# geneLevelCutOffCV <- quantile(cvGene,prob=seq(0.1,0.3,0.1)) %>% round(2)
# geneLevelCutOff <- expand.grid(TPM=geneLevelCutOffTPM,CV=geneLevelCutOffCV)
# 
# geneLevelFiltCounts <- lapply(1:nrow(geneLevelCutOff), function(i) FilterCountMatrix(countMatrix = SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged,
#                                                                                           meanAbundances = meanAbundancesGene,
#                                                                                           cv = cvGene,
#                                                                                           cvCutOffAbsolute = geneLevelCutOff[i,1],
#                                                                                           meanCutOffAbsolute = geneLevelCutOff[i,2],
#                                                                                           mappingTable = gtfTables$geneTable) %>%
#                                                                                           {ParseBiotypeTable(.$categoricalCounts)} %>%
#                                                                                           {table(.$gene_biotype)}) %>% {do.call(rbind,.)} %>% {cbind(geneLevelCutOff,.)}
# transcriptLevelFiltCounts <- lapply(1:nrow(geneLevelFiltCounts), function(i) FilterCountMatrix(countMatrix = SalmonTPM_Combat_ExpCorrected$SalmonTPM_Transcript_Combat_Merged,
#                                                                                           meanAbundances = meanAbundancesTranscript,
#                                                                                           cv = cvTranscript,
#                                                                                           cvCutOffAbsolute = geneLevelCutOff[i,1],
#                                                                                           meanCutOffAbsolute = geneLevelCutOff[i,2],
#                                                                                           mappingTable = gtfTables$transcriptTable) %>%
#                                                                                           {ParseBiotypeTable(.$categoricalCounts)} %>%
#                                                                                           {table(.$source)}) %>% {do.call(rbind,.)} %>% {cbind(geneLevelCutOff,.)}
# library(gridExtra)
# grid.table(round(geneLevelFiltCounts,3))
# grid.table(round(transcriptLevelFiltCounts,3))
# 
# 
# #Plot resulting boxplot before batch correction according to mapping rate. 
# TPMVsReads <- ggplot(stack(as.data.frame(SalmonTPM_Gene_Merged)) %>% 
#                dplyr::left_join(CollectMetadata(SalmonTPM_Gene_Merged,full=F),by=c('ind' = 'Sample_Name'))) + 
#   geom_boxplot(aes(x = ind, y = values,fill=numReads),outlier.shape = NA) + scale_x_discrete(name = "Sample") +
#   scale_y_continuous(name = "TPM",limits = c(0,4)) + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.6,size =15)) + scale_colour_gradientn(colours = rainbow(7))
# 
# source('pca_analysis.R')
# library(egg)
# #Plot PCA including expType and not including expType
# PCANoExpType <- CalcPCA(SalmonTPM_Combat$SalmonTPM_Gene_Combat_Merged)
# PCAPlotNoExpType_Instrument <- autoplot(PCANoExpType$PCA, data = PCANoExpType$Df, colour = 'Instrument',size=4,shape=F) + 
#   ggtitle('PCA of Non-Exp Type Corrected - Instrument')
# PCAPlotNoExpType_expType <- autoplot(PCANoExpType$PCA, data = PCANoExpType$Df, colour = 'expType',size=4,shape=F) + 
#                                    ggtitle('PCA of Non-Exp Type Corrected - Exp Type')
# PCAPlotNoExpType <- ggarrange(PCAPlotNoExpType_Instrument,PCAPlotNoExpType_expType,ncol=2)
# 
# PCAExpType <- CalcPCA(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged)
# PCAPlotExpType_Instrument <- autoplot(PCAExpType$PCA, data = PCAExpType$Df, colour = 'Instrument',size=4,shape=F) + ggtitle('PCA of Exp Type Corrected - Instrument')
# PCAPlotExpType_ExpType <- autoplot(PCAExpType$PCA, data = PCAExpType$Df, colour = 'expType',size=4,shape=F) + ggtitle('PCA of Exp Type Corrected - ExpType')
# PCAPlotExpType <- ggarrange(PCAPlotExpType_Instrument,PCAPlotExpType_ExpType,ncol=2)
# 
# #Boxplot of batch corrected  TPM
# TPMExpCorrected<- ggplot(stack(as.data.frame(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged)) %>% 
#                        dplyr::left_join(CollectMetadata(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged,full=F),by=c('ind' = 'Sample_Name'))) + 
#   geom_boxplot(aes(x = ind, y = values,fill=expType),outlier.shape = NA) + scale_x_discrete(name = "Sample") +
#   scale_y_continuous(name = "log(TPM+1)",limits = c(-2,4)) + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.6,size =15)) + scale_colour_gradientn(colours = rainbow(7)) + ggtitle('With Exp Type Correction')
# TPMNonExpCorrected<- ggplot(stack(as.data.frame(SalmonTPM_Combat$SalmonTPM_Gene_Combat_Merged)) %>% 
#                            dplyr::left_join(CollectMetadata(SalmonTPM_Combat$SalmonTPM_Gene_Combat_Merged,full=F),by=c('ind' = 'Sample_Name'))) + 
#   geom_boxplot(aes(x = ind, y = values,fill=expType),outlier.shape = NA) + scale_x_discrete(name = "Sample") +
#   scale_y_continuous(name = "log(TPM+1)",limits = c(-2,4)) + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.6,size =15)) + scale_colour_gradientn(colours = rainbow(7)) + ggtitle('No Exp Type Correction')
# TPMExpCorVsNonExpCor <-  ggarrange(TPMExpCorrected,TPMNonExpCorrected,nrow=2)
# 
# #Scatterplot of median exp of genes Galatro vs. Gosselin before batch correction 
# GalatroSamples <- log(SalmonTPM_Gene_Merged[,sapply(colnames(SalmonTPM_Gene_Merged),function(x) substr(x,1,3) == 'GSM')]+1)
# GosselinSamples <- log(SalmonTPM_Gene_Merged[,sapply(colnames(SalmonTPM_Gene_Merged),function(x) substr(x,1,3) == 'SRR')]+1)
# GalatroVsGosselin <- data.frame(Gosselin = apply(GosselinSamples,1,function(x) median(x)), Galatro = apply(GalatroSamples,1,function(x) median(x)))
# 
# library(devtools)
# source_gist("524eade46135f6348140",filename = "ggplot_smooth_func.R")
# GalatroVsGosselinPlot <- ggplot(GalatroVsGosselin, aes(x=Gosselin,y=Galatro)) +
#   geom_point() +   
#   geom_smooth(method=lm,   # Add linear regression line
#               se=FALSE) + stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) + 
#   theme_bw() + labs(x = 'Gosselin log(TPM+1)',y='Galatro log(TPM+1)') + ggtitle('Median Exp - Galatro vs Gosselin') 
# 
# # GalatroVsRandom <- data.frame(Galatro = GalatroVsGosselin$Galatro,Random = GalatroVsGosselin$Gosselin[sample(1:nrow(GalatroVsGosselin),size=nrow(GalatroVsGosselin),replace = F)])
# # GalatroVsRandomPlot <-ggplot(GalatroVsRandom, aes(x=Random,y=Galatro)) +
# #   geom_point() +   
# #   geom_smooth(method=lm,   # Add linear regression line
# #               se=FALSE,formula = formula) + stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) + 
# #   theme_bw() + labs(x = 'Gosselin Permuted log(TPM+1)',y='Galatro log(TPM+1)') + ggtitle('Median Exp - Galatro vs Gosselin') 
# # ggarrange(GalatroVsGosselinPlot,GalatroVsRandomPlot,ncol=2)
# 
# #Scatterplot of median exp of genes Galatro vs. Gosselin after batch correction. 
# GalatroSamplesCorrected <- SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged[,sapply(colnames(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged),function(x) substr(x,1,3) == 'GSM')]
# GosselinSamplesCorrected <- SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged[,sapply(colnames(SalmonTPM_Combat_ExpCorrected$SalmonTPM_Gene_Combat_Merged),function(x) substr(x,1,3) == 'SRR')]
# GalatroVsGosselinCorrected <- data.frame(Gosselin = apply(GosselinSamplesCorrected,1,function(x) median(x)), Galatro = apply(GalatroSamplesCorrected,1,function(x) median(x)))
# GalatroVsGosselinPlotCorrected <- ggplot(GalatroVsGosselinCorrected, aes(x=Gosselin,y=Galatro)) +
#   geom_point() +   
#   geom_smooth(method=lm,   # Add linear regression line
#               se=FALSE) + stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) + 
#   theme_bw() + labs(x = 'Gosselin log(TPM+1)',y='Galatro log(TPM+1)') + ggtitle('Batch Corrected Median Exp - Galatro vs Gosselin') 
# 
# # GalatroVsRandomCorrected <- data.frame(Galatro = GalatroVsGosselinCorrected$Galatro,Random = GalatroVsGosselinCorrected$Gosselin[sample(1:nrow(GalatroVsGosselinCorrected),size=nrow(GalatroVsGosselinCorrected),replace = F)])
# # GalatroVsRandomPlotCorrected <- ggplot(GalatroVsRandomCorrected, aes(x=Random,y=Galatro)) +
# #   geom_point() +   
# #   geom_smooth(method=lm,   # Add linear regression line
# #               se=FALSE,formula = formula) + stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) + 
# #   theme_bw() + labs(x = 'Gosselin Permuted log(TPM+1)',y='Galatro log(TPM+1)') + ggtitle('Batch Corrected Median Exp - Galatro vs Gosselin') 
# ggarrange(GalatroVsGosselinPlot,GalatroVsGosselinPlotCorrected,ncol=2)
# 
# ####################################################COMPARE WITH WHOLE BRAIN###############################################################
# 
# #Boxplot of Brain Samples
# load('../Count_Data/Galatro/SalmonTPM_Gene_WholeBrain.rda')
# BulkBrainMetadata <- CollectMetadata(SalmonTPM_Gene_WholeBrain)
# BulkBrainTPM <- ggplot(stack(as.data.frame(SalmonTPM_Gene_WholeBrain)) %>% 
#                        dplyr::left_join(CollectMetadata(SalmonTPM_Gene_WholeBrain),by=c('ind' = 'Sample_Name')) %>% dplyr::mutate(percIntergenic = intergenicTags/totalTags)) + 
#   geom_boxplot(aes(x = ind, y = values,fill=exonTags > 4e6 & percIntergenic < 0.25),outlier.shape = NA) + scale_x_discrete(name = "Sample") +
#   scale_y_continuous(name = "TPM",limits = c(0,4)) + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.6,size =15)) + scale_colour_gradientn(colours = rainbow(7))
# BulkBrainPCA <- CalcPCA(SalmonTPM_Gene_WholeBrain)
# BulkBrainPCA$Df$percIntergenic <- BulkBrainPCA$Df$intergenicTags / BulkBrainPCA$Df$totalTags
# BulkBrainPCAPlot <- autoplot(BulkBrainPCA$PCA, data = BulkBrainPCA$Df,colour='percIntergenic',size=4,shape=F) + ggtitle('PCA of Bulk Brain') 
# 
# #Remove bad Samples, and compare with Galatro and Gosselin
# SalmonTPM_Gene_WholeBrain <- SalmonTPM_Gene_WholeBrain[,BulkBrainMetadata %>%
#                                                          dplyr::mutate(percIntergenic = intergenicTags/totalTags) %>%
#                                                          dplyr::filter(percIntergenic < 0.25) %>% dplyr::select(Sample_Name) %>% t() %>% as.vector()]
# SalmonTPM_Gene_WholeBrain <- SalmonTPM_Gene_WholeBrain[rownames(GalatroSamples),]
# WholeBrainVsGalatro <- data.frame(Galatro = apply(GalatroSamples,1,function(x) median(x)),WholeBrain = apply(log(SalmonTPM_Gene_WholeBrain+1),1,function(x) median(x)))
# WholeBrainVsGosselin <- data.frame(Gosselin = apply(GosselinSamples,1,function(x) median(x)),WholeBrain = apply(log(SalmonTPM_Gene_WholeBrain+1),1,function(x) median(x)))
# WholeBrainVsGalatroPlot <- ggplot(WholeBrainVsGalatro, aes(x=WholeBrain,y=Galatro)) +
#   geom_point() +   
#   geom_smooth(method=lm,   # Add linear regression line
#               se=FALSE) + stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) + 
#   theme_bw() + labs(x = 'Whole Brain log(TPM+1)',y='Galatro log(TPM+1)') + ggtitle('Median Exp - Galatro vs WholeBrain') 
# 
# WholeBrainVsGosselinPlot <- ggplot(WholeBrainVsGosselin, aes(x=WholeBrain,y=Gosselin)) +
#   geom_point() +   
#   geom_smooth(method=lm,   # Add linear regression line
#               se=FALSE) + stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) + 
#   theme_bw() + labs(x = 'Whole Brain log(TPM+1)',y='Gosselin log(TPM+1)') + ggtitle('Median Exp - Gosselin vs WholeBrain') 
# 
# ggarrange(WholeBrainVsGalatroPlot,WholeBrainVsGosselinPlot,ncol=2)
# 
# #PCA plot comparing within Galatro
# PCAWholeBrainVsGalatro <- CalcPCA(cbind(log(SalmonTPM_Gene_WholeBrain+1),GalatroSamples))
# PCAWholeBrainVsGalatro$Df$Dataset <- ifelse(PCAWholeBrainVsGalatro$Df$readLength == 202 & PCAWholeBrainVsGalatro$Df$BulkBrain == 0,'GalatroIllumina',ifelse(PCAWholeBrainVsGalatro$Df$readLength == 252 & PCAWholeBrainVsGalatro$Df$BulkBrain == 0,'GalatroTakara','Whole Brain'))
# autoplot(PCAWholeBrainVsGalatro$PCA, data = PCAWholeBrainVsGalatro$Df, colour = 'Dataset',size=4,shape=F) + 
#   ggtitle('PCA of Tissue Type within Galatro')
# 
# #PCA plot adjusting according to Batch1 and Batch2
# WholeBrainandGalatro <- cbind(log(SalmonTPM_Gene_WholeBrain+1),GalatroSamples)
# WholeBrainandGalatro <- WholeBrainandGalatro[!apply(WholeBrainandGalatro,1,function(x) all(x==0)),]
# BatchCorrectedWholeBrainGalatro <- RunCombat(WholeBrainandGalatro,Samples='Galatro')
# PCABatchCorrected <- CalcPCA(BatchCorrectedWholeBrainGalatro$SalmonTPM_Gene_Combat_Merged)
# PCABatchCorrected$Df$Dataset <- ifelse(PCAWholeBrainVsGalatro$Df$readLength == 202 & PCAWholeBrainVsGalatro$Df$BulkBrain == 0,'GalatroIllumina',ifelse(PCAWholeBrainVsGalatro$Df$readLength == 252 & PCAWholeBrainVsGalatro$Df$BulkBrain == 0,'GalatroTakara','Whole Brain'))
# autoplot(PCABatchCorrected$PCA, data = PCABatchCorrected$Df, colour = 'Dataset',size=4,shape=F) + 
#   ggtitle('PCA of Tissue Type within Galatro - Batch Corrected')

