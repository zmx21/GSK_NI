##########################################################################################################
#Filters count matrix based on TPM and CV, and generates various plots. 
##########################################################################################################

#Load required data
library(dplyr)
# source('import_salmon.R')
# tx2gene <- ImportTx2gene()
# SalmonTPM <- ImportSalmonCounts('/local/data/public/zmx21/zmx21_private/GSK/Galatro/Salmon_aligned/Salmon_aligned_merged',tx2gene)
# SalmonTPM_Gene <- SalmonTPM$geneLevel$abundance
# SalmonTPM_Transcript <- SalmonTPM$transcriptLevel$abundance
load(file = '../Count_Data/SalmonTPM_Transcript_Microglia.rda')
load(file = '../Count_Data/SalmonTPM_Gene_Microglia.rda')
load(file = '../Count_Data/gtfTables.rda')
readDist <- read.table(file='/local/data/public/zmx21/zmx21_private/GSK/Galatro/STAR_aligned_merged/multiqc_data/multiqc_rseqc_read_distribution.txt',
                       header = T,stringsAsFactors = F)
readDist$Sample <- sapply(readDist$Sample,function(x) unlist(strsplit(x,".",fixed = T))[1])
readData <- cbind(data.frame(Sample=readDist$Sample),
                  data.frame(PercExonic = round(as.numeric(readDist$cds_exons_tag_count/readDist$total_tags),2),
                             ExonicTags = round(as.numeric(readDist$cds_exons_tag_count),2),
                             PercIntergenic = round(as.numeric((readDist$other_intergenic_tag_count)/readDist$total_tags),2)))

#Returns a filtered matrix based on cutoff of TPM and CV.
#If a mapping table to category is included, then returns the category of each remaining transcript/gene also
FilterCountMatrix <- function(countMatrixInput,meanAbundances,cv,cvCutOffAbsolute,meanCutOffAbsolute,mappingTable=NULL){
  #Filter for genes according to cutoff
  belowMeanCutOff <- which(meanAbundances <= meanCutOffAbsolute)
  belowCVCutoff <- which(cv <= cvCutOffAbsolute | is.nan(cv))
  rowsToFilter <- union(belowMeanCutOff,belowCVCutoff)
  
  filteredCountMatrix <- countMatrixInput[setdiff(1:nrow(countMatrixInput),rowsToFilter),]
  
  #Decide whether to return categorial counts.
  if(is.null(mappingTable)){
    return(list(countMatrix = filteredCountMatrix,numFiltered = length(rowsToFilter)))
  }else{
    type <- ifelse(grepl('ENST',rownames(countMatrixInput)[1]),'transcript_id','gene_id')
    allFilteredNames <- data.frame(ids = rownames(filteredCountMatrix),stringsAsFactors = F)
    colnames(allFilteredNames) <- type
    categoricalCounts <- as.data.frame(dplyr::left_join(allFilteredNames,mappingTable,by=c(type,type)))
    return(list(countMatrix = filteredCountMatrix,numFiltered = length(rowsToFilter),categoricalCounts = categoricalCounts))
  }
}

source('load_GTF.R')
gtfTables <- LoadBiotypeMapping()
#Get Gene Level Stats
meanAbundancesGene <- rowSums(SalmonTPM_Gene) / ncol(SalmonTPM_Gene)
cvGene <- apply(SalmonTPM_Gene,1,function(x) sd(x) / mean(x))
#Get Transcript Level Stats
meanAbundancesTranscript <- rowSums(SalmonTPM_Transcript) / ncol(SalmonTPM_Transcript)
cvTranscript <- apply(SalmonTPM_Transcript,1,function(x) sd(x) / mean(x))

########################################Plots ##########################################
#Box Plot of TPM and CV of each sample, at Transcript and Gene Level
par(mfrow=c(2,2))
boxplot(meanAbundancesGene,ylim=c(0,25),main='Mean Abundance - Gene Level',cex.axis=2)
boxplot(cvGene,main='Coefficient of Variation - Gene Level',cex.axis=2)
boxplot(meanAbundancesTranscript,ylim=c(0,10),main='Mean Abundance - Transcript Level',cex.axis=2)
boxplot(cvTranscript,main='Coefficient of Variation - Transcript Level',cex.axis=2)

#PCA plot of samples, after filtering based on Percent Exonic and Intergenic.
source('pca_analysis.R')
filtered1_PCA <- CalcPCA(SalmonTPM_Gene[,readData$PercExonic > 0.1],tables)
filtered1 <- autoplot(filtered1_PCA$PCA, data = filtered1_PCA$Df, colour = 'readLength',size=4,shape=F) + ggtitle('PCA - PercExonic > 0.1  & ExonicTags > 0') 
filtered2_PCA <- CalcPCA(SalmonTPM_Gene[,readData$PercExonic > 0.15],tables)
filtered2 <- autoplot(filtered2_PCA$PCA, data = filtered2_PCA$Df, colour = 'readLength',size=4,shape=F) + ggtitle('PCA - PercExonic > 0.15 & ExonicTags > 0') 
filtered3_PCA <- CalcPCA(SalmonTPM_Gene[,readData$PercExonic > 0.1 & readData$ExonicTags> 5000000],tables)
filtered3 <- autoplot(filtered3_PCA$PCA, data = filtered3_PCA$Df, colour = 'readLength',size=4,shape=F) + ggtitle('PCA - PercExonic > 0.1 & ExonicTags > 5e+06') 
filtered4_PCA <- CalcPCA(SalmonTPM_Gene[,readData$PercExonic > 0.15 & readData$ExonicTags> 5000000],tables)
filtered4 <- autoplot(filtered4_PCA$PCA, data = filtered4_PCA$Df, colour = 'readLength',size=4,shape=F) + ggtitle('PCA - PercExonic > 0.15 & ExonicTages > 5e+06') 
ggarrange(filtered1,filtered2,filtered3,filtered4,ncol=2,nrow=2)

#Categorical bar plot of biotype, without cutoff.
GeneLevelNoCutOff<- FilterCountMatrix(countMatrix = SalmonTPM_Gene,meanAbundances = meanAbundancesGene,cv = cvGene,cvCutOffAbsolute = 0,meanCutOffAbsolute = 0,mappingTable = gtfTables$geneTable)
ggplot(data=GeneLevelNoCutOff$categoricalCounts, aes(x=as.factor(gene_biotype))) + 
  geom_bar() + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15),axis.title = element_text(size=15),axis.text =  element_text(size=15)) + ggtitle('Gene Biotype - No Cutoff')+
  xlab('Gene Biotype') +  scale_y_continuous("Number of Genes",scales::pretty_breaks(n = 10),sec.axis = sec_axis(~./nrow(SalmonTPM_Gene), name = 'Fraction of All Genes',breaks = scales::pretty_breaks(n = 10)))
TranscriptLevelNoCutOff<- FilterCountMatrix(countMatrix = SalmonTPM_Transcript,meanAbundances = meanAbundancesTranscript,cv = cvTranscript,cvCutOffAbsolute = 0,meanCutOffAbsolute = 0,mappingTable = gtfTables$transcriptTable)
ggplot(data=TranscriptLevelNoCutOff$categoricalCounts, aes(x=as.factor(source))) + 
  geom_bar() + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15),axis.title = element_text(size=15),axis.text =  element_text(size=15)) + ggtitle('Transcript Source - No Cutoff')+
  xlab('Transcript Source') +  scale_y_continuous("Number of Transcripts",scales::pretty_breaks(n = 10),sec.axis = sec_axis(~./nrow(SalmonTPM_Transcript), name = 'Fraction of All Transcripts',breaks = scales::pretty_breaks(n = 10)))
NullGeneLevel <- ggplot(data=gtfTables$geneTable, aes(x=as.factor(gene_biotype))) + 
  geom_bar() + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15),axis.title = element_text(size=15),axis.text =  element_text(size=15)) + ggtitle('Gene Biotype - Null Dist')+
  xlab('Gene Biotype') +  scale_y_continuous("Number of Genes",scales::pretty_breaks(n = 10),sec.axis = sec_axis(~./nrow(SalmonTPM_Gene), name = 'Fraction of All Genes',breaks = scales::pretty_breaks(n = 10)))
NullTranscriptLevel <- ggplot(data=gtfTables$transcriptTable, aes(x=as.factor(source))) + 
  geom_bar() + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15),axis.title = element_text(size=15),axis.text =  element_text(size=15)) + ggtitle('Transcript Source - Null Dist')+
  xlab('Gene Biotype') +  scale_y_continuous("Number of Transcripts",scales::pretty_breaks(n = 10),sec.axis = sec_axis(~./nrow(SalmonTPM_Gene), name = 'Fraction of All Transcripts',breaks = scales::pretty_breaks(n = 10)))

#Mean and CV cutoff trends
meanCutoffStatsGene <- lapply(seq(0,20,1),function(x) FilterCountMatrix(countMatrix = SalmonTPM_Gene,
                                                               meanAbundances = meanAbundancesGene,
                                                               cv = cvGene,
                                                               cvCutOffAbsolute = 0,
                                                               meanCutOffAbsolute = x,mappingTable = gtfTables$geneTable) %>%
                                                               {ParseBiotypeTable(.$categoricalCounts)} %>%
                                                               {cbind(data.frame(cutoff=rep(x,nrow(.))),.)})
cvCutoffStatsGene <- lapply(seq(0,5,0.2),function(x) FilterCountMatrix(countMatrix = SalmonTPM_Gene,
                                                                          meanAbundances = meanAbundancesGene,
                                                                          cv = cvGene,
                                                                          cvCutOffAbsolute = x,
                                                                          meanCutOffAbsolute = 0,mappingTable = gtfTables$geneTable)%>%
                                                                          {ParseBiotypeTable(.$categoricalCounts)}%>%
                                                                          {cbind(data.frame(cutoff=rep(x,nrow(.))),.)})
meanCutoffStatsTranscript <- lapply(seq(0,10,0.5),function(x) FilterCountMatrix(countMatrix = SalmonTPM_Transcript,
                                                                            meanAbundances = meanAbundancesTranscript,
                                                                            cv = cvTranscript,
                                                                            cvCutOffAbsolute = 0,
                                                                            meanCutOffAbsolute = x,mappingTable = gtfTables$transcriptTable)%>%
                                                                            {ParseBiotypeTable(.$categoricalCounts)}%>%
                                                                            {cbind(data.frame(cutoff=rep(x,nrow(.))),.)})
cvCutoffStatsTranscript <- lapply(seq(0,5,0.2),function(x) FilterCountMatrix(countMatrix = SalmonTPM_Transcript,
                                                                          meanAbundances = meanAbundancesTranscript,
                                                                          cv = cvTranscript,
                                                                          cvCutOffAbsolute = x,
                                                                          meanCutOffAbsolute = 0,mappingTable = gtfTables$transcriptTable)%>%
                                                                          {ParseBiotypeTable(.$categoricalCounts)}%>%
                                                                          {cbind(data.frame(cutoff=rep(x,nrow(.))),.)})

#Num of Gene/Transcripts vs Cutoffs
ggplot(do.call(rbind,meanCutoffStatsGene), 
       aes(fill=as.factor(gene_biotype), x=as.factor(cutoff))) + 
  geom_bar(position="dodge") + theme(axis.text.x = element_text(size=15),axis.title = element_text(size=15),axis.text =  element_text(size=15)) + ggtitle('Gene Biotype vs TPM Cutoff')+
  labs(x='Cutoff (TPM)',fill='Gene Biotype') +  scale_y_continuous("Number of Genes",scales::pretty_breaks(n = 10),sec.axis = sec_axis(~./nrow(SalmonTPM_Gene), name = 'Fraction of All Genes',breaks = scales::pretty_breaks(n = 10)))

ggplot(do.call(rbind,cvCutoffStatsGene), 
       aes(fill=as.factor(gene_biotype), x=as.factor(cutoff))) + 
  geom_bar(position="dodge") + theme(axis.text.x = element_text(size=15),axis.title = element_text(size=15),axis.text =  element_text(size=15)) + ggtitle('Gene Biotype vs CV Cutoff')+
  labs(x='Cutoff (CV)',fill='Gene Biotype') +  scale_y_continuous("Number of Genes",scales::pretty_breaks(n = 10),sec.axis = sec_axis(~./nrow(SalmonTPM_Gene), name = 'Fraction of All Genes',breaks = scales::pretty_breaks(n = 10)))

ggplot(do.call(rbind,meanCutoffStatsTranscript), 
       aes(fill=as.factor(source), x=as.factor(cutoff))) + 
  geom_bar(position="dodge") + theme(axis.text.x = element_text(size=15),axis.title = element_text(size=15),axis.text =  element_text(size=15)) + ggtitle('Source vs TPM Cutoff')+
  labs(x='Cutoff (TPM)',fill='Source') +  scale_y_continuous("Number of Transcripts",scales::pretty_breaks(n = 10),sec.axis = sec_axis(~./nrow(SalmonTPM_Transcript), name = 'Fraction of All Transcripts',breaks = scales::pretty_breaks(n = 10)))

ggplot(do.call(rbind,cvCutoffStatsTranscript), 
       aes(fill=as.factor(source), x=as.factor(cutoff))) + 
  geom_bar(position="dodge") + theme(axis.text.x = element_text(size=15),axis.title = element_text(size=15),axis.text =  element_text(size=15)) + ggtitle('Source vs CV Cutoff')+
  labs(x='Cutoff (CV)',fill='Source') +  scale_y_continuous("Number of Transcripts",scales::pretty_breaks(n = 10),sec.axis = sec_axis(~./nrow(SalmonTPM_Transcript), name = 'Fraction of All Transcripts',breaks = scales::pretty_breaks(n = 10)))

#HeatMaps. X and Y are TPM and CV cutoff. Z is the number of gene/transcripts remaining
geneLevelSearchSpace <- expand.grid(TPMCutoff = seq(0,20,0.5),CVCutOff=seq(0,1,0.025))
heatMapGeneLevel <- mclapply(1:nrow(geneLevelSearchSpace),function(i) FilterCountMatrix(countMatrix = SalmonTPM_Gene,
                                                                        meanAbundances = meanAbundancesGene,
                                                                        cv = cvGene,
                                                                        cvCutOffAbsolute = geneLevelSearchSpace$CVCutOff[i],
                                                                        meanCutOffAbsolute = geneLevelSearchSpace$TPMCutoff[i],
                                                                        mappingTable = gtfTables$geneTable) %>%
                                                                        {ParseBiotypeTable(.$categoricalCounts)} %>%
                                                                        {table(.$gene_biotype)},mc.cores = 50) %>% {as.data.frame(do.call(rbind,.))}
heatMapGeneLevel <- cbind(geneLevelSearchSpace,heatMapGeneLevel)
library(gridExtra)
p1Gene <- ggplot(heatMapGeneLevel,aes(x = as.numeric(TPMCutoff), y = as.numeric(CVCutOff),fill = protein_coding)) + 
  geom_tile() + scale_fill_gradient() 
p2Gene <- ggplot(heatMapGeneLevel,aes(x = as.numeric(TPMCutoff), y = as.numeric(CVCutOff),fill = lncRNA)) + 
  geom_tile() + scale_fill_gradient() 
p3Gene <- ggplot(heatMapGeneLevel,aes(x = as.numeric(TPMCutoff), y = as.numeric(CVCutOff),fill = pseudogene)) + 
  geom_tile() + scale_fill_gradient() 
p4Gene <- ggplot(heatMapGeneLevel,aes(x = as.numeric(TPMCutoff), y = as.numeric(CVCutOff),fill = sncRNA)) + 
  geom_tile() + scale_fill_gradient() 
grid.arrange(p1Gene, p2Gene,p3Gene,p4Gene,top = textGrob("Cutoff Heat Map - Gene Level",gp=gpar(fontsize=20)))
transcriptLevelSearchSpace <- expand.grid(TPMCutoff = seq(0,10,0.5),CVCutOff=seq(0,2,0.025))
heatMapTranscriptLevel <- mclapply(1:nrow(transcriptLevelSearchSpace),function(i) FilterCountMatrix(countMatrix = SalmonTPM_Transcript,
                                                                                        meanAbundances = meanAbundancesTranscript,
                                                                                        cv = cvTranscript,
                                                                                        cvCutOffAbsolute = transcriptLevelSearchSpace$CVCutOff[i],
                                                                                        meanCutOffAbsolute = transcriptLevelSearchSpace$TPMCutoff[i],
                                                                                        mappingTable = gtfTables$transcriptTable) %>%
                                                                                        {ParseBiotypeTable(.$categoricalCounts)} %>%
                                                                                        {table(.$source)},mc.cores = 50) %>% {as.data.frame(do.call(rbind,.))}
heatMapTranscriptLevel <- cbind(transcriptLevelSearchSpace,heatMapTranscriptLevel)
library(gridExtra)
p1Transcript <- ggplot(heatMapTranscriptLevel,aes(x = as.numeric(TPMCutoff), y = as.numeric(CVCutOff),fill = protein_coding)) + 
  geom_tile() + scale_fill_gradient() 
p2Transcript <- ggplot(heatMapTranscriptLevel,aes(x = as.numeric(TPMCutoff), y = as.numeric(CVCutOff),fill = lncRNA)) + 
  geom_tile() + scale_fill_gradient() 
p3Transcript <- ggplot(heatMapTranscriptLevel,aes(x = as.numeric(TPMCutoff), y = as.numeric(CVCutOff),fill = pseudogene)) + 
  geom_tile() + scale_fill_gradient() 
p4Transcript <- ggplot(heatMapTranscriptLevel,aes(x = as.numeric(TPMCutoff), y = as.numeric(CVCutOff),fill = sncRNA)) + 
  geom_tile() + scale_fill_gradient() 
grid.arrange(p1Transcript, p2Transcript,p3Transcript,p4Transcript,top = textGrob("Cutoff Heat Map - Transcript Level",gp=gpar(fontsize=20)))

#Different Filtering Levels, get categorical counts. 
geneLevelCutOffTPM <- c(5,5,quantile(meanAbundancesGene)[4],quantile(meanAbundancesGene)[4])
geneLevelCutOffCV <- c(0.3,0.5,0.3,0.5)
transcriptLevelCutOffTPM <- c(5,5,quantile(meanAbundancesGene)[4],quantile(meanAbundancesGene)[4])
transcriptLevelCutOffCV <- c(0.3,0.5,0.3,0.5)

geneLevelFiltCounts <- lapply(1:length(geneLevelCutOffTPM), function(i) FilterCountMatrix(countMatrix = SalmonTPM_Gene,
                                     meanAbundances = meanAbundancesGene,
                                     cv = cvGene,
                                     cvCutOffAbsolute = geneLevelCutOffCV[i],
                                     meanCutOffAbsolute = geneLevelCutOffTPM[i],
                                     mappingTable = gtfTables$geneTable) %>%
                                     {ParseBiotypeTable(.$categoricalCounts)} %>%
                                     {table(.$gene_biotype)})

transLevelFiltCounts <- lapply(1:length(geneLevelCutOffTPM), function(i) FilterCountMatrix(countMatrix = SalmonTPM_Transcript,
                                     meanAbundances = meanAbundancesTranscript,
                                     cv = cvTranscript,
                                     cvCutOffAbsolute = transcriptLevelCutOffCV[i],
                                     meanCutOffAbsolute = transcriptLevelCutOffTPM[i],
                                     mappingTable = gtfTables$transcriptTable) %>%
                                     {ParseBiotypeTable(.$categoricalCounts)} %>%
                                     {table(.$source)})


#Do final filtering
badSamples <- as.character(setdiff(readData$Sample,subset(readData,ExonicTags >= 4e+06 & PercIntergenic <= 0.25)$Sample))
SalmonTPM_Gene <- SalmonTPM_Gene[,!colnames(SalmonTPM_Gene)%in%badSamples]
SalmonTPM_Transcript <- SalmonTPM_Transcript[,!colnames(SalmonTPM_Transcript)%in%badSamples]
meanAbundancesGene <- rowSums(SalmonTPM_Gene) / ncol(SalmonTPM_Gene)
cvGene <- apply(SalmonTPM_Gene,1,function(x) sd(x) / mean(x))
meanAbundancesTranscript <- rowSums(SalmonTPM_Transcript) / ncol(SalmonTPM_Transcript)
cvTranscript <- apply(SalmonTPM_Transcript,1,function(x) sd(x) / mean(x))


SalmonTPM_Gene_Filt <-FilterCountMatrix(countMatrixInput = SalmonTPM_Gene,
                                        meanAbundances = meanAbundancesGene,
                                        cv = cvGene,cvCutOffAbsolute = 0.3,
                                        meanCutOffAbsolute = quantile(meanAbundancesGene)[4])$countMatrix
SalmonTPM_Transcript_Filt <- FilterCountMatrix(countMatrixInput = SalmonTPM_Transcript,
                                               meanAbundances = meanAbundancesTranscript,
                                               cv = cvTranscript,cvCutOffAbsolute = 0.3,
                                               meanCutOffAbsolute =quantile(meanAbundancesGene)[4])$countMatrix
SalmonTPM_Gene_Df <- dplyr::left_join(stack(as.data.frame(SalmonTPM_Gene_Filt)),
                                      readData,by=c('ind'='Sample')) %>% mutate(ind=substr(ind,9,10))
#Boxplot of samples after final filtering. 
p1 <- ggplot(SalmonTPM_Gene_Df) + 
  geom_boxplot(aes(x = ind, y = values,fill=ExonicTags >= 4e+06 & PercIntergenic <= 0.25),outlier.shape = NA) + scale_x_discrete(name = "Sample") +
  scale_y_continuous(name = "TPM",limits = quantile(SalmonTPM_Gene_Df$values, c(0.1, 0.9))) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.6,size =15)) 
p2 <- ggplot(SalmonTPM_Gene_Df) +
  geom_boxplot(aes(x = ind, y = values,fill=ExonicTags >= 5e+06 & PercIntergenic <= 0.25),outlier.shape = NA) + scale_x_discrete(name = "Sample") +
  scale_y_continuous(name = "TPM",limits = quantile(SalmonTPM_Gene_Df$values, c(0.1, 0.9))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.6,size =15))
ggarrange(p1,p2,ncol=1)



#Housekeeping Genes
HouseKeepingGenes <- c(GAPDH='ENSG00000111640',TREM2='ENSG00000095970',SYK='ENSG00000165025')
HouseKeepingExp <- t(SalmonTPM_Gene[HouseKeepingGenes,])
HouseKeepingDf <- cbind(data.frame(sample=rep(rownames(HouseKeepingExp),3)),
                        tpm = c(HouseKeepingExp[,1],HouseKeepingExp[,2],HouseKeepingExp[,3]),
                        gene=c(rep('GAPDH',ncol(SalmonTPM_Gene)),rep('TREM2',ncol(SalmonTPM_Gene)),rep('SYK',ncol(SalmonTPM_Gene))))

ggplot(HouseKeepingDf,aes(x=sample,y=as.numeric(tpm),fill=gene)) +
  geom_bar(position = 'dodge',stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15))