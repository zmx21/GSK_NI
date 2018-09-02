##########################################################################################################
# Fig. 2.1.3.1
##########################################################################################################

#Import read dist data from salmon output
readDist_Galatro <- read.table(file='../../Galatro/STAR_aligned_merged/multiqc_data/multiqc_rseqc_read_distribution.txt',header=T,stringsAsFactors = F) 
readDist_Gosselin <- read.table(file='../../Gosselin/multiqc_rseqc_read_distribution.txt',header=T,stringsAsFactors = F)
readDist_Olah <- read.table(file='../../Olah/multiqc_rseqc_read_distribution.txt',header=T,stringsAsFactors = F)
readDist_Galatro$Sample <- sapply(readDist_Galatro$Sample,function(x) strsplit(x,split = '[.]')[[1]][1])

#Generate plots
library(reshape2)
library(RColorBrewer)

Gosselin_melted <- melt(readDist_Gosselin %>% dplyr::select(Sample,contains("_tags_kb")), id='Sample')
Gosselin_melted$variable <- sapply(Gosselin_melted$variable,function(x) gsub('_tags_kb','',x))
Gosselin_melted$value <- sapply(Gosselin_melted$value,function(x) x/1e6)

GosselinPlot <- ggplot(Gosselin_melted) + aes(x=Sample,y = value) + 
  geom_bar(aes(fill=variable),stat = 'identity') + coord_flip() +  
  scale_fill_brewer(palette = 'Spectral') + 
  labs(y='Number of Tags (Millions)',x='Sample Name',fill='Tag Type')
ggsave(GosselinPlot,file = '../../FinalFigures/GosselinReadDist.pdf',device = 'pdf',units = 'in',height = 5,width = 7)

Galatro_melted <- melt(readDist_Galatro %>% dplyr::select(Sample,contains("_tags_kb")), id='Sample')
Galatro_melted$variable <- sapply(Galatro_melted$variable,function(x) gsub('_tags_kb','',x))
Galatro_melted$value <- sapply(Galatro_melted$value,function(x) x/1e6)

GalatroPlot <- ggplot(Galatro_melted) + aes(x=Sample,y = value) + 
  geom_bar(aes(fill=variable),stat = 'identity') + coord_flip() +  
  scale_fill_brewer(palette = 'Spectral') + 
  labs(y='Number of Tags (Millions)',x='Sample Name',fill='Tag Type')
ggsave(GalatroPlot,file = '../../FinalFigures/GalatroReadDist.pdf',device = 'pdf',units = 'in',height = 5,width = 7)


Olah_melted <- melt(readDist_Olah %>% dplyr::select(Sample,contains("_tags_kb")), id='Sample')
Olah_melted$variable <- sapply(Olah_melted$variable,function(x) gsub('_tags_kb','',x))
Olah_melted$value <- sapply(Olah_melted$value,function(x) x/1e6)

OlahPlot <- ggplot(Olah_melted) + aes(x=Sample,y = value) + 
  geom_bar(aes(fill=variable),stat = 'identity') + coord_flip() +  
  scale_fill_brewer(palette = 'Spectral') + 
  labs(y='Number of Tags (Millions)',x='Sample Name',fill='Tag Type')
ggsave(OlahPlot,file = '../../FinalFigures/OlahReadDist.pdf',device = 'pdf',units = 'in',height = 3,width = 7)