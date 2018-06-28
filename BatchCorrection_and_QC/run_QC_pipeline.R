source('filter_by_mapping_metrics.R')
source('filter_by_TPM.R')
source('batch_adjustment.R')
source('plot_pairwise_TPM.R')
source('filter_by_CV.R')
#Remove samples which does not pass mapping criteria. Remove rRNA
FilterAllSamplesByMappingRate(plots=T)
#Remove genes which are below 75 percentile is all samples.
FilterByTPM(plots = T,inputPath = '../../Count_Data/Read_Filtered/',
            outputPath = '../../Count_Data/TPM_Filtered/',percentile=0.75)
#Batch correct individual dataset, and global dataset correction
RunBatchForAllData(plots = F)
#Save pairwise TPM comparison plots
# SavePairwiseTPMPlots()
#Keep genes which CV above median, according to distribution of coding and non-coding.
FilterByCV(plots=T,list('../../Count_Data/Batch_Corrected/','../../Count_Data/TPM_Filtered/'),
           outputPath = '../../Count_Data/CV_Filtered/',0.5)
