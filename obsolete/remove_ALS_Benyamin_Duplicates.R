library(dplyr)
ALS.Benyamin <- data.table::fread('/local/data/public/zmx21/zmx21_private/GSK/GWAS/all_studies/ALS.Benyamin.onlyrs.txt')
duplicatedRows <- ALS.Benyamin[which(duplicated(dplyr::select(ALS.Benyamin,CHR,BP))),]
ALS.Benyamin.Filt <- ALS.Benyamin
for(i in 1:nrow(duplicatedRows)){
  allRows <- dplyr::filter(ALS.Benyamin,CHR==duplicatedRows$CHR[i],BP==duplicatedRows$BP[i])
  leastSignificant <- allRows[setdiff(1:nrow(allRows),which.min(allRows$P)),]
  ALS.Benyamin.Filt <- dplyr::filter(ALS.Benyamin.Filt,!SNP%in%leastSignificant$SNP)
}
write.table(ALS.Benyamin.Filt,'/local/data/public/zmx21/zmx21_private/GSK/GWAS/all_studies/ALS.Benyamin.noduplicate.txt',row.names = F,quote = F,sep='\t')

