#Distribution of correlation coefficients.
load('../../Count_Data/Correlation_Matrices/MicrogliaGeneCodingCorr.rda')
CodingUpperDiag <- MicrogliaGeneCodingCorr$r[upper.tri(MicrogliaGeneCodingCorr$r,diag = F)]
load('../../Count_Data/Correlation_Matrices/MicrogliaGeneAllCorr.rda')
AllUpperDiag <- MicrogliaGeneAllCorr$r[upper.tri(MicrogliaGeneAllCorr$r,diag = F)]

par(mfrow=c(1,3),cex.axis=1.5,cex.axis=1.5,cex.main=3,cex.lab=1.5)
hist(CodingUpperDiag,xlab = 'Pearson Correlation Coefficients',main='',cex=5)
title('a)',adj=0)
hist(AllUpperDiag,xlab = 'Pearson Correlation Coefficients',main='',cex=5)
title('b)',adj=0)
hist(AllUpperDiag,xlab = 'Pearson Correlation Coefficients',main='',cex=5)
title('c)',adj=0)
