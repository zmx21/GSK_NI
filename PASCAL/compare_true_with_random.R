allStudiesMicrogliaPath <- GetPathToPASCALResults('../../GWAS/PASCAL_results/microglia_gene/')
JoinedDfMicroglia <- ParsePASCALFile(allStudiesMicrogliaPath,'../../Louvain_results/microglia_gene_clusters.gmt')
JoinedDfMicroglia <- AppendAdjustedPValue(JoinedDfMicroglia)

allStudiesRandomMicrogliaPath <- GetPathToPASCALResults('../../GWAS/PASCAL_results/randommicroglia_gene/')
JoinedDfMicrogliaRandom <- ParsePASCALFile(allStudiesRandomMicrogliaPath,'../../Louvain_results/randommicroglia_gene_clusters.gmt')
JoinedDfMicrogliaRandom <- AppendAdjustedPValue(JoinedDfMicrogliaRandom)

library(WRS)
pValDiff <- WRS::qcomhd(JoinedDfMicroglia$adjPval,JoinedDfMicrogliaRandom$adjPval,q=seq(0.001,0.005,0.001),nboot = 50)