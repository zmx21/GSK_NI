# load('../Count_Data/Batch_Corrected/SalmonTPM_Gene_Combat.rda')
# load('../Count_Data/Batch_Corrected/SalmonTPM_Gene_PCAdj.rda')
# load('../Count_Data/Batch_Corrected/SalmonTPM_Transcript_Combat.rda')
# load('../Count_Data/Batch_Corrected/SalmonTPM_Transcript_PCAdj.rda')
# 
# SalmonCor_Gene_Combat <- cor(t(SalmonTPM_Gene_Combat))
# SalmonCor_Transcript_Combat <- cor(t(SalmonTPM_Transcript_Combat))
# SalmonCor_Gene_PCAdj1 <- cor(t(SalmonTPM_Gene_PCAdj$OnePC))
# SalmonCor_Gene_PCAdj2 <- cor(t(SalmonTPM_Gene_PCAdj$TwoPC))
# SalmonCor_Transcript_PCAdj1 <- cor(t(SalmonTPM_Transcript_PCAdj$OnePC))
# SalmonCor_Transcript_PCAdj2 <- cor(t(SalmonTPM_Transcript_PCAdj$TwoPC))
# 
# save(SalmonCor_Gene_Combat,SalmonCor_Transcript_Combat,file = '../Count_Data/Correlation_Matrices/Combat_Cor.rda')
# save(SalmonCor_Gene_PCAdj1,SalmonCor_Gene_PCAdj2,SalmonCor_Transcript_PCAdj1,SalmonCor_Transcript_PCAdj2,file = '../Count_Data/Correlation_Matrices/PCAdj_Cor.rda')

load('../Count_Data/Correlation_Matrices/geneFriendsCorMatrix.rda')
load('../Count_Data/Correlation_Matrices/PCAdj_Cor.rda')

upper_GeneFriends <- geneFriendsCorMatrix[upper.tri(geneFriendsCorMatrix)]
upper_Combat <- SalmonCor_Gene_Combat[upper.tri(SalmonCor_Gene_Combat)]
upper_PCAdj1 <- SalmonCor_Gene_PCAdj1[upper.tri(SalmonCor_Gene_PCAdj1)]
upper_PCAdj2 <- SalmonCor_Gene_PCAdj2[upper.tri(SalmonCor_Gene_PCAdj2)]


upper_Combat <- upper_Combat[!is.na(upper_GeneFriends)]
upper_GeneFriends <- upper_GeneFriends[!is.na(upper_GeneFriends)]
upper_Random <- upper_Combat[sample(1:length(upper_Combat),replace = F,size=length(upper_Combat))]
upper_PCAdj1 <- upper_PCAdj1[!is.na(upper_GeneFriends)]
upper_PCAdj2 <- upper_PCAdj2[!is.na(upper_GeneFriends)]

randSample <- sample(1:length(upper_Combat),size=1000000,replace=F)
compDf <- data.frame(GeneFriends = upper_GeneFriends[randSample],
                     Combat = upper_Combat[randSample],
                     Random = upper_Random[randSample],
                     PCAdj1 = upper_PCAdj1[randSample],
                     PCAdj2 = upper_PCAdj2[randSample])
p1 <- ggplot(compDf, aes(x=GeneFriends, y=Combat) ) +
  geom_bin2d(aes(fill=log(..count..)),bins = 500)+
  theme_bw() + labs(x = 'Gene Coexpression Database Cor Value',y='Combat Cor Value') + 
  ggtitle('Combat Cor Matrix \n vs. Existing Gene Co-exp Cor Matrix') + xlim(-1,1) + ylim(-1,1)

p2 <- ggplot(compDf, aes(x=GeneFriends, y=Random) ) +
  geom_bin2d(aes(fill=log(..count..)),bins = 500)+
  theme_bw() + labs(x = 'Gene Coexpression Database Cor Value',y='Random Cor Value') + 
  ggtitle('Randomly Permuted Cor Matrix \n vs. Existing Gene Co-exp Cor Matrix') + xlim(-1,1) + ylim(-1,1)
ggarrange(p1,p2,ncol=2)

p3 <- ggplot(compDf, aes(x=GeneFriends, y=PCAdj1) ) +
  geom_bin2d(aes(fill=log(..count..)),bins = 500)+
  theme_bw() + labs(x = 'Gene Coexpression Database Cor Value',y='PCAdj 1 Cor Value') + 
  ggtitle('PCAdj 1 Cor Matrix \n vs. Existing Gene Co-exp Cor Matrix') + xlim(-1,1) + ylim(-1,1)
p4 <- ggplot(compDf, aes(x=GeneFriends, y=PCAdj2) ) +
  geom_bin2d(aes(fill=log(..count..)),bins = 500)+
  theme_bw() + labs(x = 'Gene Coexpression Database Cor Value',y='PCAdj 1 and 2 Cor Value') + 
  ggtitle('PCAdj 1 and 2 Cor Matrix \n vs. Existing Gene Co-exp Cor Matrix') + xlim(-1,1) + ylim(-1,1)
ggarrange(p3,p4,ncol=2)

