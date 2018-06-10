#!/bin/bash
curSample=`echo $1 | awk -F "_" '{print$1}'`
directory="../STAR_aligned_whole_brain/$curSample"
mkdir $directory
#Exact Param on Galatro et al
STAR --runThreadN 2 --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignEndsType EndToEnd --sjdbScore 1 --genomeDir ../../STAR_Index --outFileNamePrefix "$directory/" --readFilesIn $1 $2
