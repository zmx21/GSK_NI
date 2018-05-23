#!/bin/bash
curSample=`echo $1 | awk -F "_" '{print$1}'`
mkdir ../STAR_aligned_merged/$curSample
STAR --runThreadN 10 --genomeDir ../../STAR_Index --outFileNamePrefix ../STAR_aligned_merged/$curSample/ --readFilesIn $1 $2 
