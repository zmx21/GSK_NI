#!/bin/bash
curSample=`echo $1 | awk -F "_" '{print$1}'`
mkdir ../Salmon_aligned/Salmon_aligned_fmd/$curSample
salmon quant -i ../../Salmon_Index/GRCh37_salmon_index_fmd/ --gcBias -l A -1 $1 -2 $2 -p 12 -o ../Salmon_aligned/Salmon_aligned_fmd/$curSample
