#!/bin/bash
curSample=`echo $1 | awk -F "_" '{print$1}'`
directory="../salmon_aligned_whole_brain/$curSample"
mkdir $directory
nice -n 19 salmon quant -i ../../Salmon_Index/GRCh37_salmon_index -l A --gcBias -1 $1 -2 $2 -p 12 -o $directory
