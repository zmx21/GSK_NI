#!/bin/bash
curSample=`echo $1 | awk -F "_" '{print$1}'`
directory="../Salmon_aligned/Salmon_aligned_merged/$curSample"
mkdir $directory
salmon quant -i ../../Salmon_Index/GRCh37_salmon_index_k19/ -l A --gcBias -1 $1 -2 $2 -p 10 --writeUnmappedNames -o $directory --writeMappings="$directory/mapping"
