#!/bin/bash
trim_galore -o ../fastq_trimmed --stringency 5 --paired $1 $2
