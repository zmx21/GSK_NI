#!/bin/bash
salmon index -t ../GRCh37_Ensembl75/Homo_sapiens.GRCh37.75.cdna.ncrna.fa -i ../Salmon_Index/GRCh37_salmon_index_k19 -k 19
salmon index -t ../GRCh37_Ensembl75/Homo_sapiens.GRCh37.75.cdna.ncrna.fa -i ../Salmon_Index/GRCh37_salmon_index_fmd --type fmd
