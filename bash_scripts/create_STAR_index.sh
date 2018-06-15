#!/bin/bash
./STAR/source/STAR --runThreadN 30 --runMode genomeGenerate --genomeDir ./Genome_Index \
--genomeFastaFiles Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
--sjdbGTFfile Homo_sapiens.GRCh37.75.gtf
#GTF_FILE http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
#FASTA_FILE http://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz 
