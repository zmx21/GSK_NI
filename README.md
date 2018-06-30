# Integrating transcriptional networks and GWAS for identification of neurodegeneration-driving gene modules in microglia.
## MPhil Computational Biology Thesis Project (University of Cambridge), conducted at GSK

#### bash_scripts: 
	Scripts to run Salmon and STAR

#### AlignmentProcessing: 
	Importing alignment results into R

#### BatchCorrection_QC: 
	Removing technical variability, and filtering genes based on TPM and CV. 
		In order, should be performed as 
		1. Mapping statistics filtering (incl rRNA gene removal)
		2. TPM gene filtering (Where a gene is filtered if it is lowly expressed in all samples)
		3. Batch Correction
		4. CV Filtering (Seperating coding and non-coding genes).

#### Louvain: 
	Scripts for louvain clustering, from correlation network. 
		- Uses python i-graph package implementation of Louvain (https://github.com/vtraag/louvain-igraph/)
		- For each round, Louvain is run recursively on every cluster remaining (above minimum cluster size)
		- The algorithm stops when convergence is reached (ie no clusters are further divided)
#### PASCAL: 
	Scripts for co-expression cluster enrichment analysis, which takes in output from Louvain. 

#### HotNet2: 
	Scripts to run HotNet2 from co-expression network. 

