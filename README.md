# Integrating transcriptional networks and GWAS for identification of neurodegeneration-driving gene modules in microglia.
## MPhil Computational Biology Thesis Project (University of Cambridge), conducted at GSK

###### bash_scripts: scripts to run Salmon and STAR

###### AlignmentProcessing: 
                    importing alignment results into R

###### BatchCorrection_QC: 
                    Removing technical variability, and filtering genes based on TPM and CV. 
                    In order, should be performed as 
                    1. Mapping statistics filtering (incl rRNA gene removal)
                    2. TPM gene filtering (Where a gene is filtered if it is lowly expressed in all samples)
                    3. Batch Correction
                    4. CV Filtering (Seperating coding and non-coding genes).

###### Louvain: 
                    scripts for louvain clustering, from correlation network. 

###### PASCAL: 
                    scripts for co-expression cluster enrichment analysis, which takes in output from Louvain. 

###### HotNet2: 
                    scripts to run HotNet2 from co-expression network. 
