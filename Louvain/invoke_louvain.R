#Import python function into R, can directly call. 
library(reticulate)
louvainModule <- reticulate::import_from_path('louvain','/home/zmx21/miniconda2/lib/python2.7/site-packages/')
source_python("run_louvain.py")
