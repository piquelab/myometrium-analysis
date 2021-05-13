

library(Seurat)
library(Matrix)
library(tidyverse)

library(tidyverse)
library(dplyr)
library(stringr)



outFolder <- paste0("./14_deconvolution_analysis/")
system(paste0("mkdir -p ",outFolder))

######################################################
# bulk data
######################################################
load("TLTNLmyoToCaseWest.rdata")

eset_bulk<-eset
colnames(eset_bulk)<-unlist(strsplit(colnames(eset),"_"))[seq(2,2*length(colnames(eset)),by=2)]
rownames(eset_bulk)<-unlist(strsplit(rownames(eset),"_"))[seq(2,2*length(rownames(eset)),by=2)]

write_rds(eset_bulk,file=paste0(outFolder,"eset_bulk.rds"))


######################################################
# single cell data
######################################################
