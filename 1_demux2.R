# 
# Processing the demuxlet data (parturition)

##library("rhdf5")
library(Matrix)
library(tidyverse)
library(parallel)
library(data.table)
### 

rm(list=ls())


outdir <- "./1_demux2_output/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)

##
## demux2 for downstream analysis

#################################
### 1. folders of counts data ###
#################################
#basefolder <- "/nfs/rprdata/scilab/labor2/counts.covid.2020-08-15/demuxlet/demuxlet/"
basefolder <- "../counts_2020-12-18/demuxlet/demuxlet/"
basefolder <- "../counts_7only/demuxlet/demuxlet/"

demuxfn <- list.files(basefolder,pattern="*.out.best")

#* all the libraries from the membrane placenta 
expNames <- gsub(".out.best", "", demuxfn) 
## remove bad emulsions
expNames
 


############################
### 2, read demuxlet data ###
############################

#* new column "NEW_BARCODE" was added to the data frame associated to each experiment (experiment name+barcode) 
# "-1" was removed from end of barcode

demux <- mclapply(expNames,function(ii){
  cat("#Loading ", ii, "\n")
  fn <- paste0(basefolder, ii,".out.best")
  dd <- data.frame(fread(fn,header=T))
  dd <- dd%>%mutate(NEW_BARCODE=paste0(ii,"_", gsub("-1","",BARCODE)),EXP=ii)
},mc.cores=10)


#* merging all the experiments  (12 experiments)
demux <- do.call(rbind,demux)
###

#* saving the final demux file obtained from 12 experiments into "1_demux_New.ALL.rds" file
### output
#opfn <- "./1_demux2_output/1_demux_New.ALL.rds" 
opfn <- "./1_demux2_output/1_demux_New.ALL_7only.rds" 
write_rds(demux, opfn)

# DROPLET.TYPE: "SNG" "AMB"
#opfn <- "./1_demux2_output/1_demux_New.SNG.rds" 
opfn <- "./1_demux2_output/1_demux_New.SNG_7only.rds" 

demux %>% filter(DROPLET.TYPE=="SNG") %>%
    write_rds(opfn)
