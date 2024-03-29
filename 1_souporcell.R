
#############################################################
# Processing the souporcell data for downstream analysis

#############################################################

library(Matrix)
library(tidyverse)
library(parallel)
library(data.table)
rm(list=ls())


outdir <- "./1_souporcell_output/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)


#################################
### 1. folders of souporcell data 
#################################


basefolder <- "/wsu/home/groups/prbgenomics/labor_myo/counts_7only/souporcell/"
basefolder <- "/wsu/home/groups/prbgenomics/labor_myo/counts_2020-12-18/souporcell/"

demuxfn <- list.files(basefolder,pattern="^Labor*",include.dirs=TRUE)
expNames<-demuxfn
 
############################
### 2, read souporcell data 
############################
demux <- mclapply(expNames,function(ii){
    ##ii <- expNames[1]  
    cat("#Loading ", ii, "\n")
    fn <- paste0(basefolder, ii,"/clusters.tsv")
    dd <- data.frame(fread(fn,header=T))
    dd <- dd %>%
      mutate(barcode=paste0(ii,"_", gsub("-1","",barcode)),EXP=ii) %>%
      mutate(assignment=paste0(ii,"_",assignment)) %>%
      select(barcode,status,assignment,log_prob_singleton,log_prob_doublet)             
    ##
    fn <- paste0("cat ",basefolder, ii,"/plink2.kin0 | grep Labor | grep HUZ")
    cat(fn,"\n")
    kin <- data.frame(fread(cmd=fn,header=F))
    kin <- kin %>% group_by(V2) %>% top_n(wt=V6,n=1)
    nn <- kin$V1
    names(nn) <- kin$V2
    nn[kin$V6<0] <- NA
    ##
    dd$assig2 <- nn[dd$assignment]
    ## head(dd)
    ## head(kin)    
    dd  
},mc.cores=10)


demux2 <- do.call(rbind,demux) #"/wsu/home/groups/prbgenomics/labor_myo/counts_7only/souporcell/"
demux1 <- do.call(rbind,demux) #"/wsu/home/groups/prbgenomics/labor_myo/counts_2020-12-18/souporcell/"

#combine 
demux<-rbind(demux2,demux1)


### output
opfn <- "./1_souporcell_output/1_souporcell.ALL.rds" 
write_rds(demux, opfn)


#output after filtering doublets
opfn <- "./1_demux2_output/1_souporcell.SNG.rds" 
demux %>% filter(status=="singlet") %>%
    write_rds(opfn)

table(demux$status)
table(demux$assig2)
