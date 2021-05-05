#
library(tidyverse)
library(parallel)
##library(data.table)

### 
outFolder="./1_demux2_output/"

#opfn <- "./1_demux2_output/1_demux_New.SNG.rds" 
opfn <- "./1_demux2_output/1_demux_New.SNG_7only.rds" 
demux <- read_rds(opfn)

# filter 
aa <- demux %>% dplyr::filter(NUM.READS>10,NUM.SNPS>10,NUM.READS<20000) %>%
    select(NEW_BARCODE,NUM.READS,NUM.SNPS,EXP,Sample_ID=SNG.BEST.GUESS)

## 


##SNG.BEST.GUESS
##
#cc <- read_tsv("ParturitionAndMyoSamples_old.txt") 
cc <- read.csv("ParturitionAndMyoSamples.csv",stringsAsFactors = FALSE) 
#%>%mutate(Sample_ID=paste(Pregnancy_ID,Origin,sep="-"))
head(cc)

aa <- aa %>% select(EXP,Sample_ID) %>% left_join(cc) %>%
    select(EXP,Sample_ID,Pregnancy_ID,Origin,Condition)


#aa$Sample_ID[which(!aa$Sample_ID %in% cc$Sample_ID[which(cc$Study =="Myometrium" )])]

cell.counts <- aa %>% group_by(EXP,Sample_ID,Pregnancy_ID,Origin,Condition) %>%
    summarize(n=n()) 


cc.wider <- cell.counts %>% ungroup() %>% select(EXP,Sample_ID,Pregnancy_ID,Origin,Condition,n) %>%
    group_by(EXP,Sample_ID,Pregnancy_ID,Origin,Condition) %>%
    pivot_wider(names_from=Origin,values_from=n,values_fill=0,names_prefix="n") %>%
    mutate(nT=nM+nF) %>% 
    filter(nT>100) %>% select(-nNA)

#write_tsv(cc.wider,paste0(outFolder,"cc.wider_withsamples.tsv"))
write_tsv(cc.wider,paste0(outFolder,"cc.wider_withsamples_7only.tsv"))


cc.wider <- cell.counts %>% ungroup() %>% select(EXP,Pregnancy_ID,Origin,Condition,n) %>%
    group_by(EXP,Pregnancy_ID,Origin,Condition) %>%
    pivot_wider(names_from=Origin,values_from=n,values_fill=0,names_prefix="n") %>%
    mutate(nT=nM+nF) %>% 
    filter(nT>100) %>% select(-nNA)

#write_tsv(cc.wider,paste0(outFolder,"cc.wider.tsv"))
write_tsv(cc.wider,paste0(outFolder,"cc.wider_7only.tsv"))


#write_tsv(cell.counts,paste0(outFolder,"cell.counts.tsv"))
write_tsv(cell.counts,paste0(outFolder,"cell.counts_7only.tsv"))

