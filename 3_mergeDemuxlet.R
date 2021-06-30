
library(Seurat)
library(Matrix)
library(tidyverse)
rm(list = ls())

#############################################################
# Filtering based on percent.mt, nFeature_RNA
#############################################################

sc<- read_rds("/wsu/home/groups/prbgenomics/labor_myo/myometrium_analysis/2_kb_Output/kb_Seurat.obj.rds")
sc@meta.data$Library <- (gsub("_[ACGT]{6,}","",colnames(sc)))

table(sc$Library)

expNames <- unique(sc@meta.data$Library)
expNames


outFolder= "./3_MergeDemux_Output/"
system(paste0("mkdir -p ",outFolder))
##setwd(outFolder)

future::plan(strategy = 'multicore', workers = 20)
options(future.globals.maxSize = 20 * 1024 ^ 3)




#%>%mutate(Sample_ID=paste(Pregnancy_ID,Origin,sep="-"))
#head(cc)

## Demuxlet output. 
dd2 <- read_rds("./1_demux2_output/1_demux_New.ALL_7only.rds") %>% filter(DIFF.LLK.BEST.NEXT>1) %>%
    select(BARCODE=NEW_BARCODE,SNG.BEST.GUESS,DROPLET.TYPE,NUM.SNPS,NUM.READS,DIFF.LLK.BEST.NEXT,EXP)


dd1 <- read_rds("./1_demux2_output/1_demux_New.ALL.rds") %>% filter(DIFF.LLK.BEST.NEXT>1) %>%
    select(BARCODE=NEW_BARCODE,SNG.BEST.GUESS,DROPLET.TYPE,NUM.SNPS,NUM.READS,DIFF.LLK.BEST.NEXT,EXP)

#dd1 <- dd1 %>% left_join(cc) 


dd<-rbind(dd1,dd2)

table(dd$DROPLET.TYPE)
table(dd$EXP)

md = sc@meta.data %>% rownames_to_column("BARCODE") %>% left_join(dd) %>%
    column_to_rownames("BARCODE")

stopifnot(identical(rownames(md),rownames(sc@meta.data)))



# ## Repeat with soup or cell
dd <- read_rds("./1_souporcell_output/1_souporcell.ALL.rds")

md = md %>% rownames_to_column("barcode") %>% left_join(dd) %>%
    column_to_rownames("barcode")

md$bc <- rownames(md)


# sample information
cc <- read.csv("/wsu/home/groups/prbgenomics/labor_myo/myometrium_analysis/ParturitionAndMyoSamples.csv",stringsAsFactors = FALSE)%>% mutate(SNG.BEST.GUESS=Sample_ID)
head(cc)
# # 

#cc <- read_tsv("/nfs/rprdata/scilab/labor2/Covid19.Samples.txt") %>%  mutate(SNG.BEST.GUESS=paste(Pregnancy_ID,Origin,sep="-"))

# joint with meta data
md <- md %>% left_join(cc) 
head(md)



sc@meta.data <- md  %>% column_to_rownames("bc")

## Annotate genes with annotables or other and filtering chrM
cmd <- paste0("zcat ",
              "/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_31/gencode.v31.annotation.gff3.gz",
              " | awk '$3~/gene/'",
              " | sed 's/ID=.*gene_id=//;s/;gene_type=/\\t/;s/;gene_name=/\\t/;s/;.*//'")
cat(cmd,"\n")

## Check 0 or 1-based coordinates. 

aux <- data.table::fread(cmd=cmd) %>% mutate(TSS=ifelse(V7=="+",V4,V5)) %>%
    select(Chr=V1,Min=V4,Max=V5,kbid=V9,TSS,Strand=V7,Type=V10,gene_name=V11) 
##anno <- tibble(kbid=rownames(sc)) %>% mutate(ensgene=gsub("\\..*","",kbid)) %>% left_join(grch38)

anno <- tibble(kbid=rownames(sc),rs=rowSums(sc@assays$RNA@data)) %>% filter(rs>0) %>% left_join(aux) %>%
    filter(!is.na(Chr))

table(is.na(anno$Chr))

table(anno$Chr)

table(anno$Type)

head(anno)

head(aux)

write_rds(anno,paste0(outFolder,"anno.rds"))


## Subset to genes in anno 
sc <- sc[anno$kbid,]

## Calculate features 
sc[["percent.mt"]] <- PercentageFeatureSet(sc,features=anno[anno$Chr=="chrM",]$kbid)

sc[["percent.mt"]] %>% summary()

## Filter sc for things matching genotype. chrM or number of RNAs.  

sc[["percent.Y"]] <- PercentageFeatureSet(sc,features=anno[anno$Chr=="chrY",]$kbid)
sc[["percent.Y"]] %>% summary()

##anno[anno$gene_name=="XIST",]$kbid 
sc[["percent.XIST"]] <- PercentageFeatureSet(sc,features=anno[anno$gene_name=="XIST",]$kbid )
sc[["percent.XIST"]] %>% summary()

table(sc[["percent.XIST"]]>0.01,md$SNG.BEST.GUESS)


## Save merged object without any filter.
fname=paste0(outFolder,"scFullSeurat.Rdata")
cat(fname,"\n")
save(sc,anno,file=fname)

##
table(sc[["percent.mt"]]<25,sc@meta.data$nFeature_RNA > 100) 

table(sc@meta.data$Library,sc[["percent.mt"]]<25 & sc@meta.data$DIFF.LLK.BEST.NEXT > 3 & sc@meta.data$nFeature_RNA > 100 & sc$status=="singlet") 

table(sc@meta.data$SNG.BEST.GUESS,sc[["percent.mt"]]<25 & sc@meta.data$DIFF.LLK.BEST.NEXT > 3 & sc@meta.data$nFeature_RNA > 100 & sc$status=="singlet") 

table(sc$status,sc$nFeature_RNA > 8000)

# DIFF.LLK.BEST.NEXT
sc <- subset(sc, subset = nFeature_RNA > 100 & nFeature_RNA < 10000 & DIFF.LLK.BEST.NEXT > 3 & percent.mt < 25)


table(sc@meta.data$Library, sc@meta.data$SNG.BEST.GUESS) 


## Save merged object after filter.
fname=paste0(outFolder,"scFilteredSeurat.Rdata")
cat(fname,"\n")
save(sc,anno,file=fname)




## Run doubletFinder earlier?
### END- HERE ###
########################################################


