## options(repos = c(CRAN = "http://cran.rstudio.com"))
##   This uses updated Seurat package 3 - starts with merged counts/demux from step 2

library(Seurat)
library(Matrix)
library(tidyverse)

library(future)

library(harmony)

#################
library(SingleR)



future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)


###########################################
## Testing sc transform           
## 2and3_Diem_Output
## adata <- read_rds("./kb_diem_Output/kb_diem_Seurat.list.rds")

load("3_MergeDemux_Output/scFilteredSeurat.Rdata")
sc1 <- sc

#ref
#/nfs/rprdata/scilab/endometrium/otherdata/RVT-2021/endometrium_all.h5ad
#


#this does not work
#sc3<-ReadH5AD("/nfs/rprdata/scilab/endometrium/otherdata/RVT-2021/endometrium_all.h5ad")

library(reticulate)
#pathsc3<-Convert("/nfs/rprdata/scilab/endometrium/otherdata/RVT-2021/endometrium_all.h5ad", dest="h5seurat")
pathsc3<-"/nfs/rprdata/scilab/endometrium/otherdata/RVT-2021/endometrium_all.h5seurat"
sc3 <- LoadH5Seurat(pathsc3)


#sc3
# Subset a Seurat object
sc3 <- subset(sc3, subset = nFeature_RNA > 100)
#sc3<-subset(sc3, subset = Cell.type != "Other")
sc3<-sc3[ ,sc3$Cell.type!="Other"]

##sc1 <- read_rds("/nfs/rprdata/scilab/novogene/otherdata/roser/Analysis/20200226/3_scTransferLabel_scLabor/ST_Integrated.scLabor.obj.rds")
##sc1 <- subset(sc1, subset = nFeature_RNA > 100)



outFolder="./4_harmony_cellClass_RVT-2021/"
system(paste0("mkdir -p ", outFolder))
##setwd(outFolder)

### Merge
# to identify the cell types of sc1, another study with known celltypes will be merged to this study and the cell types will be identified



##sc <- merge(sc1,list(sc2,sc3,sc4))

sc <- merge(sc1,list(sc3))


dim(sc)

table(sc$Library)

## table(sc$Location) 

##table(sc$sclabor.tlabel)

##table(sc$Location,sc$sclabor.tlabel)

## Harmony

DefaultAssay(sc) <- "RNA"

sc <- NormalizeData(sc, verbose=TRUE) 

sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 3000)

sc <- ScaleData(sc, verbose = TRUE) 

sc <- RunPCA(sc,pc.genes = sc@var.genes, npcs = 100, verbose = TRUE)

##sc <- RunHarmony(sc,c("Location","percent.mt","Rep"),reduction="pca")
sc <- RunHarmony(sc,c("Library"),reduction="pca")

sc <- RunUMAP(sc,reduction = "harmony", dims = 1:30)


###### Cluster

sc <- FindNeighbors(sc, reduction = "harmony", dims = 1:30, verbose = TRUE)

sc <- FindClusters(sc, verbose = TRUE,resolution=0.8)



################

he <- t(sc@reductions$harmony@cell.embeddings[,1:30])

#unknown cell types
query.he <- he[,is.na(sc@meta.data$Cell.type)]

#known cell types
ref.he <- he[,!is.na(sc@meta.data$Cell.type)]

ref.labels <- sc@meta.data$Cell.type[!is.na(sc@meta.data$Cell.type)]

pred.labels <- SingleR(test = query.he, ref = ref.he, labels = ref.labels)

##table(pred.labels)

table(pred.labels$pruned.labels)

sum(is.na(pred.labels$pruned.labels))

##sum(is.na(pred.labels$labels))

##md <- sc@meta.data %>% select(scLaborLabs=pruned.labels)


fname=paste0(outFolder,"sc.NormByLocation.ref.Harmony.singler.res0.8.rds")
write_rds(pred.labels,fname)

md <- pred.labels %>% as.data.frame() %>% 
    rownames_to_column("BARCODES") %>%
    left_join(sc@meta.data %>% rownames_to_column("BARCODES"))


fname=paste0(outFolder,"sc.NormByLocation.ref.Harmony.singler.res0.8.csv")
write_csv(md,fname)

## save object.



#sc<-read_rds("/4_harmony_cellClass/sc.NormIntegrated.ref.Harmony.rds")

 fname=paste0(outFolder,"sc.NormIntegrated.ref.Harmony.res0.8.rds")
 
 write_rds(sc,fname)




### END- HERE ###
########################################################

