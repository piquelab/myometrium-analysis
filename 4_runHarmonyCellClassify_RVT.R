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


##sc3 <- read_rds("/nfs/rprdata/scilab/novogene/Analyses/Roger_20200218/3_scTransferLabel_scLabor/ST_Integrated.scLabor.obj.rds")
##sc3 <- subset(sc3, subset = nFeature_RNA > 100)
sc2 <- read_rds("/nfs/rprdata/scilab/novogene/otherdata/roser/Analysis/20200226/3_scTransferLabel_RoserPrep/ST_Integrated_RVT.obj.rds")
##sc2 <- read_rds("/nfs/rprdata/scilab/novogene/otherdata/roser/Analysis/20200226/3_scTransferLabel_scLabor/")
sc2 <- subset(sc2, subset = nFeature_RNA > 200)

outFolder="./4_harmony_cellClass_RVT/"
system(paste0("mkdir -p ", outFolder))
##setwd(outFolder)

### Merge



##sc <- merge(sc1,list(sc2,sc3,sc4))

sc <- merge(sc1,list(sc2))


dim(sc)

table(sc$Library)

## table(sc$Location) 

## table(sc$sclabor.tlabel)

## table(sc$Location,sc$sclabor.tlabel)

## Harmony

DefaultAssay(sc) <- "RNA"

sc <- NormalizeData(sc, verbose=TRUE) 

sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 3000)

sc <- ScaleData(sc, verbose = TRUE) 

sc <- RunPCA(sc,pc.genes = sc@var.genes, npcs = 100, verbose = TRUE)

##sc <- RunHarmony(sc,c("Location","percent.mt","Rep"),reduction="pca")

sc <- RunHarmony(sc,c("Library"),reduction="pca",seed=TRUE)

sc <- RunUMAP(sc,reduction = "harmony", dims = 1:30)


###### Cluster

sc <- FindNeighbors(sc, reduction = "harmony", dims = 1:30, verbose = TRUE)

sc <- FindClusters(sc, verbose = TRUE, resolution=0.8)



################

he <- t(sc@reductions$harmony@cell.embeddings[,1:30])

query.he <- he[,is.na(sc@meta.data$annotation)]

ref.he <- he[,!is.na(sc@meta.data$annotation)]

##sc@meta.data$annotation

##ref.labels <- sc@meta.data$annotation[!is.na(sc@meta.data$annotation)]

ref.labels <- sc@meta.data$final_cluster[!is.na(sc@meta.data$annotation)]

pred.labels <- SingleR(test = query.he, ref = ref.he, labels = ref.labels)

##table(pred.labels)

table(pred.labels$pruned.labels)

sum(is.na(pred.labels$pruned.labels))

##sum(is.na(pred.labels$labels))

##md <- sc@meta.data %>% select(scLaborLabs=pruned.labels)


fname=paste0(outFolder,"sc.NormByLibrary.ref.Harmony.singler.RVT.res0.8.rds")
write_rds(pred.labels,fname)

md <- pred.labels %>% as.data.frame() %>% 
    rownames_to_column("BARCODES") %>%
    left_join(sc1@meta.data %>% rownames_to_column("BARCODES"))


fname=paste0(outFolder,"sc.NormByLibrary.refRVT.Harmony.singler.res0.8.csv")
write_csv(md,fname)

## save object.
fname=paste0(outFolder,"sc.NormByLibFullIntegrated.refRVT.Harmony.res0.8.rds")
write_rds(sc,fname)




### END- HERE ###
########################################################


