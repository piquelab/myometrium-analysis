#############################################################
### cell type classification using SingleR 
### reference: Vento-Tormo et al. Single-cell reconstruction of the early maternal-fetal interface in humans. Nature 563, 347-353 (2018).
#############################################################
library(Seurat)
library(Matrix)
library(tidyverse)
library(future)
library(harmony)
library(SingleR)

future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)


###########################################
outFolder="./4_harmony_cellClass_RVT/"
system(paste0("mkdir -p ", outFolder))

          
# reference
sc2 <- read_rds("/nfs/rprdata/scilab/novogene/otherdata/roser/Analysis/20200226/3_scTransferLabel_RoserPrep/ST_Integrated_RVT.obj.rds")
sc2 <- subset(sc2, subset = nFeature_RNA > 200)


# query: sc 
load("3_MergeDemux_Output/scFilteredSeurat.Rdata")
sc1 <- sc

### Merge
sc <- merge(sc1,list(sc2))
dim(sc)
table(sc$Library)


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

ref.labels <- sc@meta.data$final_cluster[!is.na(sc@meta.data$annotation)]

pred.labels <- SingleR(test = query.he, ref = ref.he, labels = ref.labels)

##table(pred.labels)

table(pred.labels$pruned.labels)

sum(is.na(pred.labels$pruned.labels))

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


