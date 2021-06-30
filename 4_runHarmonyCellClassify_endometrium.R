#############################################################
### cell type classification using SingleR 
### reference: W. Wang, F. Vilella, P. Alama, I. Moreno, M. Mignardi, A. Isakova, W. Pan, C. Simon, S. R. Quake, Single-cell transcriptomic atlas of the human endometrium during the menstrual cycle. Nat Med 26, 1644-1653 (2020).

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
outFolder="./4_harmony_cellClass_endometrium/"
system(paste0("mkdir -p ", outFolder))

#reference: sc,anno
load("/nfs/rprdata/scilab/endometrium/otherdata/GSE111976/gs10x/Analysis/4_harmony_cellClass/scFullSeurat.Rdata")
sc2<-sc
sc2 <- subset(sc2, subset = nFeature_RNA > 200)



#query: sc 
load("3_MergeDemux_Output/scFilteredSeurat.Rdata")
sc1 <- sc

### Merge
# to identify the cell types of sc1, another study with known celltypes will be merged to this study and the cell types will be identified
sc <- merge(sc1,list(sc2))


dim(sc)

table(sc$Library)

## Harmony

DefaultAssay(sc) <- "RNA"

sc <- NormalizeData(sc, verbose=TRUE)

sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 3000)

sc <- ScaleData(sc, verbose = TRUE)

sc <- RunPCA(sc,pc.genes = sc@var.genes, npcs = 100, verbose = TRUE)

sc <- RunHarmony(sc,c("Library"),reduction="pca")

sc <- RunUMAP(sc,reduction = "harmony", dims = 1:30)


###### Cluster

sc <- FindNeighbors(sc, reduction = "harmony", dims = 1:30, verbose = TRUE)

sc <- FindClusters(sc, verbose = TRUE,resolution=0.8)

################

he <- t(sc@reductions$harmony@cell.embeddings[,1:30])

#unknown cell types
#query.he <- he[,is.na(sc@meta.data$FinalName)]
query.he <- he[,is.na(sc@meta.data$cell_type)]

#known cell types
#ref.he <- he[,!is.na(sc@meta.data$FinalName)]
ref.he <- he[,!is.na(sc@meta.data$cell_type)]


ref.labels <- sc@meta.data$cell_type[!is.na(sc@meta.data$cell_type)]

pred.labels <- SingleR(test = query.he, ref = ref.he, labels = ref.labels)

##table(pred.labels)

table(pred.labels$pruned.labels)

sum(is.na(pred.labels$pruned.labels))


fname=paste0(outFolder,"sc.NormByLocation.ref.Harmony.singler.rds")
write_rds(pred.labels,fname)

#ref
md <- pred.labels %>% as.data.frame() %>%
    rownames_to_column("BARCODES") %>%
    left_join(sc@meta.data %>% rownames_to_column("BARCODES"))


fname=paste0(outFolder,"sc.NormByLocation.ref.Harmony.singler.csv")
write_csv(md,fname)


## save object
fname=paste0(outFolder,"sc.NormByLibFullIntegrated.refendometrium.Harmony.res0.8.rds")
write_rds(sc,fname)

