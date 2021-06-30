#############################################################
### cell type classification using SingleR 
### reference: R. Pique-Regi, R. Romero, A. L. Tarca, E. D. Sendler, Y. Xu, V. Garcia-Flores, Y. Leng, F. Luca, S. S. Hassan, N. Gomez-Lopez, Single cell transcriptional signatures of the human placenta in term and preterm parturition. Elife 8,  (2019).
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
outFolder="./4_harmony_cellClass_parturition_elife/"
system(paste0("mkdir -p ", outFolder))





# reference
sc3 <- read_rds("/nfs/rprdata/scilab/novogene/Analyses/Roger_20200218/3_scTransferLabel_scLabor/ST_Integrated.scLabor.obj.rds")
# Subset a Seurat object
sc3 <- subset(sc3, subset = nFeature_RNA > 100)



# query: sc 
load("3_MergeDemux_Output/scFilteredSeurat.Rdata")
sc1 <- sc



### Merge
# to identify the cell types of sc1, another study with known celltypes will be merged to this study and the cell types will be identified

sc <- merge(sc1,list(sc3))
dim(sc)
table(sc$Library)

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
#sc <- FindClusters(sc, verbose = TRUE,resolution=0.6)


################

he <- t(sc@reductions$harmony@cell.embeddings[,1:30])

#unknown cell types
query.he <- he[,is.na(sc@meta.data$FinalName)]

#known cell types
ref.he <- he[,!is.na(sc@meta.data$FinalName)]

ref.labels <- sc@meta.data$FinalName[!is.na(sc@meta.data$FinalName)]

pred.labels <- SingleR(test = query.he, ref = ref.he, labels = ref.labels)

##table(pred.labels)

table(pred.labels$pruned.labels)

sum(is.na(pred.labels$pruned.labels))


fname=paste0(outFolder,"sc.NormByLocation.ref.Harmony.singler.res0.8.rds")
write_rds(pred.labels,fname)

md <- pred.labels %>% as.data.frame() %>% 
    rownames_to_column("BARCODES") %>%
    left_join(sc@meta.data %>% rownames_to_column("BARCODES"))


fname=paste0(outFolder,"sc.NormByLocation.ref.Harmony.singler.res0.8.csv")
write_csv(md,fname)

## save object.

fname=paste0(outFolder,"sc.NormIntegrated.ref.Harmony.res0.8.rds")
write_rds(sc,fname)

