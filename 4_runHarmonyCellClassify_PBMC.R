# install.packages("Seurat")
# remotes::install_github("mojaveazure/seurat-disk")
# install.packages('remotes')
# 

# remotes::install_version(package = 'Seurat', version = package_version('4.0.0'))
# install.packages("Seurat", version='4.0.0')



outFolder="./4_harmony_cellClass_PBMC/"
system(paste0("mkdir -p ", outFolder))


library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)


##############################################
# loading the PBMC reference 
##############################################
#library(reticulate)
#reference <- LoadH5Seurat("Seurat_PBMC_multimodal/pbmc_multimodal.h5seurat")
#save(reference,file="Seurat_PBMC_multimodal/pbmc_multimodal.RData")
load("Seurat_PBMC_multimodal/pbmc_multimodal.RData")


DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()



# seurat query example
#library(SeuratData)
# InstallData('pbmc3k')
# data("pbmc3k")
#pbmc3k <- SCTransform(pbmc3k, verbose = TRUE)



##############################################
# loading the query
##############################################

load("3_MergeDemux_Output/scFilteredSeurat.Rdata")
sc1 <- sc
DefaultAssay(sc1) <- "RNA"
 






##############################################
# Converting ensamble to symbol (qurery)
##############################################



library(tidyverse)
## Annotate genes with annotables or other and filtering chrM?
cmd <- paste0("zcat ",
              "/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_31/gencode.v31.annotation.gff3.gz",
              " | awk '$3~/gene/'",
              " | sed 's/ID=.*gene_id=//;s/;gene_type=/\\t/;s/;gene_name=/\\t/;s/;.*//'")
cat(cmd,"\n")

## Check 0 or 1-based coordinates.

aux <- data.table::fread(cmd=cmd) %>% mutate(TSS=ifelse(V7=="+",V4,V5)) %>%
    select(Chr=V1,Min=V4,Max=V5,kbid=V9,TSS,Strand=V7,Type=V10,gene_name=V11)
##anno <- tibble(kbid=rownames(sc)) %>% mutate(ensgene=gsub("\\..*","",kbid)) %>% left_join(grch38)

anno <- tibble(kbid=rownames(sc1),rs=rowSums(sc1@assays$RNA@data)) %>% filter(rs>0) %>% left_join(aux) %>%
    filter(!is.na(Chr))



# RenameGenesSeurat <- function(obj, newnames ) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
#     print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
#     RNA <- obj@assays$RNA
#     
#     if (nrow(RNA) == length(newnames)) {
#         if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
#         if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
#         if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
#     } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
#     obj@assays$RNA <- RNA
#     return(obj)
# }


sc1 <- sc1[anno$kbid,]

cn<-sc1@assays$RNA@counts
md<-sc1@meta.data

cn<-cn[anno$kbid,]
rownames(cn)<-anno$gene_name

selectedgene <- intersect(rownames(cn),rownames(reference))
cn<-cn[selectedgene,]
newquery <- CreateSeuratObject(counts = cn, meta.data = md)



#sc1<-RenameGenesSeurat(sc1,newnames=anno$gene_name)
sc1 <- SCTransform(newquery, verbose = TRUE)

## Subset to genes in anno





# sc1<-sc1[selectedgene,]
# reference<-reference[selectedgene,]

# sc1@assays$RNA@data<-sc1@assays$RNA@data[selectedgene,]
# reference@assays$SCT@data<-reference@assays$SCT@data[selectedgene,]



# 
 
#  sc1 <- FindVariableFeatures(sc1, selection.method = "vst")
#  sc1 <- ScaleData(sc1, verbose = TRUE) 
#  sc1 <- RunPCA(sc1,pc.genes = sc1@var.genes, npcs = 100, verbose = TRUE)
# # DefaultAssay(sc1) <- "RNA"
#  sc1 <- RunUMAP(object = sc1, reduction = "pca", dims = 1:30)
#  
#  sc1 <- FindNeighbors(object = sc1, reduction = "pca", dims = 1:30)
#  
#  sc1 <- FindClusters(sc1, resolution=0.8)
 


##############################################
#  finding anchors between reference and query
##############################################


anchors <- FindTransferAnchors(
    reference = reference,
    query = sc1, #pbmc3k, 
    normalization.method = "SCT",
    reference.reduction = "spca",
    recompute.residuals=FALSE,
    dims = 1:50)


##############################################
#transfer cell type labels and protein data from the reference to the query
##############################################
sc1 <- MapQuery(
    anchorset = anchors,
    query = sc1,
    reference = reference,
    refdata = list(
        celltype.l1 = "celltype.l1",
        celltype.l2 = "celltype.l2",
        predicted_ADT = "ADT"
    ),
    reference.reduction = "spca", 
    reduction.model = "wnn.umap"
)

#Explore the mapping results

fname=paste0(outFolder,"UMAP_predicted.celltypes.l1.l2.png");
png(fname,width=1000,height=1000)
p1 = DimPlot(sc1, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(sc1, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
p1 + p2
dev.off()


# predicted.celltype.l1<-sc1$predicted.celltype.l1
# predicted.celltype.l2<-sc1$predicted.celltype.l2




predicted.celltype<-sc1@meta.data %>% select(predicted.celltype.l1,predicted.celltype.l2,predicted.celltype.l1.score,predicted.celltype.l2.score)
    
    
fname=paste0(outFolder,"sc.NormByLocation.ref.Anchors.rds")
write_rds(predicted.celltype,fname)

#ref
 md <- predicted.celltype %>% as.data.frame() %>%
    rownames_to_column("BARCODES") %>%
     left_join(sc1@meta.data %>% rownames_to_column("BARCODES"))

 fname=paste0(outFolder,"sc.NormByLocation.ref.Anchors.csv")
 write_csv(md,fname)

 
fname=paste0(outFolder,"sc.Integrated.Harmony.rds")
write_rds(sc1,fname)



