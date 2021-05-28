library(MuSiC)
library(Seurat)
library(Matrix)
library(tidyverse)
library(xbioc)
library(tidyverse)
library(dplyr)
library(stringr)
library(Biobase)
library(reshape)
library(cowplot)
library(pheatmap)

#outFolder <- paste0("./14_deconvolution_analysis_music/")
outFolder <- paste0("./14_deconvolution_analysis_music/merge_subcelltypes")
system(paste0("mkdir -p ",outFolder))



# MuSiC uses two types of input data:Bulk expression obtained from RNA sequencing, Multi-subject single cell expression obtained from single-cell RNA sequencing (scRNA-seq)

# # example esets
# EMTAB.eset = readRDS('14_deconvolution_analysis_musci/EMTABesethealthy.rds')
# XinT2D.eset=readRDS("14_deconvolution_analysis_musci/XinT2Deset.rds")
# XinT2D.eset@phenoData@data$cellType
# GSE50244.bulk.eset = readRDS('14_deconvolution_analysis_musci/GSE50244bulkeset.rds')




# construct ExpressionSet for bulk and single cell data 

# ene_exprs.matrix is your gene expression matrix (genes by cells)
# pheno.matrix is a data frame of phenotype annotation (rownames must match the column names of gene_exprs.matrix). 
#metadata <- data.frame(labelDescription= c("Sample ID", "Subject Name", "Cell Type ID", "Cell Type Name"), row.names=c("sampleID", "SubjectName", "cellTypeID", "cellType"))
#SC.eset = ExpressionSet(assayData = data.matrix(gene_exprs.matrix), phenoData =  new("AnnotatedDataFrame", data = pheno.matrix, varMetadata = metadata) )

#ExpressionSet(assayData, phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE), featureData=annotatedDataFrameFrom(assayData, byrow=TRUE), experimentData=MIAME(), annotation=character(), protocolData=annotatedDataFrameFrom(assayData, byrow=FALSE), ...)

############################################################
# bulk 
load("TLTNLmyoToCaseWest.rdata")
eset_bulk<-eset
colnames(eset_bulk)<-unlist(strsplit(colnames(eset),"_"))[seq(2,2*length(colnames(eset)),by=2)]
rownames(eset_bulk)<-unlist(strsplit(rownames(eset),"_"))[seq(2,2*length(rownames(eset)),by=2)]

tl<-table(rownames(eset_bulk))
nonrepeat<-eset_bulk[rownames(eset_bulk) %in% names(tl)[tl==1],]
repeats<-sapply(names(tl)[tl>1], function(x) {
  print(x)
  mat<-eset_bulk[which(rownames(eset_bulk)==x),]
  mat[which(rowSums(mat)==max(rowSums(mat))),]
  
})

repeats<-t(repeats)
rownames(repeats)<-names(tl)[tl>1]
bulk_matrix<-rbind(nonrepeat,repeats)
colnames(bulk_matrix)<-colnames(eset)

bulk_matrix<-apply(bulk_matrix,c(1,2),function(x) {return(2 ^ as.numeric(x))})


metadata <- data.frame(labelDescription= c("SubjectName", "Disease","sampleID"), row.names=c("SubjectName", "Disease","sampleID"))
gene_exprs.matrix<-bulk_matrix
Disease<-unlist(strsplit(colnames(eset),"_"))[seq(2,2*length(colnames(eset)),by=2)]
SubjectName<-unlist(strsplit(colnames(eset),"_"))[seq(1,2*length(colnames(eset)),by=2)]
sampleID<-colnames(bulk_matrix)
pheno.matrix<-cbind(Disease,SubjectName,sampleID)
rownames(pheno.matrix)<-colnames(bulk_matrix)
bulk.eset = ExpressionSet(assayData = data.matrix(gene_exprs.matrix), phenoData =  new("AnnotatedDataFrame", data = as.data.frame(pheno.matrix), varMetadata = metadata) )


############################################################
#single cell 
# sc2 <- read_rds("4_harmony/sc.NormByLibrary.Harmony.StringentFiltering.res0.8.rds")
# clust2Name<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
#               "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast",
#               "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
# clust2Name<-paste0(c(0:23),"_",clust2Name)
# names(clust2Name)<-c(0:23)

# merging sub cell types
clust2Name<-c("Stromal","Macrophage","Macrophage","Endothelial","Monocyte","T-cell","Decidual","T-cell","LED","Stromal","ILC","NK-cell","Smooth muscle cells","Myofibroblast","Macrophage","Endothelial","DC","Smooth muscle cells","EVT","Plasmablast","Smooth muscle cells","Macrophage","B-cell","Unciliated Epithelial")
#clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
#clust2Names<-paste0(c(0:23),"_",clust2Names)
names(clust2Name)<-c(0:23)




load("3_MergeDemux_Output/scFilteredSeurat.Rdata")
data<-sc@assays$RNA@data
data<-data[anno$kbid,]
rownames(data)<-anno$gene_name
cl<-colnames(data)
cl<-cl[which(cl%in% names(sc2$seurat_clusters) )]
data<-data[,cl]

tl<-table(rownames(data))
nonrepeat<-data[rownames(data) %in% names(tl)[tl==1],]
repeats<-sapply(names(tl)[tl>1], function(x) {
  print(x)
  mat<-data[which(rownames(data)==x),]
  mat[which(rowSums(mat)==max(rowSums(mat)))[1],]
  
})

repeats<-t(repeats)
sc_matrix<-rbind(nonrepeat,repeats)
colnames(sc_matrix)<-colnames(data)

# sc_matrix[which(rownames(sc_matrix)=="RF00006"),]

metadata <- data.frame(labelDescription= c("sampleID", "CellTypeID", "CellTypeName"), row.names=c("sampleID", "CellTypeID", "CellTypeName"))
gene_exprs.matrix<-sc_matrix

CellTypeID<-names(clust2Name[sc2$seurat_clusters[cl]])
CellTypeName<-as.character(clust2Name[sc2$seurat_clusters[cl]])
# colnames  
# libraries<-sc2@meta.data$Pregnancy_ID #sc@meta.data$Library
# names(libraries)<-rownames(sc2@meta.data)
# libraries<-libraries[colnames(gene_exprs.matrix)]
# SampleID<- libraries
sampleID<- colnames(sc_matrix)
  
pheno.matrix<-cbind(sampleID,CellTypeID,CellTypeName)
rownames(pheno.matrix)<-colnames(gene_exprs.matrix)
sc.eset = ExpressionSet(assayData = data.matrix(gene_exprs.matrix), phenoData =  new("AnnotatedDataFrame", data = as.data.frame(pheno.matrix), varMetadata = metadata) )



# Estimate cell type proportions
Est.prop.myo = music_prop(bulk.eset = bulk.eset, sc.eset = sc.eset, clusters = 'CellTypeName',samples = 'sampleID', verbose = T)
names(Est.prop.myo)
save(Est.prop.myo,file=paste0(outFolder,"Est.prop.myo.RData"))





# prop.weighted<-Est.prop.myo$Est.prop.weighted
# rownames(prop.weighted)<-unlist(strsplit(rownames(prop.weighted),"_"))[seq(2,2*length(rownames(prop.weighted)),by=2)]
# 
# prop.weighted<-t(prop.weighted)
# 
# 
# 
# 
# 
# pvalues_celltypes<-sapply(rownames(prop.weighted), function(x){
#   pvalue=t.test(prop.weighted[x,colnames(prop.weighted)=="TL"],prop.weighted[x,colnames(prop.weighted)=="TNL"])$p.value 
#   pvalue
#   })
# 
# pvalues_celltypes<-cbind(pvalues_celltypes,p.adjust(pvalues_celltypes,"fdr"))
# colnames(pvalues_celltypes)<-c("pvalue","pvalue.fdr")
# 
# write.csv(pvalues_celltypes,paste0(outFolder,"prop.weighted_pvalues.csv"))
# 
# 
# 
# 
# prop.allgene<-Est.prop.myo$Est.prop.allgene
# rownames(prop.allgene)<-unlist(strsplit(rownames(prop.allgene),"_"))[seq(2,2*length(rownames(prop.allgene)),by=2)]
# 
# prop.allgene<-t(prop.allgene)
# 
# pvalues_celltypes<-sapply(rownames(prop.allgene), function(x){
#   pvalue=t.test(prop.allgene[x,colnames(prop.allgene)=="TL"],prop.allgene[x,colnames(prop.allgene)=="TNL"])$p.value 
#   pvalue
# })
# 
# pvalues_celltypes<-cbind(pvalues_celltypes,p.adjust(pvalues_celltypes,"fdr"))
# colnames(pvalues_celltypes)<-c("pvalue","pvalue.fdr")
# 
# write.csv(pvalues_celltypes,paste0(outFolder,"prop.allgene_pvalues.csv"))
# 
# 
# 
# 
# fname=paste0(outFolder,"heatmap_proportion_music.pdf");
# pdf(fname,width=7,height=7)
# pheatmap(Est.prop.myo$Est.prop.weighted,cluster_rows=TRUE,scale="none")
# dev.off()
# 
# 
# 
