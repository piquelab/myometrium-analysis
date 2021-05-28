
library(tidyverse)
library(dplyr)
library(Seurat)
library(Matrix)
library(tidyverse)
library(tidyverse)
library(dplyr)
library(stringr)
library(DESeq2)
library(enrichplot)
library(dplyr)
library(org.Hs.eg.db)
library(ReactomePA)
library(clusterProfiler)
set.seed(0)

outFolder <- paste0("./14_deconvolution_analysis/merge_subcelltypes/")
system(paste0("mkdir -p ",outFolder))



anno <- read_rds("3_MergeDemux_Output/anno.rds")
# gene_symbol <- anno$gene_name
# names(gene_symbol) <- anno$kbid

ref_data <- read.delim("myometrium_term_TL-TNL_ALLList.txt")
ref_data<-ref_data %>% dplyr::select(SYMBOL,logFC,P.Value,adj.P.Val,ENTREZ,t )
colnames(ref_data)<-c("R.gene_name","R.Log2FC","Rpvalue","Rpadj","ENTREZID","Rt")
ref_data <- ref_data %>% filter(!is.na(R.Log2FC) & !is.na(ENTREZID)  & !is.na(Rpadj))
ref_data$ENTREZID<-as.character(ref_data$ENTREZID)

selected_genes<-ref_data %>% filter(abs(R.Log2FC)>=0.4 & Rpadj<=0.064) %>% dplyr::select ("R.gene_name") %>% unlist
selected_genes<-unique(selected_genes)


write.table(selected_genes, file=paste0(outFolder,'selected_1000genes.txt'), quote=FALSE, sep='\t',col.names = FALSE,row.names = FALSE)


selected_genes<-ref_data %>% filter(Rpadj<=0.3) %>% dplyr::select ("R.gene_name") %>% unlist

write.table(selected_genes, file=paste0(outFolder,'selected_6448genes.txt'), quote=FALSE, sep='\t',col.names = FALSE,row.names = FALSE)


selected_genes<-ref_data %>% dplyr::select ("R.gene_name") %>% unlist

write.table(selected_genes, file=paste0(outFolder,'selected_genes_all.txt'), quote=FALSE, sep='\t',col.names = FALSE,row.names = FALSE)


######################################################
# bulk data
######################################################
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

#write_rds(eset_bulk,file=paste0(outFolder,"eset_bulk.rds"))

#bulk_matrix<-as.matrix(eset_bulk)
#write_rds(bulk_matrix,file=paste0(outFolder,"bulk_matrix.rds"))

#bulk_matrix<-bulk_matrix[rw,]






bulk_matrix<-cbind(rownames(bulk_matrix), bulk_matrix)
colnames(bulk_matrix)[1]<-"Gene symbol"
bulk_matrix<-rbind(colnames(bulk_matrix), bulk_matrix)


 rownames(bulk_matrix)<-NULL
 colnames(bulk_matrix)<-NULL
 bulk_matrix[-1,-1]<-apply(bulk_matrix[-1,-1],c(1,2),function(x) {return(2 ^ as.numeric(x))})
 
#write.csv(bulk_matrix,file=paste0(outFolder,"bulk_matrix.csv"))
#write.table(bulk_matrix, file=paste0(outFolder,'bulk_matrix.txt'), quote=FALSE, sep='\t',col.names = FALSE,row.names = FALSE)
#write_tsv(as.data.frame(bulk_matrix), file=paste0(outFolder,'bulk_matrix.tsv'),col_names = FALSE )

#bulk_matrix[,1][which(bulk_matrix[,1] %in% c('ABCF2', 'AHRR', 'ARMCX5-GPRASP2', 'ATXN7', 'C2orf27A', 'CCDC39', 'DIABLO', 'GGT1', 'GOLGA8M', 'HSPA14', 'ITFG2-AS1', 'LINC01238', 'MATR3', 'PDE11A', 'PINX1', 'POLR2J3', 'POLR2J4', 'RF00003', 'RF00004', 'RF00006', 'RF00012', 'RF00015', 'RF00017', 'RF00019', 'RF00026', 'RF00045', 'RF00056', 'RF00072', 'RF00090', 'RF00091', 'RF00096', 'RF00100', 'RF00156', 'RF00191', 'RF00322', 'RF00402', 'RF00409', 'RF00410', 'RF00416', 'RF00422', 'RF00425', 'RF00432', 'RF00554', 'RF00560', 'RF00561', 'RF00568', 'RF00586', 'RF00601', 'RF00614', 'RF01210', 'RF01225', 'RF01241', 'RF02106', 'RF02271', 'SOD2', 'TBCE', 'TMSB15B' ))]

#which(table(bulk_matrix[,1])>1)


write.table(bulk_matrix, file=paste0(outFolder,'bulk_matrix_exp.txt'), quote=FALSE, sep='\t',col.names = FALSE,row.names = FALSE)


#bulk_matrix2<-read_delim(paste0(outFolder,"bulk_matrix_exp.txt"),delim = "\t")


######################################################
# single cell data
######################################################

######################################################
# # single cell signature matrix 
# 
# m2 = read_tsv("5_harmonyClustersDGE/ClusterDEG.tsv")
# 
# m3 <- m2 %>% filter(p_val_adj<0.1,avg_log2FC>0.5) %>%
#   group_by(gene) %>%
#   mutate(H=log2(length(cluster))) %>%
#   filter(H<=1) %>%
#   ungroup()
# 
# table(m3$cluster)
# 
# m3$cluster<-clust2Name[as.character(m3$cluster)]
# top20 <- m3 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) %>% ungroup()
# 
# top20<-top20 %>% select(cluster,symbol,avg_log2FC)
# 
# singlecell_signature<-matrix(-1,nrow=length(unique(top20$symbol)),ncol=length(unique(top20$cluster)))
# rownames(singlecell_signature)<-unique(top20$symbol)
# colnames(singlecell_signature)<-unique(top20$cluster)
# 
# # this could not be used- check the format of signature matrix in example below
# write.table(top20, file='singlecell.singature.txt', quote=FALSE, sep='\t', row.names = TRUE, col.names = TRUE)
# 
# # example
# single_cell_signature_example_LM22<-read_delim(paste0(outFolder,"single_cell_signature_example_LM22.txt"),delim = "\t")

####################################
# single cell signature // another way average
# provide the average gene expression per cell-type calculated manually after normalization

sc2 <- read_rds("4_harmony/sc.NormByLibrary.Harmony.StringentFiltering.res0.8.rds")


# # all cell types
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
#sc <- NormalizeData(sc, verbose=TRUE) 
data<-sc@assays$RNA@data
data<-data[anno$kbid,]
rownames(data)<-anno$gene_name
cl<-colnames(data)
cl<-cl[which(cl%in% names(sc2$seurat_clusters) )]
data<-data[,cl]
colnames(data)<-clust2Name[sc2$seurat_clusters[cl]] #


#rn_singlecell<-rownames(data)
#rn_bulk<-bulk_matrix[-1,1]
#rw<-intersect(rn_bulk,rn_singlecell)

#data<-data[rw,]


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
data<-sc_matrix
# average per cell type

#data<-apply(data,c(1,2),function(x) {return(2 ^ as.numeric(x))})


mdata<-lapply(unique(colnames(data)),function(x){
  d<-data[,which(colnames(data)==x)]
  dat<-apply(d,1,mean)
  dat

  })

singlecell.singature.avg <- do.call(cbind,mdata)





singlecell.singature.avg<-cbind(rownames(singlecell.singature.avg), singlecell.singature.avg)
colnames(singlecell.singature.avg)<-c("Gene symbol",unique(colnames(data)))
#write.table(singlecell.singature.avg, file=paste0(outFolder,'singlecell.singature.avg.txt'), quote=FALSE, sep='\t', col.names = TRUE)
#singlecell.singature.avg<-as.data.frame(singlecell.singature.avg)
#write_tsv(singlecell.singature.avg,file=paste0(outFolder,"singlecell.singature.avg.tsv")) #,delim = "\t"



#rownames(singlecell.singature.avg)[which(rownames(singlecell.singature.avg) %in% c('ABCF2', 'AHRR', 'ARMCX5-GPRASP2', 'ATXN7', 'C2orf27A', 'CCDC39', 'DIABLO', 'GGT1', 'GOLGA8M', 'HSPA14', 'ITFG2-AS1', 'LINC01238', 'MATR3', 'PDE11A', 'PINX1', 'POLR2J3', 'POLR2J4', 'RF00003', 'RF00004', 'RF00006', 'RF00012', 'RF00015', 'RF00017', 'RF00019', 'RF00026', 'RF00045', 'RF00056', 'RF00072', 'RF00090', 'RF00091', 'RF00096', 'RF00100', 'RF00156', 'RF00191', 'RF00322', 'RF00402', 'RF00409', 'RF00410', 'RF00416', 'RF00422', 'RF00425', 'RF00432', 'RF00554', 'RF00560', 'RF00561', 'RF00568', 'RF00586', 'RF00601', 'RF00614', 'RF01210', 'RF01225', 'RF01241', 'RF02106', 'RF02271', 'SOD2', 'TBCE', 'TMSB15B' ))]


write.table(singlecell.singature.avg, file=paste0(outFolder,'singlecell.singature.avg_exp.txt'), quote=FALSE, sep='\t', col.names = TRUE)




anno <- read_rds("3_MergeDemux_Output/anno.rds")
sc <- read_rds("4_harmony/sc.NormByLibrary.Harmony.StringentFiltering.res0.8.rds")

clust2Name<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
              "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast",
              "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Name<-paste0(c(0:23),"_",clust2Name)
names(clust2Name)<-c(0:23)



# md <- read_rds("./4_harmony_cellClass_PBMC/sc.NormByLocation.ref.Anchors.rds") %>%
#   as.data.frame %>%
#   rownames_to_column("BARCODES") %>%
#   select(BARCODES,scLabor_ID=predicted.celltype.l2,scLabor_Score=predicted.celltype.l2.score)
# md <- sc@meta.data %>% rownames_to_column("BARCODES") %>%
#   left_join(md) 
# identical(md$BARCODES,rownames(sc@meta.data))

######################################################
# count data
counts<-sc@assays$RNA@counts
rw<-anno$gene_name[which(anno$kbid %in% rownames(counts))]
names(rw)<-anno$kbid[which(anno$kbid %in% rownames(counts))]
counts<-counts[names(rw),]
rownames(counts)<-as.character(rw)
colnames(counts)<-as.character(clust2Name[sc$seurat_clusters[colnames(counts)]])

write.table(counts, file='singlecell.scale.counts.txt', quote=FALSE, sep='\t', col.names = TRUE,row.names = TRUE)

counts<-as.data.frame(counts)
counts<-cbind(rownames(counts), counts)
colnames(counts)[1]<-"GeneSymbol"
cl<-colnames(counts)
counts<-rbind(cl, counts)


write.csv(counts,file=paste0(outFolder,"singlecell.scale.counts.csv"),col.names = cl)
write.table(counts, file='singlecell.scale.counts.txt', quote=FALSE, sep='\t', col.names = TRUE,row.names = TRUE)


###########################################
#random sampled  


table(sc$seurat_clusters)
random_samples<-sapply(unique(sc$seurat_clusters),function(x){
  index<-names(sc$seurat_clusters) [which( sc$seurat_clusters==x)]
  ln<-min(length(index),1000)
  random_samples_celltype<-index[sample(ln)]
  random_samples_celltype
})

random_samples<-unlist(random_samples)


counts<-sc@assays$RNA@counts
counts<-counts[,random_samples]

rw<-anno$gene_name[which(anno$kbid %in% rownames(counts))]
names(rw)<-anno$kbid[which(anno$kbid %in% rownames(counts))]
counts<-counts[names(rw),]
rownames(counts)<-as.character(rw)
colnames(counts)<-as.character(clust2Name[sc$seurat_clusters[colnames(counts)]])

write.table(counts, file='singlecell.random_cells.counts.txt', quote=FALSE, sep='\t', col.names = TRUE,row.names = TRUE)


counts<-as.data.frame(counts)
counts<-cbind(rownames(counts), counts)
colnames(counts)[1]<-"GeneSymbol"
cl<-colnames(counts)
counts<-rbind(cl, counts)

######################################################
# # scaled data
# 
# scale.data<-sc@assays$RNA@scale.data
# rw<-anno$gene_name[which(anno$kbid %in% rownames(scale.data))]
# names(rw)<-anno$kbid[which(anno$kbid %in% rownames(scale.data))]
# scale.data<-scale.data[names(rw),]
# rownames(scale.data)<-as.character(rw)
# colnames(scale.data)<-clust2Name[sc$seurat_clusters[colnames(scale.data)]]
# 
# scale.data<-cbind(rownames(scale.data), scale.data)
# colnames(scale.data)[1]<-"GeneSymbol"
# cl<-colnames(scale.data)
# scale.data<-rbind(cl, scale.data)
# rownames(scale.data)<-NULL
# colnames(scale.data)<-NULL
# 
# 
# 
# write.csv(counts,file=paste0(outFolder,"singlecell.scale.data.csv"),col.names = cl)
# write.table(counts, file='singlecell.scale.data.txt', quote=FALSE, sep='\t', row.names = TRUE, col.names = TRUE)
# 




########################
# bulk<-as.matrix(bulk_matrix[-1,-1])
# rownames(bulk)<-bulk_matrix[-1,1]
# 
# bulk_mean<-apply(bulk,1,function(x) mean(as.numeric(x)))
# single_cell_mean<-as.matrix(singlecell.singature.avg[-1,-1])
# rownames(single_cell_mean)<-singlecell.singature.avg[-1,1]
# 
# 
# data<-cbind(single_cell_mean,bulk_mean)
# 
# library(pheatmap)
# data<-apply(data,c(1,2),as.numeric)
# cor_matrix<-cor(data)
# 
# fname=paste0(outFolder,"heatmap_bulk_singlecell.pdf");
# pdf(fname,width=7,height=7)
# # paletteLength<-30
# # my_palette <- colorRampPalette(colors = c("#333399", "white", "#A50021"))(n = paletteLength)
# # myBreaks <- c(seq(min(cor_matrix), 0, length.out=ceiling(paletteLength/2) + 1), seq(max(cor_matrix)/paletteLength, max(cor_matrix), length.out=floor(paletteLength/2)))
# pheatmap(cor_matrix,cluster_rows=TRUE,scale="none")#,breaks=myBreaks)
# dev.off()
# 
# 
# dens <- density(single_cell)
# # plot density
# plot(dens, frame = FALSE, , main = "Density plot- single cell data")
# hist(single_cell, col = "steelblue", frame = FALSE, breaks = 30, main = "Histogram- sinlge cell data")


########################################################################################################################
# cibersortx output
########################################################################################################################

# outFolder <- paste0("cibersortx_results/CIBERSORTx_Job12_output/")
# 
# system(paste0("cat cibersortx_results/CIBERSORTx_Job12_output/README_Group.txt"))
# list.files(outFolder)
# 
# Fractions<-read_delim(paste0(outFolder,"CIBERSORTxGEP_Job12_Fractions.txt"),delim = "\t")
# GEP<-read_delim(paste0(outFolder,"CIBERSORTxGEP_Job12_GEPs.txt"),delim = "\t")
# GEP_Weights<-read_delim(paste0(outFolder,"CIBERSORTxGEP_Job12_Weights.txt"),delim = "\t")
# #SM_GEPs_Filtered<-read_delim(paste0(outFolder,"CIBERSORTxGEP_Job12_SM_GEPs_Filtered.txt"),delim = "\t")
# GEPs_Filtered<-read_delim(paste0(outFolder,"CIBERSORTxGEP_Job12_GEPs_Filtered.txt"),delim = "\t")
# 
# GEP_lognorm<-GEP
# GEP_lognorm[-1,-1]<-apply(GEP[-1,-1],c(1,2),function(x)return(log2(x)))



#######################
# job 15: high resolution
#######################


# file enumerating the fractions of the different cell types in bulks samples. 

# *_Fractions.txt:
#CIBERSORTxGEP_Job15_Fractions.txt


# CIBERSORTx Group Mode imputes representative cell type-specific expression 
# profiles.In doing so, it generates a set of regression coefficients that 
# represent the average expression value of each gene within each cell type across
# the set (i.e., "group") of input mixture files.
#*_GEPs.txt
#CIBERSORTxGEP_Job15_GEPs.txt



# high resolution results
#CIBERSORTxHiRes_Job15_Plasmablast_Window20.txt // and other cell types

HiResfiles <- list.files(paste0(outFolder,"CIBERSORTx_Job17_output/"),pattern="*_Window20.txt")

#HiResfiles <- list.files(paste0(outFolder,"CIBERSORTx_Job15_output/"),pattern="*_Window20.txt")

HiResfiles<-HiResfiles[c(-6,-8)]

resList<-lapply(HiResfiles,function(x)
  {
  

  celltype_gexp_samples<-read_delim(paste0("14_deconvolution_analysis/merge_subcelltypes/CIBERSORTx_Job17_output/",x),delim = "\t")
  celltype<-unlist(strsplit(x,"_"))[3]
  rw<-celltype_gexp_samples$GeneSymbol
  cl<-colnames(celltype_gexp_samples)[-1]
  celltype_gexp_samples<-as.matrix(celltype_gexp_samples)
  celltype_gexp_samples<-celltype_gexp_samples[,-1]
  rownames(celltype_gexp_samples)<-rw
  colnames(celltype_gexp_samples)<-cl
  conditions<-unlist(strsplit(cl,"_"))[seq(2,2*length(cl),by=2)]
  samples<-colnames(celltype_gexp_samples)
  cvt<-data.frame(Indiv=samples, Group=conditions)
  
  cleaned<-na.omit(celltype_gexp_samples)
  
  cleaned<-apply(cleaned, c(1,2),as.numeric)
 

  cat("########## ",celltype,"\n")
  
  if(!all(cleaned==1))
  {
    
    filterrow<-apply(cleaned, 1, function(x) if (all(x==1)) return (TRUE) else return(FALSE))
    cleaned<-cleaned[!filterrow,]
    if( nrow(cleaned)>0)
    {
      
      #cleaned<-apply(cleaned, c(1,2),log2)
      dds <- DESeqDataSetFromMatrix(round(cleaned),cvt, ~ Group)    
      dds <- DESeq(dds,parallel=TRUE)
      res <- results(dds)
      #write_tsv(dds,paste0(outFolder,))
      res$celltype<-rep(celltype,nrow(res))
      res$symbol<-rownames(res)
      as.data.frame(res)
    }
   
  }
  
  
  })

outFolder<-paste0(outFolder,"CIBERSORTx_Job17_DEGs/")
#outFolder<-paste0(outFolder,"CIBERSORTx_Job15_DEGs/")
system(paste0("mkdir -p ",outFolder))

res <- do.call(rbind,resList)

colnames(res)[8]<-"gene_name"
res<-as.data.frame(res)
             
   
eg = bitr(res$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
names(eg)[1]="gene_name"
head(eg)
e2g <- eg$gene_name
names(e2g) <- eg$ENTREZID
res <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))




res %>% filter(padj<0.1,abs(log2FoldChange)>0) %>% dplyr::count(celltype) %>%
  write_tsv(paste0(outFolder,"Summary.FDR.tsv"))

res %>% write_tsv(paste0(outFolder,"ALL.combined.tsv"))

res %>% filter(padj<0.1,abs(log2FoldChange)>0.0) %>%
  write_tsv(paste0(outFolder,"SIG.combined.tsv"))

res_deconBulk<-res %>% filter(padj<0.1,abs(log2FoldChange)>0) 

####################################################################
#  single vs bulk deconvoluted
####################################################################


res <- read_tsv("7_outputs_DESeq_ConditionsByCluster/SIG.combined.2021-02-17.tsv")


eg = bitr(res$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
names(eg)[1]="gene_name"
head(eg)
e2g <- eg$gene_name
names(e2g) <- eg$ENTREZID
res <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))



res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
res <- res %>% filter(!is.na(pvalue))
clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(c(0:23),"_",clust2Names)
names(clust2Names)<-c(0:23)

res$Cell_type<-clust2Names[res$Cell_type]

res_deconvBulk_SMC<-res_deconBulk %>% filter (celltype== "Smoothmusclecells" )
res_single_cell_SMC<-res %>% filter (Cell_type== "12_Smooth muscle cells-1" )

length(which(res_deconvBulk_SMC$ENTREZID %in% res_single_cell_SMC$ENTREZID))
length(which(!res_deconvBulk_SMC$ENTREZID %in% res_single_cell_SMC$ENTREZID))
length(which(! res_single_cell_SMC$ENTREZID %in% res_deconvBulk_SMC$ENTREZID))



##############################################################################
# ORA
##############################################################################
genes <- unique(c(res_deconvBulk_SMC$ENTREZID, res_single_cell_SMC$ENTREZID))



genes <- unique(c(res_deconBulk$ENTREZID, res_single_cell_SMC$ENTREZID))



geneUniv <- res %>% dplyr::select(ENTREZID) %>% unlist %>% unique
geneUniv<-unique(c(geneUniv,genes))


ego<- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP")
res_df_enrichGO<-ego@result
res_df_enrichGO<-res_df_enrichGO%>% filter(p.adjust<0.1)
res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})

pdf(paste0(outFolder,"enrichGO_decov_combined_singlecell_SMC_DotPlot.pdf"),width=10,height=4)
ggplot(res_df_enrichGO, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Description)) +
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  #scale_colour_gradient(limits=c(0, 0.10), low="red") +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
  ylab(NULL)+
  xlab("GeneRatio")
dev.off()


enrichKEGG.res <- enrichKEGG(gene=genes,universe=geneUniv,organism="hsa")

res_df_enrichKEGG<-enrichKEGG.res@result

res_df_enrichKEGG<-res_df_enrichKEGG%>% filter(qvalue<0.1)


# res_df_enrichKEGG$GeneRatio<-sapply(res_df_enrichKEGG$GeneRatio, function(x){
#   numden<-unlist(strsplit(x,"/"))
#   return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# enrichRPath.res <- enrichPathway(gene=genes,universe=geneUniv)
# 
# res_dfRPath<-enrichRPath.res@result
# res_dfRPath<-res_dfRPath%>% filter(qvalue<0.1)
# res_dfRPath$GeneRatio<-sapply(res_dfRPath$GeneRatio, function(x){
#   numden<-unlist(strsplit(x,"/"))
#   return (as.numeric(numden[1])/as.numeric(numden[2]))
# })





pdf(paste0(outFolder,"enrichRPath_combined_singlecell_SMC_DotPlot.pdf"),width=10,height=4)
ggplot(res_dfRPath, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Description)) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  #scale_colour_gradient(limits=c(0, 0.10), low="red") +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
  ylab(NULL)+ 
  xlab("GeneRatio") 
dev.off()




library(magrittr)
library(clusterProfiler)

gene<-genes
#gene <- filter(res,padj<0.1,abs(log2FoldChange)>=0.5) %>% dplyr::select(ENTREZID) %>% unlist %>% unique

wp2gene <- read.gmt("13_sample_investigation_plots/wikipathways-20191210-gmt-Homo_sapiens.gmt")
#wp2gene <-read.gmt(wpgmtfile)

wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME





#"IGFBP5"  "GUCY1A1"      "JUN" 

ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
head(ewp)


res_dfewp<-ewp@result
res_dfewp<-res_dfewp%>% filter(qvalue<0.1)
res_dfewp$GeneRatio<-sapply(res_dfewp$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})



pdf(paste0(outFolder,"enrichWiki_decov_combined_singlecell_SMC_DotPlot.pdf"),width=10,height=10)
ggplot(res_dfewp, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Description)) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  #scale_colour_gradient(limits=c(0, 0.10), low="red") +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
  ylab(NULL)+ 
  xlab("GeneRatio") 
dev.off()



res_dfewp$Description[1]

genes_Myometrialpathway<-unlist(strsplit(res_dfewp$geneID[1],"/"))

eg$gene_name[which(eg$ENTREZID %in% genes_Myometrialpathway)]

 
##############################################################################
# GSEA
##############################################################################
res <- read_tsv(paste0(outFolder,"ALL.combined.tsv"))
res<-res %>% filter(celltype=="Smoothmusclecells")
geneList <- -log10(res$pvalue)
names(geneList) <- res$ENTREZID
geneList = sort(geneList, decreasing = TRUE)

gseGO.res <- gseGO(geneList,  OrgDb=org.Hs.eg.db,ont="BP")
res_df_gseGO<-gseGO.res@result%>% filter(qvalues<=0.1)

message("gsePathway")
gseRPath.res <- gsePathway(geneList,pvalueCutoff = 1)
print(head(gseRPath.res))

res_df<-gseRPath.res@result
which(res_df$ID =="R-HSA-445355")


gseKEGG.res <-gseKEGG( geneList)
res_df<-gseKEGG.res@result %>% filter(qvalues<=0.1) #[1:10,]

