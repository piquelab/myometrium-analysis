#############################################################
### Deconvolution analysis using Cibersortx

#############################################################

library(limma)
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
library(magrittr)


#############################################################
# Cibersortx Input data input preparation
#############################################################


set.seed(0)
system(paste0("mkdir -p ",outFolder))



# bulk data- log2fc  
ref_data <- read.delim("myometrium_term_TL-TNL_ALLList.txt")
ref_data<-ref_data %>% dplyr::select(SYMBOL,logFC,P.Value,adj.P.Val,ENTREZ,t )
colnames(ref_data)<-c("R.gene_name","R.Log2FC","Rpvalue","Rpadj","ENTREZID","Rt")
ref_data <- ref_data %>% filter(!is.na(R.Log2FC) & !is.na(ENTREZID)  & !is.na(Rpadj))
ref_data$ENTREZID<-as.character(ref_data$ENTREZID)



######################################################
# 1. "Mixture" input file for cibersortx

# bulk data- samples  
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



bulk_matrix<-cbind(rownames(bulk_matrix), bulk_matrix)
colnames(bulk_matrix)[1]<-"Gene symbol"
bulk_matrix<-rbind(colnames(bulk_matrix), bulk_matrix)
rownames(bulk_matrix)<-NULL
colnames(bulk_matrix)<-NULL
bulk_matrix[-1,-1]<-apply(bulk_matrix[-1,-1],c(1,2),function(x) {return(2 ^ as.numeric(x))})
 write.table(bulk_matrix, file=paste0(outFolder,'bulk_matrix_exp.txt'), quote=FALSE, sep='\t',col.names = FALSE,row.names = FALSE)




######################################################
# 2. "Signature Genes" input file for cibersortx
 
# single cell data
 
######################################################

# single cell data 
anno <- read_rds("3_MergeDemux_Output/anno.rds")
sc2 <- read_rds("4_harmony/sc.NormByLibrary.Harmony.StringentFiltering.res0.8.rds")

# merging sub cell types to one main cell type 
clust2Name<-c("Stromal","Macrophage","Macrophage","Endothelial","Monocyte","T-cell","Decidual","T-cell","LED","Stromal","ILC","NK-cell","Smooth muscle cells","Myofibroblast","Macrophage","Endothelial","DC","Smooth muscle cells","EVT","Plasmablast","Smooth muscle cells","Macrophage","B-cell","Unciliated Epithelial")
names(clust2Name)<-c(0:23)


# single cell signature 
load("3_MergeDemux_Output/scFilteredSeurat.Rdata")
data<-sc@assays$RNA@data
data<-data[anno$kbid,]
rownames(data)<-anno$gene_name
cl<-colnames(data)
cl<-cl[which(cl%in% names(sc2$seurat_clusters) )]
data<-data[,cl]
colnames(data)<-clust2Name[sc2$seurat_clusters[cl]] #


tl<-table(rownames(data))
nonrepeat<-data[rownames(data) %in% names(tl)[tl==1],]


# The row with maximum exp is selected if there are multiple replicates for a gene
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


# Average based on cell types

mdata<-lapply(unique(colnames(data)),function(x){
  d<-data[,which(colnames(data)==x)]
  dat<-apply(d,1,mean)
  dat

  })

singlecell.singature.avg <- do.call(cbind,mdata)


singlecell.singature.avg<-cbind(rownames(singlecell.singature.avg), singlecell.singature.avg)
colnames(singlecell.singature.avg)<-c("Gene symbol",unique(colnames(data)))

write.table(singlecell.singature.avg, file=paste0(outFolder,'singlecell.singature.avg_exp.txt'), quote=FALSE, sep='\t', col.names = TRUE)

singlecell.singature.avg<-read.delim("14_deconvolution_analysis/merge_subcelltypes/singlecell.singature.avg_exp.txt")


anno <- read_rds("3_MergeDemux_Output/anno.rds")
sc <- read_rds("4_harmony/sc.NormByLibrary.Harmony.StringentFiltering.res0.8.rds")

clust2Name<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
              "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast",
              "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Name<-paste0(c(0:23),"_",clust2Name)
names(clust2Name)<-c(0:23)




##################################################################################################################################
# Processing Cibersortx output
# High resolution results 

##################################################################################################################################


outFolder <- paste0("./14_deconvolution_analysis/merge_subcelltypes/")
HiResfiles <- list.files(paste0(outFolder,"CIBERSORTx_Job17_output/"),pattern="*_Window20.txt")

HiResfiles<-HiResfiles[c(-6,-8)]

# processing the cibersortx output

# The "1" values in the expression matrix txt files are genes with insufficient evidence of expression (these genes are either not expressed or have inadequate statistical power to be imputed).
# The NA values are genes that have inadequate statistical power to be imputed.


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
      
      cleaned<-apply(cleaned, c(1,2),log2)
      pDat1<-matrix(conditions,nrow=ncol(cleaned),ncol=1)
      colnames(pDat1)<-"title"
      rownames(pDat1)<-colnames(cleaned)
      aDFrame <- new("AnnotatedDataFrame",  data = as.data.frame(pDat1))
      rownames(aDFrame)<-colnames(cleaned)

      #deseq
      ##cleaned<-apply(cleaned, c(1,2),log2)
      dds <- DESeqDataSetFromMatrix(round(cleaned),cvt, ~ Group)
      dds <- DESeq(dds,parallel=TRUE)
      res <- results(dds)
      res$celltype<-rep(celltype,nrow(res))
      res$symbol<-rownames(res)
      as.data.frame(res)
    }
   
  }
  
  
  })




########################################################
# Differentially expressed genes 
########################################################
outFolder<-paste0("./14_deconvolution_analysis/merge_subcelltypes/CIBERSORTx_Job17_limma_DEGs/")
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

colnames(res)[which(colnames(res)=="adj.P.Val")]<-"padj"
colnames(res)[which(colnames(res)=="logFC")]<-"log2FoldChange"
colnames(res)[which(colnames(res)=="P.Value")]<-"pvalue"

save(res,file=paste0(outFolder,"ALL.combined.RData"))   


res %>% filter(padj<0.1,abs(log2FoldChange)>0) %>% dplyr::count(celltype) %>%
  write_tsv(paste0(outFolder,"Summary.FDR.tsv"))

res %>% write_tsv(paste0(outFolder,"ALL.combined.tsv"))

res %>% filter(padj<0.1,abs(log2FoldChange)>0.0) %>%
  write_tsv(paste0(outFolder,"SIG.combined.tsv"))

res_deconBulk<-res %>% filter(padj<0.1,abs(log2FoldChange)>0) 

####################################################################
#  single vs bulk deconvoluted
####################################################################


load(paste0("14_deconvolution_analysis/merge_subcelltypes/genes_Myometrialpathway_decov_all_singlecell_SMC1_ORAwiki_myogenes.RData"))

outFolder<-paste0("./14_deconvolution_analysis/merge_subcelltypes/CIBERSORTx_Job17_limma_DEGs/")

res <- read_tsv(paste0(outFolder,"ALL.combined.tsv"))

res_deconBulk<-res %>% filter(padj<0.1,abs(log2FoldChange)>0) 

res_deconBulk_genes<-res_deconBulk %>% filter (ENTREZID %in% genes_Myometrialpathway) 
write.csv(res_deconBulk_genes,file=paste0("14_deconvolution_analysis/merge_subcelltypes/CIBERSORTx_Job17_deseq_DEGs/genes_bulk-deconv_wikiMyometrialpathway.csv"))

#res <- read_tsv("7_outputs_DESeq_ConditionsByCluster/SIG.combined.2021-02-17.tsv")
res <- read_tsv("7_outputs_DESeq_ConditionsByCluster_bath_library/SIG.combined.2021-10-18.tsv")

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

res_singlecell_genes<-res %>% filter (ENTREZID %in% genes_Myometrialpathway) 
write.csv(res_singlecell_genes,file=paste0("14_deconvolution_analysis/merge_subcelltypes/CIBERSORTx_Job17_deseq_DEGs/genes_singlecell_wikiMyometrialpathway.csv"))

res_deconvBulk_SMC<-res_deconBulk %>% filter (celltype== "Smoothmusclecells" )
res_single_cell_SMC<-res %>% filter (Cell_type== "12_Smooth muscle cells-1" )

length(which(res_deconvBulk_SMC$ENTREZID %in% res_single_cell_SMC$ENTREZID))
length(which(!res_deconvBulk_SMC$ENTREZID %in% res_single_cell_SMC$ENTREZID))
length(which(! res_single_cell_SMC$ENTREZID %in% res_deconvBulk_SMC$ENTREZID))



##############################################################################
# ORA - single cell and deconv bulk
##############################################################################

#res <- read_tsv("7_outputs_DESeq_ConditionsByCluster/SIG.combined.2021-02-17.tsv")
res <- read_tsv("7_outputs_DESeq_ConditionsByCluster_bath_library/SIG.combined.2021-10-18.tsv")

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
res_single_cell_SMC<-res %>% filter (Cell_type== "12_Smooth muscle cells-1" )


outFolder<-paste0("./14_deconvolution_analysis/merge_subcelltypes/CIBERSORTx_Job17_deseq_DEGs/")

res <- read_tsv(paste0(outFolder,"ALL.combined.tsv"))
res_deconBulk<-res %>% filter(padj<0.1,abs(log2FoldChange)>0) 
res_deconvBulk_SMC<-res_deconBulk %>% filter (celltype== "Smoothmusclecells" )


res_deconBulk<-res %>% filter(padj<0.1,abs(log2FoldChange)>0) 
genes_deconvBulk_SMC<-res_deconvBulk_SMC %>%filter(padj <0.1 & abs(log2FoldChange)>0.5)%>% dplyr::select(gene_name) %>% unlist %>% unique


genes <- unique(c(res_deconBulk$ENTREZID, res_single_cell_SMC$ENTREZID))
geneUniv <- res %>% dplyr::select(ENTREZID) %>% unlist %>% unique
geneUniv<-unique(c(geneUniv,genes))

################################
# ORA GO
################################
ego<- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP")
res_df_enrichGO<-ego@result
res_df_enrichGO<-res_df_enrichGO%>% filter(p.adjust<0.1)
res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})

outFolder<-"14_deconvolution_analysis_plots_batch_correction/"
system(paste0("mkdir -p ",outFolder))


pdf(paste0(outFolder,"enrichGO_decov_combined_singlecell_SMC_DotPlot.pdf"),width=10,height=10)
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



################################
# ORA KEGG
################################

enrichKEGG.res <- enrichKEGG(gene=genes,universe=geneUniv,organism="hsa")

res_df_enrichKEGG<-enrichKEGG.res@result

res_df_enrichKEGG<-res_df_enrichKEGG%>% filter(qvalue<0.1)


res_df_enrichKEGG$GeneRatio<-sapply(res_df_enrichKEGG$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})


pdf(paste0(outFolder,"enrichKEGG_decov_combined_singlecell_SMC_DotPlot.pdf-1_DotPlot.pdf"),width=6,height=4)
ggplot(res_df_enrichKEGG, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Description)) +
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw() +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab", limit = c(0.00001, 0.1))+
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
  ylab(NULL)+
  xlab("GeneRatio")
dev.off()




################################
# ORA Reactome
################################

enrichRPath.res <- enrichPathway(gene=genes,universe=geneUniv)

res_dfRPath<-enrichRPath.res@result
res_dfRPath<-res_dfRPath%>% filter(qvalue<0.1)
res_dfRPath$GeneRatio<-sapply(res_dfRPath$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})

pdf(paste0(outFolder,"enrichRPath_combined_singlecell_SMC_DotPlot.pdf"),width=10,height=4)
ggplot(res_dfRPath, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Description)) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
  ylab(NULL)+ 
  xlab("GeneRatio") 
dev.off()





################################
# ORA WikiPathways
################################
gene<-genes
#gene <- filter(res,padj<0.1,abs(log2FoldChange)>=0.5) %>% dplyr::select(ENTREZID) %>% unlist %>% unique

wp2gene <- read.gmt("13_sample_investigation_plots/wikipathways-20191210-gmt-Homo_sapiens.gmt")

wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
head(ewp)


res_dfewp<-ewp@result
res_dfewp<-res_dfewp%>% filter(qvalue<0.1)
res_dfewp$GeneRatio<-sapply(res_dfewp$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})


outFolder<-"./14_deconvolution_analysis_plots_batch_correction/"
pdf(paste0(outFolder,"enrichWiki_decov_SMC_singlecell_SMC_DotPlot_0.05.pdf"),width=10,height=10)
res_dfewp<-res_dfewp %>% filter(qvalue<0.1)
ggplot(res_dfewp, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Description)) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
  ylab(NULL)+ 
  xlab("GeneRatio") 
dev.off()



res_dfewp$Description[1]

genes_Myometrialpathway<-unlist(strsplit(res_dfewp$geneID[which(res_dfewp$Description=="Myometrial Relaxation and Contraction Pathways")],"/"))
save(genes_Myometrialpathway,file=paste0(outFolder,"genes_Myometrialpathway_decov_all_singlecell_SMC1_ORAwiki_myogenes.RData"))




##############################################################################
# ORA- Myofibroblast
##############################################################################

outFolder<-paste0("./14_deconvolution_analysis/merge_subcelltypes/CIBERSORTx_Job17_limma_DEGs/")

res <- read_tsv(paste0(outFolder,"ALL.combined.tsv"))

res_deconBulk<-res %>% filter(padj<0.1,abs(log2FoldChange)>0) 
res_deconvBulk_Myofibroblast<-res_deconBulk %>% filter (celltype== "Myofibroblast" )

genes <- as.character(unique(res_deconvBulk_Myofibroblast$ENTREZID))
geneUniv <- as.character(res %>% dplyr::select(ENTREZID) %>% unlist %>% unique)


################################
# ORA WikiPathways
################################


ego<- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP")
res_df_enrichGO<-ego@result
res_df_enrichGO<-res_df_enrichGO%>% filter(p.adjust<0.1)
res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})

pdf(paste0(outFolder,"enrichGO_Myofibroblas_DotPlot.pdf"),width=10,height=10)

ggplot(res_df_enrichGO,
       aes(x = GeneRatio, y = Description)) +
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
  ylab(NULL)+
  xlab("GeneRatio")
dev.off()



################################
# ORA KEGG
################################


enrichKEGG.res <- enrichKEGG(gene=genes,universe=geneUniv,organism="hsa")
res_df_enrichKEGG<-enrichKEGG.res@result

res_df_enrichKEGG<-res_df_enrichKEGG%>% filter(qvalue<0.1)

res_df_enrichKEGG$GeneRatio<-sapply(res_df_enrichKEGG$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})


pdf(paste0(outFolder,"enrichKEGG_Myofibroblas_DotPlot.pdf"),width=6,height=4)
ggplot(res_df_enrichKEGG, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Description)) +
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw() +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab", limit = c(0.00001, 0.1))+
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
  ylab(NULL)+
  xlab("GeneRatio")
dev.off()



################################
# ORA Reactome
################################

enrichRPath.res <- enrichPathway(gene=genes,universe=geneUniv)

res_dfRPath<-enrichRPath.res@result
res_dfRPath<-res_dfRPath%>% filter(qvalue<0.1)
res_dfRPath$GeneRatio<-sapply(res_dfRPath$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})


pdf(paste0(outFolder,"enrichRPath_Myofibroblas_DotPlot.pdf"),width=10,height=4)
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




################################
# ORA WikiPathways
################################


gene<-genes
wp2gene <- read.gmt("13_sample_investigation_plots/wikipathways-20191210-gmt-Homo_sapiens.gmt")


wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME



ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
head(ewp)


res_dfewp<-ewp@result
res_dfewp<-res_dfewp%>% filter(qvalue<0.1)
res_dfewp$GeneRatio<-sapply(res_dfewp$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})


pdf(paste0(outFolder,"enrichWiki_decov_SMC_singlecell_SMC_DotPlot_0.05.pdf"),width=10,height=10)
#pdf(paste0(outFolder,"enrichWiki_decov_combined_singlecell_SMC_DotPlot_0.05.pdf"),width=10,height=10)
ggplot(res_dfewp, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Description)) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
  ylab(NULL)+ 
  xlab("GeneRatio") 
dev.off()



##############################################################################
# ORA- SMC
##############################################################################

outFolder<-paste0("./14_deconvolution_analysis/merge_subcelltypes/CIBERSORTx_Job17_limma_DEGs/")

res <- read_tsv(paste0(outFolder,"ALL.combined.tsv"))

res_deconBulk<-res %>% filter(padj<0.1,abs(log2FoldChange)>0) 
res_deconvBulk_SMC<-res_deconBulk %>% filter (celltype== "Smoothmusclecells" )

genes <- as.character(unique(res_deconvBulk_SMC$ENTREZID))

geneUniv <- as.character(res %>% dplyr::select(ENTREZID) %>% unlist %>% unique)


################################
# ORA GO
################################

ego<- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP")
res_df_enrichGO<-ego@result
res_df_enrichGO<-res_df_enrichGO%>% filter(p.adjust<0.1)
res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})

pdf(paste0(outFolder,"enrichGO_SMC_DotPlot.pdf"),width=10,height=10)

ggplot(res_df_enrichGO, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Description)) +
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
  ylab(NULL)+
  xlab("GeneRatio")
dev.off()



################################
# ORA KEGG
################################

enrichKEGG.res <- enrichKEGG(gene=genes,universe=geneUniv,organism="hsa")

res_df_enrichKEGG<-enrichKEGG.res@result

res_df_enrichKEGG<-res_df_enrichKEGG%>% filter(qvalue<0.1)


res_df_enrichKEGG$GeneRatio<-sapply(res_df_enrichKEGG$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})


pdf(paste0(outFolder,"enrichKEGG_SMC_DotPlot.pdf"),width=6,height=4)
ggplot(res_df_enrichKEGG, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Description)) +
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw() +
  #scale_colour_gradient(limits=c(0, 0.10), low="red") +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab", limit = c(0.00001, 0.1))+
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
  ylab(NULL)+
  xlab("GeneRatio")
dev.off()


################################
# ORA REACTOME
################################

enrichRPath.res <- enrichPathway(gene=genes,universe=geneUniv)

res_dfRPath<-enrichRPath.res@result
res_dfRPath<-res_dfRPath%>% filter(qvalue<0.1)
res_dfRPath$GeneRatio<-sapply(res_dfRPath$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})

pdf(paste0(outFolder,"enrichRPath_SMC_DotPlot.pdf"),width=10,height=4)
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




################################
# ORA Wikipathways
################################
gene<-genes

wp2gene <- read.gmt("13_sample_investigation_plots/wikipathways-20191210-gmt-Homo_sapiens.gmt")

wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
head(ewp)


res_dfewp<-ewp@result
res_dfewp<-res_dfewp%>% filter(qvalue<0.1)
res_dfewp$GeneRatio<-sapply(res_dfewp$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})


pdf(paste0(outFolder,"enrichWiki_decov_SMC_singlecell_SMC_DotPlot_0.05.pdf"),width=10,height=10)
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



##############################################################################
# GSEA SMC / Myofibroblast
##############################################################################

outFolder<-paste0("./14_deconvolution_analysis/merge_subcelltypes/CIBERSORTx_Job17_limma_DEGs/")

res <- read_tsv(paste0(outFolder,"ALL.combined.tsv"))
#res<-res %>% filter(celltype=="Smoothmusclecells")
res<-res %>% filter(celltype=="Myofibroblast")

res_string<-res %>% filter(padj<=0.05) 

write.csv(res_string,file=paste0(outFolder,"res_smc_limma.csv"))

write.csv(res_string,file=paste0(outFolder,"res_myofibroblast_limma.csv"))

geneList <- -log10(res$pvalue)
names(geneList) <- res$ENTREZID
geneList = sort(geneList, decreasing = TRUE)


##############################################################################
# GSEA GO
##############################################################################


gseGO.res <- gseGO(geneList,  OrgDb=org.Hs.eg.db,ont="BP")
res_df_gseGO<-gseGO.res@result%>% filter(qvalues<=0.1)

fname<-paste0(outFolder,"gseGO_Myofibroblast_DotPlot.png")
p2<-ggplot(res_df_gseGO, 
           aes(x = enrichmentScore, y = Description)) + 
  geom_point(aes(size = enrichmentScore, color = p.adjust)) +
  theme_bw(base_size = 14) +
  #scale_colour_gradient(limits=c(0, 0.10), low="red") +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term" #enrichmentScore
  ylab(NULL)+ 
  theme_bw()
ggsave(fname,p2,width=6,height=5)



##############################################################################
# GSEA Reactome
##############################################################################


message("gsePathway")
gseRPath.res <- gsePathway(geneList,pvalueCutoff = 1)
print(head(gseRPath.res))

gseRPath.res<-gseRPath.res@result %>% filter(qvalues<=0.1) #[1:10,]

fname<-paste0(outFolder,"/gsePathway_Myofibroblast.png")
#fname<-paste0(outFolder,"/gsePathway_SMC_DotPlot.png")
p3<-ggplot(gseRPath.res, # you can replace the numbers to the row number of pathway of your interest
           aes(x = enrichmentScore, y = Description)) + 
  geom_point(aes(size = enrichmentScore, color = p.adjust)) +
  #theme_bw(base_size = 11) +
  #scale_colour_gradient(limits=c(0, 0.10), low="red") +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term"
  ylab(NULL) +
  theme_bw()
ggsave(fname,p3,width=8,height=2)



##############################################################################
# GSEA KEGG
##############################################################################


gseKEGG.res <-gseKEGG( geneList)
res_df<-gseKEGG.res@result %>% filter(qvalues<=0.1) #[1:10,]

fname<-paste0(outFolder,"gseKEGG.res_Myofibroblast_DotPlot.png")
p1<-ggplot(res_df, # you can replace the numbers to the row number of pathway of your interest
           aes(x = enrichmentScore, y = Description)) + 
  geom_point(aes(size = enrichmentScore, color = p.adjust)) +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term"
  ylab(NULL)+
  theme_bw()
ggsave(fname,p1,width=8,height=5)



ewp2 <- GSEA(geneList, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE,pvalueCutoff = 1 )#,minGSSize=5, maxGSSize=3000,eps =0,pvalueCutoff = 1)
head(ewp2)
res_df_gsewiki<-ewp2@result%>% filter(qvalues<=0.1)


##############################################################################
# ORA- Bulk
##############################################################################

outFolder<-paste0("./14_deconvolution_analysis/bulk_enrichment/")
system(paste0("mkdir -p ",outFolder))

anno <- read_rds("3_MergeDemux_Output/anno.rds")

ref_data <- read.delim("myometrium_term_TL-TNL_ALLList.txt")
ref_data<-ref_data %>% dplyr::select(SYMBOL,logFC,P.Value,adj.P.Val,ENTREZ,t )
colnames(ref_data)<-c("R.gene_name","R.Log2FC","Rpvalue","Rpadj","ENTREZID","Rt")
ref_data <- ref_data %>% filter(!is.na(R.Log2FC) & !is.na(ENTREZID)  & !is.na(Rpadj))
ref_data$ENTREZID<-as.character(ref_data$ENTREZID)

res_DE<-ref_data %>% filter(Rpadj<0.1,abs(R.Log2FC)>0) 
genes <- as.character(unique(res_DE$ENTREZID))

geneUniv <- as.character(ref_data %>% dplyr::select(ENTREZID) %>% unlist %>% unique)


##############################################################################
# ORA GO
##############################################################################

ego<- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP")
res_df_enrichGO<-ego@result
res_df_enrichGO<-res_df_enrichGO%>% filter(p.adjust<0.1)
res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})


res_df_enrichGO<-res_df_enrichGO%>% filter(p.adjust<0.05)
res_df_enrichGO<-res_df_enrichGO[1:30,]

pdf(paste0(outFolder,"enrichGO_bulk_DotPlot.pdf"),width=15,height=24)

ggplot(res_df_enrichGO, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Description)) +
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
  labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
  ylab(NULL)+
  xlab("GeneRatio")
dev.off()


##############################################################################
# ORA KEGG
##############################################################################

enrichKEGG.res <- enrichKEGG(gene=genes,universe=geneUniv,organism="hsa")
res_df_enrichKEGG<-enrichKEGG.res@result

res_df_enrichKEGG<-res_df_enrichKEGG%>% filter(qvalue<0.1)


res_df_enrichKEGG$GeneRatio<-sapply(res_df_enrichKEGG$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})


pdf(paste0(outFolder,"enrichKEGG_bulk_DotPlot.pdf"),width=6,height=4)
ggplot(res_df_enrichKEGG, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Description)) +
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw() +
  #scale_colour_gradient(limits=c(0, 0.10), low="red") +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab", limit = c(0.00001, 0.1))+
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
  ylab(NULL)+
  xlab("GeneRatio")
dev.off()


##############################################################################
# ORA Reactome
##############################################################################


enrichRPath.res <- enrichPathway(gene=genes,universe=geneUniv)

res_dfRPath<-enrichRPath.res@result
res_dfRPath<-res_dfRPath%>% filter(qvalue<0.1)
res_dfRPath$GeneRatio<-sapply(res_dfRPath$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})


pdf(paste0(outFolder,"enrichRPath_bulk_DotPlot.pdf"),width=10,height=4)
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




##############################################################################
# ORA WikiPathway
##############################################################################

gene<-genes
wp2gene <- read.gmt("13_sample_investigation_plots/wikipathways-20191210-gmt-Homo_sapiens.gmt")

wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME


ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
head(ewp)


res_dfewp<-ewp@result
res_dfewp<-res_dfewp%>% filter(qvalue<0.1)
res_dfewp$GeneRatio<-sapply(res_dfewp$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})


pdf(paste0(outFolder,"enrichWiki_bulk_DotPlot_0.05.pdf"),width=10,height=10)
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


##############################################################################
# GSEA bulk
##############################################################################


outFolder<-paste0("./14_deconvolution_analysis/bulk_enrichment/")
anno <- read_rds("3_MergeDemux_Output/anno.rds")

ref_data <- read.delim("myometrium_term_TL-TNL_ALLList.txt")
ref_data<-ref_data %>% dplyr::select(SYMBOL,logFC,P.Value,adj.P.Val,ENTREZ,t )
colnames(ref_data)<-c("R.gene_name","R.Log2FC","Rpvalue","Rpadj","ENTREZID","Rt")
ref_data <- ref_data %>% filter(!is.na(R.Log2FC) & !is.na(ENTREZID)  & !is.na(Rpadj))
ref_data$ENTREZID<-as.character(ref_data$ENTREZID)

geneList <- -log10(ref_data$Rpvalue)
names(geneList) <- ref_data$ENTREZID
geneList = sort(geneList, decreasing = TRUE)


##############################################################################
# GSEA GO
##############################################################################


gseGO.res <- gseGO(geneList,  OrgDb=org.Hs.eg.db,ont="BP")
res_df_gseGO<-gseGO.res@result%>% filter(qvalues<=0.1)

res_df_gseGO<-res_df_gseGO[1:30,]
fname<-paste0(outFolder,"gseGO_bulk_DotPlot.png")
p2<-ggplot(res_df_gseGO, 
           aes(x = enrichmentScore, y = Description)) + 
  geom_point(aes(size = enrichmentScore, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
  labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term" #enrichmentScore
  ylab(NULL)+ 
  theme_bw()
ggsave(fname,p2,width=10,height=15)



##############################################################################
# GSEA Reactome
##############################################################################


message("gsePathway")
gseRPath.res <- gsePathway(geneList,pvalueCutoff = 1)
print(head(gseRPath.res))

gseRPath.res<-gseRPath.res@result %>% filter(qvalues<=0.1) #[1:10,]

fname<-paste0(outFolder,"gsePathway_bulk.png")
p3<-ggplot(gseRPath.res, # you can replace the numbers to the row number of pathway of your interest
           aes(x = enrichmentScore, y = Description)) + 
  geom_point(aes(size = enrichmentScore, color = p.adjust)) +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
  labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term"
  ylab(NULL) +
  theme_bw()
ggsave(fname,p3,width=10,height=10)


##############################################################################
# GSEA KEGG
##############################################################################

gseKEGG.res <-gseKEGG( geneList)
res_df<-gseKEGG.res@result %>% filter(qvalues<=0.1) #[1:10,]

res_df<-res_df[1:30,]
fname<-paste0(outFolder,"gseKEGG.bulk_DotPlot.png")
p1<-ggplot(res_df, # you can replace the numbers to the row number of pathway of your interest
           aes(x = enrichmentScore, y = Description)) + 
  geom_point(aes(size = enrichmentScore, color = p.adjust)) +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
  labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term"
  ylab(NULL)+
  theme_bw()
ggsave(fname,p1,width=10,height=10)

ewp2 <- GSEA(geneList, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE,pvalueCutoff = 1 )#,minGSSize=5, maxGSSize=3000,eps =0,pvalueCutoff = 1)
head(ewp2)
res_df_gsewiki<-ewp2@result%>% filter(qvalues<=0.1)



##############################################################################
# Bulk after deconvolution- cell types separately  
##############################################################################


outFolder<-paste0("./14_deconvolution_analysis/merge_subcelltypes/CIBERSORTx_Job17_limma_DEGs/")

res <- read_tsv(paste0(outFolder,"ALL.combined.tsv"))
geneList <- -log10(res$pvalue)
names(geneList) <- res$ENTREZID
geneList = sort(geneList, decreasing = TRUE)



colnames(res)[which(colnames(res)=="celltype")]<-"Cell_type"


pathway_enrich<-function(res_gene=res,cname_select ,padj_cutoff=0.1,log2FoldChange_cutoff=0)
{
  wp2gene <- read.gmt("13_sample_investigation_plots/wikipathways-20191210-gmt-Homo_sapiens.gmt")

  print(cname_select)
  print(cname_select)
  
  aux <- res_gene %>% filter(Cell_type==cname_select)
  
  if(nrow(aux)>0)
  {
    gene<-genes <- filter(aux,padj<padj_cutoff,abs(log2FoldChange)>log2FoldChange_cutoff) %>% dplyr::select(ENTREZID) %>% unlist
    geneUniv <- aux %>% dplyr::select(ENTREZID) %>% unlist
    ##geneList <- aux$log2FoldChange
    geneList <- -log10(aux$pvalue)
    names(geneList) <- aux$ENTREZID
    geneList = sort(geneList, decreasing = TRUE)
    #length(genes)
    
    
    
    
    print(length(gene))
    
    if (length(gene)>0)
    {
      
      
      wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
      wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
      wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
      
      ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
      
      if (!is.null(ewp))
      {
        res_en<-ewp@result
        res_en<-res_en%>%filter(p.adjust<0.1)
        dim1<-nrow(res_en)
        #print(dim1)
        if(min(dim1,10)>0)
        {
          
          
          res_en<-res_en[1:min(dim1,10),c("ID","Description" ,"GeneRatio","p.adjust")]
          res_en$Cell_type<-rep(cname_select,min(dim1,10))
        }
        return (res_en)
      }
      
    }
    
  }
  
}


##############################################################################
# ORA Wikipathway
##############################################################################


res_enrichwikilist<-lapply(unique(res$Cell_type), function(x) {pathway_enrich(cname_select=x)})  
res_df_enrichwiki <- do.call(rbind,res_enrichwikilist)

res_df_enrichwiki$GeneRatio<-sapply(res_df_enrichwiki$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})



pdf(paste0(outFolder,"enrich_wikipathways_cname_DotPlot.pdf"),width=12,height=12)
ggplot(res_df_enrichwiki, # you can replace the numbers to the row number of pathway of your interest
       aes(x = Cell_type, y = Description)) +
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
  ylab(NULL)+
  xlab(NULL)+
  coord_fixed(ratio = 1)+
  #theme_black()+
  theme_bw()+
  theme(text = element_text(size=30)) +
  theme(axis.text.y = element_text(hjust = 1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=10)) 
dev.off()





pathway_GSEA<-function(res_gene=res,cname_select ,padj_cutoff=0.1,log2FoldChange_cutoff=0)
{
  
  wp2gene <- read.gmt("13_sample_investigation_plots/wikipathways-20191210-gmt-Homo_sapiens.gmt")
  print(cname_select)
  aux <- res_gene %>% filter(Cell_type==cname_select)
  
  if(nrow(aux)>0)
  {
    gene<-genes <- filter(aux,padj<padj_cutoff,abs(log2FoldChange)>log2FoldChange_cutoff) %>% dplyr::select(ENTREZID) %>% unlist
    geneUniv <- aux %>% dplyr::select(ENTREZID) %>% unlist
    ##geneList <- aux$log2FoldChange
    geneList <- -log10(aux$pvalue)
    names(geneList) <- aux$ENTREZID
    geneList = sort(geneList, decreasing = TRUE)
    #length(genes)
    
    print(length(gene))
    
    wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
    wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
    wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
    
    
    ewp <- GSEA(geneList, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE,pvalueCutoff = 1 )#,minGSSize=5, maxGSSize=3000,eps =0,pvalueCutoff = 1)
    
    if (!is.null(ewp))
    {
      res_en<-ewp@result
      res_en<-res_en%>%filter(p.adjust<0.1)
      dim1<-nrow(res_en)
      #print(dim1)
      if(min(dim1,10)>0)
      {
        
        
        res_en<-res_en[1:min(dim1,10),c("ID","Description" ,"enrichmentScore","p.adjust")]
        res_en$Cell_type<-rep(cname_select,min(dim1,10))
      }
      return (res_en)
    }
    
  }
  
}

##############################################################################
# GSEA Wikipathway
##############################################################################

res_gseawikilist<-lapply(unique(res$Cell_type), function(x) {pathway_GSEA(cname_select=x)})  
res_df_gseawiki <- do.call(rbind,res_gseawikilist)


pdf(paste0(outFolder,"gse_wikipathways_cname_DotPlot.pdf"),width=7,height=5)
ggplot(res_df_gseawiki, 
       aes(x = Cell_type, y = Description)) + 
  geom_point(aes(size = enrichmentScore, color = p.adjust)) +
  theme_bw(base_size = 14) +
  #scale_colour_gradient(limits=c(0, 0.10), low="red") +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
  labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term" #enrichmentScore
  ylab(NULL)+ 
  #theme_black()+
  theme_bw()+
  #theme(text = element_text(size=30)) +
  theme(axis.text.x = element_text(angle = 45))+
  #xlab(NULL) +
  theme(axis.text.y = element_text(hjust = 1))+
  theme(text = element_text(size=30)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=10)) 
dev.off()




##############################################################################
# ORA GO
##############################################################################


pathway_enrich_go<-function(res_gene=res,cname_select ,padj_cutoff=0.1,log2FoldChange_cutoff=0)
{
  
  print(cname_select)
  print(cname_select)
  
  aux <- res_gene %>% filter(Cell_type==cname_select)
  
  if(nrow(aux)>0)
  {
    gene<-genes <- filter(aux,padj<padj_cutoff,abs(log2FoldChange)>log2FoldChange_cutoff) %>% dplyr::select(ENTREZID) %>% unlist
    genes<-gene<-as.character(genes)
    geneUniv <- as.character(aux %>% dplyr::select(ENTREZID) %>% unlist)
    ##geneList <- aux$log2FoldChange
    geneList <- -log10(aux$pvalue)
     print(length(gene))
    
    if (length(gene)>0)
    {
          message(".................................")
      message("enrichGO")
      ego <- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP")
      print(head(ego))
      
      
      
      
      if (!is.null(ego) & nrow(ego@result)>0 )
      {
        res_en<-ego@result
        res_en<-res_en%>%filter(p.adjust<0.1)
        dim1<-nrow(res_en)
        #print(dim1)
        if(min(dim1,10)>0)
        {
          
          
          res_en<-res_en[1:min(dim1,10),c("ID","Description" ,"GeneRatio","p.adjust")]
          res_en$Cell_type<-rep(cname_select,min(dim1,10))
        }
        return (res_en)
      }
      
    }
    
  }
  
}

res_enrichgo<-lapply(unique(res$Cell_type), function(x) {pathway_enrich_go(cname_select=x)})  
res_df_goi <- do.call(rbind,res_enrichgo)

res_df_goi$GeneRatio<-sapply(res_df_goi$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})


pdf(paste0(outFolder,"enrichGO_cname_DotPlot.pdf"),width=8,height=10)
ggplot(res_df_goi, 
       aes(x = Cell_type, y = Description)) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  #scale_colour_gradient(limits=c(0, 0.10), low="red") +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term" #GeneRatio
  ylab(NULL)+ 
  #theme_black()+
  theme_bw()+
  theme(text = element_text(size=30)) +
  theme(axis.text.x = element_text(angle = 45))+
  #xlab(NULL) +
  #theme(text = element_text(size=30)) +
  theme(axis.text.y = element_text(hjust = 1))+
  theme(text = element_text(size=30)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=10)) 

dev.off()



##############################################################################
# ORA KEGG
##############################################################################


pathway_enrich_kegg<-function(res_gene=res,cname_select ,padj_cutoff=0.1,log2FoldChange_cutoff=0)
{
  
  print(cname_select)
  print(cname_select)
  
  aux <- res_gene %>% filter(Cell_type==cname_select)
  
  if(nrow(aux)>0)
  {
    gene<-genes <- filter(aux,padj<padj_cutoff,abs(log2FoldChange)>log2FoldChange_cutoff) %>% dplyr::select(ENTREZID) %>% unlist
    genes<-gene<-as.character(genes)
    geneUniv <- as.character(aux %>% dplyr::select(ENTREZID) %>% unlist)
    ##geneList <- aux$log2FoldChange
    geneList <- -log10(aux$pvalue)
    print(length(gene))
    
    if (length(gene)>0)
    {
      message(".................................")
      message("enrichGO")
      ego <- enrichKEGG(gene=genes,universe=geneUniv,organism="hsa")
      print(head(ego))
      
      
      
      
      if (!is.null(ego) & nrow(ego@result)>0 )
      {
        res_en<-ego@result
        res_en<-res_en%>%filter(p.adjust<0.1)
        dim1<-nrow(res_en)
        #print(dim1)
        if(min(dim1,20)>0)
        {
          
          
          res_en<-res_en[1:min(dim1,20),c("ID","Description" ,"GeneRatio","p.adjust")]
          res_en$Cell_type<-rep(cname_select,min(dim1,20))
          return (res_en)
        }
        
      }
      
    }
    
  }
  
}

res_enrichkegg<-lapply(unique(res$Cell_type), function(x) {pathway_enrich_kegg(cname_select=x)})  
res_df_kegg <- do.call(rbind,res_enrichkegg)

res_df_kegg$GeneRatio<-sapply(res_df_kegg$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})


pdf(paste0(outFolder,"enrichKEGG_cname_DotPlot.pdf"),width=8,height=10)
ggplot(res_df_kegg, 
       aes(x = Cell_type, y = Description)) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  #scale_colour_gradient(limits=c(0, 0.10), low="red") +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term" #GeneRatio
  ylab(NULL)+ 
  #theme_black()+
  theme_bw()+
  theme(text = element_text(size=30)) +
  theme(axis.text.x = element_text(angle = 45))+
  #xlab(NULL) +
  theme(text = element_text(size=30)) +
  theme(axis.text.y = element_text(hjust = 1))+
  theme(text = element_text(size=30)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=10)) 

dev.off()



##############################################################################
# ORA Reactome
##############################################################################


pathway_enrich_reactome<-function(res_gene=res,cname_select ,padj_cutoff=0.1,log2FoldChange_cutoff=0)
{
  
  print(cname_select)
  print(cname_select)
  
  aux <- res_gene %>% filter(Cell_type==cname_select)
  
  if(nrow(aux)>0)
  {
    gene<-genes <- filter(aux,padj<padj_cutoff,abs(log2FoldChange)>log2FoldChange_cutoff) %>% dplyr::select(ENTREZID) %>% unlist
    genes<-gene<-as.character(genes)
    geneUniv <- as.character(aux %>% dplyr::select(ENTREZID) %>% unlist)
    ##geneList <- aux$log2FoldChange
    geneList <- -log10(aux$pvalue)
    print(length(gene))
    
    if (length(gene)>0)
    {
      message(".................................")
      
      
      message("enrichReactome")
      erpath <- enrichPathway(gene=genes,universe=geneUniv)
      
      print(head(erpath))
      
      
      
      
      if (!is.null(erpath) & nrow(erpath@result)>0 )
      {
        res_en<-erpath@result
        res_en<-res_en%>%filter(p.adjust<0.1)
        dim1<-nrow(res_en)
        #print(dim1)
        if(min(dim1,20)>0)
        {
          
          
          res_en<-res_en[1:min(dim1,20),c("ID","Description" ,"GeneRatio","p.adjust")]
          res_en$Cell_type<-rep(cname_select,min(dim1,20))
          return (res_en)
        }
        
      }
      
    }
    
  }
  
}

res_enrichreacome<-lapply(unique(res$Cell_type), function(x) {pathway_enrich_reactome(cname_select=x)})  
res_df_reacome <- do.call(rbind,res_enrichreacome)

res_df_reacome$GeneRatio<-sapply(res_df_reacome$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})


pdf(paste0(outFolder,"enrichReactome_cname_DotPlot.pdf"),width=15,height=18)
ggplot(res_df_reacome, 
       aes(x = Cell_type, y = Description)) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term" #GeneRatio
  ylab(NULL)+ 
  theme_bw()+
  theme(text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 45))+
  theme(axis.text.y = element_text(hjust = 1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
dev.off()


##############################################################################
# GSEA GO
##############################################################################


pathway_GSEA_go<-function(res_gene=res,cname_select ,padj_cutoff=0.1,log2FoldChange_cutoff=0)
{
  
   print(cname_select)
  
  
  aux <- res_gene %>% filter(Cell_type==cname_select)
  
  if(nrow(aux)>0)
  {
    gene<-genes <- filter(aux,padj<padj_cutoff,abs(log2FoldChange)>log2FoldChange_cutoff) %>% dplyr::select(ENTREZID) %>% unlist
    geneUniv <- aux %>% dplyr::select(ENTREZID) %>% unlist
    ##geneList <- aux$log2FoldChange
    geneList <- -log10(aux$pvalue)
    names(geneList) <- aux$ENTREZID
    geneList = sort(geneList, decreasing = TRUE)
    #length(genes)
    
    print(length(gene))
    
     
      
    gseGO.res <- gseGO(geneList,  OrgDb=org.Hs.eg.db,ont="BP")
    
    
    if (!is.null(gseGO.res))
    {
      res_en<-gseGO.res@result
      res_en<-res_en%>%filter(p.adjust<0.1)
      dim1<-nrow(res_en)
      #print(dim1)
      if(min(dim1,10)>0)
      {
        
        
        res_en<-res_en[1:min(dim1,10),c("ID","Description" ,"enrichmentScore","p.adjust")]
        res_en$Cell_type<-rep(cname_select,min(dim1,10))
        return (res_en)
      }
      
    }
    
  }
  
}

res_gseago<-lapply(unique(res$Cell_type), function(x) {pathway_GSEA_go(cname_select=x)})  
res_df_gsego <- do.call(rbind,res_gseago)


pdf(paste0(outFolder,"gse_go_cname_DotPlot.pdf"),width=10,height=10)
ggplot(res_df_gsego, 
       aes(x = Cell_type, y = Description)) + 
  geom_point(aes(size = enrichmentScore, color = p.adjust)) +
  theme_bw(base_size = 14) +
  #scale_colour_gradient(limits=c(0, 0.10), low="red") +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term" #enrichmentScore
  ylab(NULL)+ 
  #theme_black()+
  theme_bw()+
  theme(text = element_text(size=30)) +
  theme(axis.text.x = element_text(angle = 45))+
  #xlab(NULL) +
  theme(text = element_text(size=30)) +
  theme(axis.text.y = element_text(hjust = 1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=10)) 
dev.off()




##############################################################################
# GSEA KEGG
##############################################################################

pathway_GSEA_kegg<-function(res_gene=res,cname_select ,padj_cutoff=0.1,log2FoldChange_cutoff=0)
{
  
  print(cname_select)
  
  
  aux <- res_gene %>% filter(Cell_type==cname_select)
  
  if(nrow(aux)>0)
  {
    gene<-genes <- filter(aux,padj<padj_cutoff,abs(log2FoldChange)>log2FoldChange_cutoff) %>% dplyr::select(ENTREZID) %>% unlist
    geneUniv <- aux %>% dplyr::select(ENTREZID) %>% unlist
    ##geneList <- aux$log2FoldChange
    geneList <- -log10(aux$pvalue)
    names(geneList) <- aux$ENTREZID
    geneList = sort(geneList, decreasing = TRUE)
    #length(genes)
    
    print(length(gene))
    
    
    
    gsekegg.res <- gseKEGG(geneList,  organism = "hsa")
    
    
    if (!is.null(gsekegg.res))
    {
      res_en<-gsekegg.res@result
      res_en<-res_en%>%filter(p.adjust<0.1)
      dim1<-nrow(res_en)
      #print(dim1)
      if(min(dim1,10)>0)
      {
        
        
        res_en<-res_en[1:min(dim1,10),c("ID","Description" ,"enrichmentScore","p.adjust")]
        res_en$Cell_type<-rep(cname_select,min(dim1,10))
        return (res_en)
      }
      
    }
    
  }
  
}

res_gseakegg<-lapply(unique(res$Cell_type), function(x) {pathway_GSEA_kegg(cname_select=x)})  
res_df_gsekegg <- do.call(rbind,res_gseakegg)


pdf(paste0(outFolder,"gse_kegg_cname_DotPlot.pdf"),width=3,height=2)
ggplot(res_df_gsekegg, 
       aes(x = Cell_type, y = Description)) + 
  geom_point(aes(size = enrichmentScore, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term" #enrichmentScore
  ylab(NULL)+ 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=10)) 
dev.off()



##############################################################################
# GSEA Reactome
##############################################################################


pathway_GSEA_reactome<-function(res_gene=res,cname_select ,padj_cutoff=0.1,log2FoldChange_cutoff=0)
{
  
  print(cname_select)
  
  
  aux <- res_gene %>% filter(Cell_type==cname_select)
  
  if(nrow(aux)>0)
  {
    gene<-genes <- filter(aux,padj<padj_cutoff,abs(log2FoldChange)>log2FoldChange_cutoff) %>% dplyr::select(ENTREZID) %>% unlist
    geneUniv <- aux %>% dplyr::select(ENTREZID) %>% unlist
    ##geneList <- aux$log2FoldChange
    geneList <- -log10(aux$pvalue)
    names(geneList) <- aux$ENTREZID
    geneList = sort(geneList, decreasing = TRUE)
    #length(genes)
    
    print(length(gene))
    
    
    
    gsereacome.res <- gsePathway(geneList)
    
    
    if (!is.null(gsereacome.res))
    {
      res_en<-gsereacome.res@result
      res_en<-res_en%>%filter(p.adjust<0.1)
      dim1<-nrow(res_en)
      #print(dim1)
      if(min(dim1,10)>0)
      {
        
        
        res_en<-res_en[1:min(dim1,10),c("ID","Description" ,"enrichmentScore","p.adjust")]
        res_en$Cell_type<-rep(cname_select,min(dim1,10))
        return (res_en)
      }
      
    }
    
  }
  
}

res_gseareacome<-lapply(unique(res$Cell_type), function(x) {pathway_GSEA_reactome(cname_select=x)})  
res_df_gsereacome <- do.call(rbind,res_gseareacome)


pdf(paste0(outFolder,"gse_reactome_cname_DotPlot.pdf"),width=7,height=10)
ggplot(res_df_gsereacome, 
       aes(x = Cell_type, y = Description)) + 
  geom_point(aes(size = enrichmentScore, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term" #enrichmentScore
  ylab(NULL)+ 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=10)) 
dev.off()






