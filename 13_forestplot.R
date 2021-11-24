library(Matrix)
library(tidyverse)
library(dplyr)
library(org.Hs.eg.db)


##################################################################
### Forest plot on genes that were previously generated for meta plots (cell type specific genes)

# Based on meta plots generated from mashr analysis 
##################################################################


outFolder<-"./13_forestPlots/"
#outFolder<-"./8_outputs_DESeq/"
system(paste0("mkdir -p ",outFolder))


genes<-list.files(paste0("8_outputs_DESeq_batch_library_Plots/meta_plots_selected"))
genes<-unlist(str_split(genes,".pdf"))[seq(1,2*length(genes),by=2)]


cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D","#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","2C3E18","#745B48","#AA5485","#4E747A","#C59A89","#C9C76F")   
names(cluster.Colors)<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte", "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")



# All DE genes from single cell analysis
res <- read_tsv("7_outputs_DESeq_ConditionsByCluster_bath_library/ALL.combined.2021-10-18.tsv")
# Adding location, cell type, and origin columns 
res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
res <- res %>% filter(!is.na(pvalue))
res <- res %>% filter(!is.na(pvalue))

#anno <- read_rds("3_MergeDemux_Output/anno.rds")
eg = bitr(res$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
names(eg)[1]="gene_name"
head(eg)
e2g <- eg$gene_name
names(e2g) <- eg$ENTREZID

#ENTREZID
res <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))


clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
names(clust2Names)<-c(0:23)
res$Cell_type<-clust2Names[res$Cell_type]



# Monocyte<-res %>% filter (Cell_type =="Monocyte"  & padj <0.05 & log2FoldChange>2) %>% select(gene_name)%>% unlist %>% unique
# others<-res %>% filter (Cell_type !="Monocyte"  & padj <0.05 & log2FoldChange < -1) %>% select(gene_name)%>%unlist %>% unique
# 
# onlyMonocyte<-Monocyte[which(!Monocyte %in% others)]
# 
# res_genes<-res[res$gene_name %in% onlyMonocyte,]

res_genes<-res[res$gene_name %in% genes,]

unique(res_genes$gene_name)

#outFolder<-"./13_forestPlots/OnlyMonocyte/"
#outFolder<-"./8_outputs_DESeq/"
system(paste0("mkdir -p ",outFolder))


allgenes<-unique(res_genes$gene_name)

sapply(allgenes, function(gene_show){
  #gene_show<-"ERRFI1"
  res_genes_select<-res_genes %>% filter (gene_name==gene_show)

  # showing all cell types (no transparency set)
  p <- res_genes_select %>%
    ggplot(aes(x=Cell_type, y = log2FoldChange, color=Cell_type)) +
    geom_errorbar(aes(ymax = log2FoldChange + 1.96*lfcSE, ymin = log2FoldChange - 1.96*lfcSE), size=0.8,width=0, position=position_dodge(width=0.7)) +
    geom_point(aes(size=padj,color=Cell_type), position=position_dodge(width=0.7)) +
    ##scale_color_manual(guide = guide_legend(reverse = TRUE) ) +
    scale_size("q-values", trans="log10", range=c(7, 1),limits=c(1E-10,1), breaks=c(1E-12,1E-6,0.001,0.01,0.1)) +
    geom_hline(yintercept=0,lty=2) + 
    scale_color_manual(values=cluster.Colors) +
    xlab("Gene") + ylab(expression(log[2](Fold~Change))) +
    theme_bw() +
    theme(axis.text=element_text(size=20),strip.text.y = element_text(angle = 0,hjust=0,vjust=0.5),strip.background=element_rect(fill="white",color="white")) +
    coord_flip()+
    xlab(NULL)
  ggsave(paste0(outFolder,gene_show,".pdf"),p,width=10,height=8)
  
})



