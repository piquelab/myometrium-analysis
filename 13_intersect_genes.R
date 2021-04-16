setwd("/wsu/home/groups/prbgenomics/labor_myo/myometrium_analysis/")
library(tidyverse)
#library(dplyr)

## load data 
load_ref_data<-function(outFolder="reference_Adi/",fl="CELLECTA.rds")
{
  
  if(fl=="supp_adi.txt") 
  {
    ref_data<-read.delim("supp_adi.txt")
    ref_data<-ref_data %>%filter(adj.P.Val<0.1 & !is.na(adj.P.Val) &!is.na(FC))
    ref_data<-ref_data %>%select("SYMBOL", "FC" ,"P.Value" ,"adj.P.Val")
    colnames(ref_data)<-c("gene_name", "R.Log2FC" ,"Rpvalue" ,"Rpadj")
    return (ref_data)
  }
  else 
    if(fl=="myometrium_bulk")
  {
    ref_data <- read_tsv("myo_bulk_TIN_TNL.txt")
    ref_data<-ref_data[,c(-1,-3)]
    colnames(ref_data) <- c("gene_name","R.Log2FC","Rpadj")
    ref_data <- ref_data %>% filter(!is.na(R.Log2FC))
    return(ref_data)
  }
  
  else 
    if (fl=="PMID31921132")
  {
    ref_data<-read.delim("PMID31921132.txt")
    ref_data<-ref_data %>% select(SYMBOL,FC,P.Value,adj.P.Val )
    colnames(ref_data)<-c("gene_name","R.Log2FC","Rpvalue","Rpadj")
    return (ref_data)
  }
  
  else 
  {
    
    d<-readRDS(paste0(outFolder,fl))
    ref_data<-d[["LaborEffect"]]
    if(fl=="PCR.rds")
      ref_data<-ref_data %>% select(SYMBOL,logFC,P.Value,adj.P.Val)
    else
      ref_data<-ref_data %>% select(SYMBOL,log2FoldChange,P.Value,adj.P.Val)
    
    colnames(ref_data)<-c("gene_name","R.Log2FC","Rpvalue","Rpadj")
    return (ref_data)
  }
  
}


# our data
clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(c(0:23),"_",clust2Names)
names(clust2Names)<-c(0:23)
res <- read_tsv("7_outputs_DESeq_ConditionsByCluster/SIG.combined.2021-02-17.tsv")
res<-res %>%filter(padj<0.1)
res<-res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
res <- res %>% filter(!is.na(pvalue))
res$Cell_type<-clust2Names[res$Cell_type]


################################################################################
#Targeted expression profiling by RNA-Seq improves detection of cellular dynamics during pregnancy and identifies a role for T cells in term parturition

ref_data<-load_ref_data(fl="RNASeq.rds")
outFolder<-"./8_outputs_DESeq_Plots/RNASeq/"
ref_data<-ref_data %>% filter(Rpadj<0.1 )
#intersected_genes<-read.csv(paste0(outFolder,"intersected_genes.csv"),stringsAsFactors = FALSE)
resJoin_RNASeq <- res %>% inner_join(ref_data) 
resJoin_RNASeq<-resJoin_RNASeq %>% select("gene_name","Cell_type","pvalue","padj","log2FoldChange","Rpadj")
colnames(resJoin_RNASeq)<-c("gene_name","Cell_type","pvalue","padj","log2FoldChange","padj_RNASeq")


outFolder<-"./8_outputs_DESeq_Plots/DriverMap/"
ref_data<-load_ref_data(fl="CELLECTA.rds")
ref_data<-ref_data %>% filter(Rpadj<0.1 )
#intersected_genes<-read.csv(paste0(outFolder,"intersected_genes.csv"),stringsAsFactors = FALSE)
ref_data <- ref_data %>% filter(!is.na(Rpadj))
resJoin_drivermap <- res %>% inner_join(ref_data) 
resJoin_drivermap<-resJoin_drivermap %>% select("gene_name","Cell_type","pvalue","padj","log2FoldChange","Rpadj")
colnames(resJoin_drivermap)<-c("gene_name","Cell_type","pvalue","padj","log2FoldChange","padj_drivermap")




################################################################################
#The Cellular Transcriptome in the Maternal Circulation During Normal Pregnancy: A Longitudinal Study
ref_data<-load_ref_data(fl="supp_adi.txt")
ref_data <- ref_data %>% filter(!is.na(Rpadj))
resJoin_GA <- res %>% inner_join(ref_data) 
resJoin_GA<-resJoin_GA %>% select("gene_name","Cell_type","pvalue","padj","log2FoldChange","Rpadj")
colnames(resJoin_GA)<-c("gene_name","Cell_type","pvalue","padj","log2FoldChange","padj_GA")


combined_intersect_matrix<-resJoin_GA %>% dplyr::full_join (resJoin_drivermap)%>% full_join (resJoin_RNASeq)
write.csv(combined_intersect_matrix,file="combined_intersect_matrix.csv")






genes<-unique(c(resJoin_GA$gene_name, resJoin_drivermap$gene_name,resJoin_RNASeq$gene_name))
matrix(nrow=length(genes), ncol=)

system(paste0("mkdir -p ", paste0(outFolder,"intersect/")))
write.csv(res_intersect,file=paste0(outFolder,"intersect/res_intersect_ref_Adi.csv"))
write.csv(gene_intersect,file=paste0(outFolder,"intersect/gene_intersect_intersect_ref_Adi.csv"))
