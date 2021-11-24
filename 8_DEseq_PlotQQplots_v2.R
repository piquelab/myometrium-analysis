library(tidyverse)
library(DESeq2)
library(qqman)
library(org.Hs.eg.db)
library(clusterProfiler)


#outFolder <- paste0("./8_outputs_DESeq_Plots/")
outFolder <- paste0("./8_outputs_DESeq_batch_library_Plots/")
system(paste0("mkdir -p ",outFolder))



########################################################
# load single cell data 

########################################################
res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster/ALL.combined.2021-02-17.tsv")
res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_bath_library/ALL.combined.2021-08-30.tsv")

res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)

## cluster colors 

clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte", "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")

names(clust2Names)<-as.character(c(0:23))
cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D","#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48","#AA5485","#4E747A","#C59A89","#C9C76F")   
names(cluster.Colors)<-clust2Names #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte", "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast", "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")


# cluster labels
res$Cell_type<-clust2Names[res$Cell_type]



#ENTREZID
eg = bitr(res$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
names(eg)[1]="gene_name"
head(eg)
e2g <- eg$gene_name
names(e2g) <- eg$ENTREZID
res <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))



res_basemean_total<-res%>% group_by(Cell_type) %>% summarise(total = n())

res_basemean_basegreater<-res%>% group_by(Cell_type) %>% summarise(basemean_greater = sum(count(baseMean>100)))

res_summary<-res_basemean_total %>% inner_join(res_basemean_basegreater) %>% mutate(percentage=round(basemean_greater/total,3))

write.csv(res_summary,file="8_outputs_DESeq_Plots/res_basemean_percentage.csv")



# Removing na pvalues
# Grouping pvalues based on the Location,Cell_type,and Origin
# Adding a column showing the rank of each pvalue devided by the number of pvalues in each group 

res2 <- res %>% filter(!is.na(pvalue)) %>%
    arrange(pvalue) %>%
    group_by(Cell_type,Origin) %>%
    mutate(r=rank(pvalue, ties.method = "random"),pexp=r/length(pvalue))

####################################################################
# load reference data 
####################################################################
load_ref_data<-function(outFolder="reference_Adi/",fl="CELLECTA.rds")
{
    
    if(fl=="TLvsTNL_blood_ENTREZ" || fl=="TL-TNL_21vs28")
    {
        ref_data <- read.csv(paste0(fl,".csv"),stringsAsFactors = FALSE)
        ref_data<-ref_data %>% dplyr::select(SYMBOL,logFC,P.Value,adj.P.Val,ENTREZ,t)
        colnames(ref_data)<-c("R.gene_name","R.Log2FC","Rpvalue","Rpadj","ENTREZID","Rt")
        ref_data$ENTREZID<-as.character(ref_data$ENTREZID)
        return(ref_data)
        
   
    if (fl=="myometrium_term_TL-TNL_ALLList"|| fl=="PMID31921132")
    {
        ref_data <- read.delim(paste0(fl,".txt"))
        ref_data<-ref_data %>% dplyr::select(SYMBOL,logFC,P.Value,adj.P.Val,ENTREZ,t )
        colnames(ref_data)<-c("R.gene_name","R.Log2FC","Rpvalue","Rpadj","ENTREZID","Rt")
        ref_data <- ref_data %>% filter(!is.na(R.Log2FC) & !is.na(ENTREZID)  & !is.na(Rpadj))
        ref_data$ENTREZID<-as.character(ref_data$ENTREZID)
        return(ref_data)
    }
    else
        if(fl=="myometrium_bulk")
        {
            ref_data <- read_tsv("myo_bulk_TIN_TNL.txt")
            ref_data<-ref_data %>% dplyr::select(SYMBOL,FoldChange,pval.fdr,ENTREZ ,t)
            colnames(ref_data) <- c("R.gene_name","R.Log2FC","Rpadj","ENTREZID","Rt")
            ref_data <- ref_data %>% filter(!is.na(R.Log2FC) & !is.na(ENTREZID)  & !is.na(Rpadj))
            ref_data$ENTREZID<-as.character(ref_data$ENTREZID)
            return(ref_data)
        }
    
    else 
    {
        
        d<-readRDS(paste0("reference_Adi/",fl,".rds"))
        ref_data<-d[["LaborEffect"]]
        if(fl=="PCR.rds")
            ref_data<-ref_data %>% dplyr::select(SYMBOL,logFC,P.Value,adj.P.Val,t)
        else
            ref_data<-ref_data %>% dplyr::select(SYMBOL,log2FoldChange,P.Value,adj.P.Val,t)
        
        colnames(ref_data)<-c("R.gene_name","R.Log2FC","Rpvalue","Rpadj","Rt")
        return (ref_data)
    }
    
    }
}
# reference data 


#experiment<-"myometrium_term_TL-TNL_ALLList"
#experiment<-"TL-TNL_21vs28"
experiment<-"TLvsTNL_blood_ENTREZ"

ref_data<-load_ref_data(fl=experiment) 

ref_data2 <- ref_data %>% filter(!is.na(Rpvalue)) %>%
    arrange(Rpvalue) %>%
    mutate(r=rank(Rpvalue, ties.method = "random"),Rpexp=r/length(Rpvalue))


####################################################################
# simple qq-plot  
####################################################################


fname=paste0(outFolder,"all.qqplot.png");
png(fname,width=800,height=800)
qq(res$pvalue)
dev.off()


####################################################################
# simple qq-plot / cell types with colors 
####################################################################

# qqplot to show the p-values splited by Origin and Location  
fname=paste0(outFolder,"split.qqplot.png");
p1 <- res2 %>%
    ggplot(aes(x=-log10(pexp),y=-log10(pvalue),color=Cell_type)) +
    geom_point() +
    scale_color_manual(values=cluster.Colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Cell Type")) +
    geom_abline(slope=1,intercept=0) +
    #facet_grid(Origin ~ Location) +
    xlab(expression(Expected -log[10](p))) +
    ylab(expression(Observed -log[10](p))) + 
    theme_bw()

ggsave(fname,p1,width=6,height=4.5)



###################################################
# Common between Single cell and bulk
###################################################

# res3<-res2 %>% filter(gene_name %in% ref_data$gene_name)
# # qqplot to show the p-values splited by Origin and Location
# fname=paste0(outFolder,"split.qqplot_overlapgenes_with_refrence.png");
# p1 <- res3 %>%
#     ggplot(aes(x=-log10(pexp),y=-log10(pvalue),color=Cell_type)) +
#     geom_point() +
#     scale_color_manual(values=cluster.Colors) +
#     guides(colour = guide_legend(override.aes = list(size=5),title="Cell Type")) +
#     geom_abline(slope=1,intercept=0) +
#     #facet_grid(Origin ~ Location) +
#     xlab(expression(Expected -log[10](p))) +
#     ylab(expression(Observed -log[10](p))) +
#     theme_bw()
# 
# ggsave(fname,p1,width=6,height=4.5)


res3<-res2 %>% filter(ENTREZID %in% ref_data2$ENTREZID)
# qqplot to show the p-values splited by Origin and Location
fname=paste0(outFolder,"split.qqplot_overlapgenes_with_refrence.png");
p1 <- res3 %>%
    ggplot(aes(x=-log10(pexp),y=-log10(pvalue),color=Cell_type)) +
    geom_point() +
    scale_color_manual(values=cluster.Colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Cell Type")) +
    geom_abline(slope=1,intercept=0) +
    #facet_grid(Origin ~ Location) +
    xlab(expression(Expected -log[10](p))) +
    ylab(expression(Observed -log[10](p))) +
    theme_bw()

ggsave(fname,p1,width=6,height=4.5)


###################################################
# Only Single cell
###################################################

#based on gene symbol
#de_gene_singlecell<-res2 %>% filter (padj<0.1) %>%  select(gene_name)
# res4<-res2 %>% filter(gene_name %in% de_gene_only_ref$gene_name & padj<0.1)
# 
# # qqplot to show the p-values splited by Origin and Location  
# fname=paste0(outFolder,"split.qqplot_overlapgenes_with_refrence.png");
# p1 <- res4 %>%
#     ggplot(aes(x=-log10(pexp),y=-log10(pvalue),color=Cell_type)) +
#     geom_point() +
#     scale_color_manual(values=cluster.Colors) +
#     guides(colour = guide_legend(override.aes = list(size=5),title="Cell Type")) +
#     geom_abline(slope=1,intercept=0) +
#     #facet_grid(Origin ~ Location) +
#     xlab(expression(Expected -log[10](p))) +
#     ylab(expression(Observed -log[10](p))) + 
#     theme_bw()
# 
# ggsave(fname,p1,width=6,height=4.5)


de_gene_singlecell<-res2 %>% filter (padj<0.1) %>%  dplyr::select(ENTREZID)
de_gene_only_ref<-ref_data2 %>% filter (Rpadj < 0.1 & !ENTREZID %in% de_gene_singlecell$ENTREZID )
res4<-res2 %>% filter(!ENTREZID %in% de_gene_only_ref$ENTREZID)

# qqplot to show the p-values splited by Origin and Location  
fname=paste0(outFolder,"split.qqplot_nonoverlapgenes_with_refrence.png");
p1 <- res4 %>%
    ggplot(aes(x=-log10(pexp),y=-log10(pvalue),color=Cell_type)) +
    geom_point() +
    scale_color_manual(values=cluster.Colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Cell Type")) +
    geom_abline(slope=1,intercept=0) +
    #facet_grid(Origin ~ Location) +
    xlab(expression(Expected -log[10](p))) +
    ylab(expression(Observed -log[10](p))) + 
    theme_bw()
ggsave(fname,p1,width=6,height=4.5)



###################################################
# Only Bulk
###################################################


system(paste0("mkdir -p ",outFolder,experiment,"/"))
# single qqplot for ref  #qq(ref_data2$Rpvalue)
fname1=paste0(outFolder,experiment,"/reference.qqplot.png");
p1 <- ref_data2 %>%
    ggplot(aes(x=-log10(Rpexp),y=-log10(Rpvalue))) +
    geom_point() +
    geom_abline(slope=1,intercept=0) +
    #facet_grid(Origin ~ Location) +
    xlab(expression(Expected -log[10](p))) +
    ylab(expression(Observed -log[10](p))) + 
    theme_bw()
ggsave(fname1,p1,width=6,height=4.5)



# multiple qqplots related to show enrichment at different cell types
# preparing data to show the enrichment across multiple cell types


colnames(ref_data)[which(colnames(ref_data)=="ENTREZID")]<-"R.ENTREZID"
colnames(ref_data)[which(colnames(ref_data)=="R.gene_name")]<-"gene_name"
resDE<-res2 %>% filter(padj<0.1 ) #single cell fdr 0.1
res_join<-ref_data %>% inner_join(resDE)
res_join <- res_join %>% filter(!is.na(Rpvalue)) %>%
    arrange(Rpvalue) %>%
    group_by(Cell_type,Origin) %>%
    mutate(r=rank(Rpvalue, ties.method = "random"),Rpexp=r/length(Rpvalue))
fname2=paste0(outFolder,experiment,"/reference.enriched.celltypes.qqplot.png");
p2 <- res_join %>%
    ggplot(aes(x=-log10(Rpexp),y=-log10(Rpvalue),color=Cell_type)) +
    geom_point() +
    scale_color_manual(values=cluster.Colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Cell Type")) +
    geom_abline(slope=1,intercept=0) +
    #facet_grid(Origin ~ Location) +
    xlab(expression(Expected -log[10](p))) +
    ylab(expression(Observed -log[10](p))) +
    theme_bw()
ggsave(fname2,p2,width=6,height=4.5)



#based on gene entrezid

de_gene_bulk<-ref_data %>% filter (Rpadj<0.1) %>%  dplyr::select(ENTREZID)
de_gene_only_singlecell<-res2 %>% filter (padj < 0.1 & !ENTREZID %in% de_gene_bulk$ENTREZID )
ref_data3<-ref_data %>% filter(!ENTREZID %in% de_gene_only_singlecell$ENTREZID)
res5<-res2 %>% filter(!ENTREZID %in% ref_data3$ENTREZID)


# qqplot to show the p-values splited by Origin and Location  
fname=paste0(outFolder,"split.qqplot_bulk_only_nonoverlapgenes_with_singlecell.png");
p1 <- res5 %>%
    ggplot(aes(x=-log10(pexp),y=-log10(pvalue),color=Cell_type)) +
    geom_point() +
    scale_color_manual(values=cluster.Colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Cell Type")) +
    geom_abline(slope=1,intercept=0) +
    #facet_grid(Origin ~ Location) +
    xlab(expression(Expected -log[10](p))) +
    ylab(expression(Observed -log[10](p))) + 
    theme_bw()
ggsave(fname,p1,width=6,height=4.5)




################################
## 
# Tarca, Adi L., Roberto Romero, Zhonghui Xu, Nardhy Gomez-Lopez, Offer Erez, Chaur-Dong Hsu, Sonia S. Hassan, and Vincent J. Carey. "Targeted expression profiling by RNA-Seq improves detection of cellular dynamics during pregnancy and identifies a role for T cells in term parturition." Scientific reports 9, no. 1 (2019): 1-13.
################################


ref_data<-load_adi(fl="PCR.rds")
outFolder<-"./8_outputs_DESeq_Plots/PCR/"
de_gene_singlecell<-res2 %>% filter (padj<0.1) %>%  dplyr::select(gene_name)
de_gene_only_ref<-ref_data %>% filter (Rpadj < 0.1 & !gene_name %in% de_gene_singlecell$gene_name )
res3<-res2 %>% filter(gene_name %in% de_gene_only_ref$gene_name)


ref_data<-load_adi(fl="CELLECTA.rds")
outFolder<-"./8_outputs_DESeq_Plots/DriverMap/"
de_gene_singlecell<-res2 %>% filter (padj<0.1) %>%  dplyr::select(gene_name)
de_gene_only_ref<-ref_data %>% filter (Rpadj < 0.1 & !gene_name %in% de_gene_singlecell$gene_name )
res3<-res2 %>% filter(gene_name %in% de_gene_only_ref$gene_name)



ref_data<-load_adi(fl="RNASeq.rds")
outFolder<-"./8_outputs_DESeq_Plots/RNASeq/"
de_gene_singlecell<-res2 %>% filter (padj<0.1) %>%  dplyr::select(gene_name)
de_gene_only_ref<-ref_data %>% filter (Rpadj < 0.1 & !gene_name %in% de_gene_singlecell$gene_name )
res3<-res2 %>% filter(gene_name %in% de_gene_only_ref$gene_name)

fname=paste0(outFolder,"split.qqplot_DE_only_refrence.png");
p1 <- res3 %>%
    ggplot(aes(x=-log10(pexp),y=-log10(pvalue),color=Cell_type)) +
    geom_point() +
    scale_color_manual(values=cluster.Colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Cell Type")) +
    geom_abline(slope=1,intercept=0) +
    #facet_grid(Origin ~ Location) +
    xlab(expression(Expected -log[10](p))) +
    ylab(expression(Observed -log[10](p))) + 
    theme_bw()
ggsave(fname,p1,width=6,height=4.5)







