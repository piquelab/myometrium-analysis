library(tidyverse)
##library(knitr)
library(DESeq2)
##library(annotables)
library(qqman)
library(org.Hs.eg.db)
library(clusterProfiler)





outFolder<-paste0("./14_deconvolution_analysis/merge_subcelltypes/CIBERSORTx_Job17_limma_DEGs/")

#outFolder<-paste0("./14_deconvolution_analysis/merge_subcelltypes/CIBERSORTx_Job17_deseq_DEGs/")

system(paste0("mkdir -p ",outFolder))




########################################################
# load single cell data 
########################################################
res <- read_tsv(paste0(outFolder,"ALL.combined.tsv"))
colnames(res)[which(colnames(res)=="celltype")]<-"Cell_type"



## cluster colors 

res$Cell_type[res$Cell_type=="Smoothmusclecells"]<-"Smooth Muscle Cells"
res$Cell_type[res$Cell_type=="UnciliatedEpithelial"]<-"Unciliated Epithelial"

cluster.Colors<-c("#DF7D99","#9B91B9","#C59A89" ,"#D14357" , "#FFC000","2C3E18" ,"#D5438E" ,"#2BA3D3" , "#F88091","#7BC791" ,"#CA8588" ,"#C9C76F"   )

names(cluster.Colors)<-c("Stromal","DC","B-cell","Decidual","Endothelial" ,"EVT" ,"LED","Monocyte","Myofibroblast","NK-cell","Smooth Muscle Cells","Unciliated Epithelial" )   



clust2Name<-c("Stromal","Macrophage","Macrophage","Endothelial","Monocyte","T-cell","Decidual","T-cell","LED","Stromal","ILC","NK-cell","Smooth Muscle Cells","Myofibroblast","Macrophage","Endothelial","DC","Smooth Muscle Cells","EVT","Plasmablast","Smooth Muscle Cells","Macrophage","B-cell","Unciliated Epithelial")
names(clust2Name)<-c(0:23)

# Removing na pvalues
# Grouping pvalues based on the Location,Cell_type,and Origin
# Adding a column showing the rank of each pvalue devided by the number of pvalues in each group 

res2 <- res %>% filter(!is.na(pvalue)) %>%
    arrange(pvalue) %>%
    group_by(Cell_type) %>%
    mutate(r=rank(pvalue, ties.method = "random"),pexp=r/length(pvalue))

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



