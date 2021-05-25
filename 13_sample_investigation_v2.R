library(Seurat)
library(Matrix)
library(tidyverse)

library(future)

library(harmony)


library(enrichplot)
library(dplyr)
library(org.Hs.eg.db)
library(ReactomePA)
#################
##library(SingleR)




future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)



outFolder="13_sample_investigation_plots_V2/"
system(paste0("mkdir -p ", outFolder))


###########################################
clust2Name<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
              "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast",
              "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")

#clust2Name<-paste0(c(0:23),"_",clust2Name)
names(clust2Name)<-c(0:23)

sc <- read_rds("4_harmony/sc.NormByLibrary.Harmony.StringentFiltering.res0.8.rds")


sc@meta.data$cluster_name <- clust2Name[sc@meta.data$seurat_clusters]


aa <- FetchData(sc,c("UMAP_1","UMAP_2","seurat_clusters","Location","Condition","Origin","status","SNG.BEST.GUESS","Pregnancy_ID","cluster_name")) 



############### Forest plot for specific genes 

# so same code as in the forest plot, but you only show those genes, and for all cell-types available, 
# and you arrange/sort by gene name and color by cell-type, this way we can see something like the meta plot but using the forest plot.
# basically you can use the same code as in the forest plot, but you color by cell-type, and you can sort/arrange/facet by gene name.


cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
                  "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","2C3E18","#745B48",
                  "#AA5485","#4E747A","#C59A89","#C9C76F")   
names(cluster.Colors)<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
                         "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast",
                         "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")



##outFolder="13_sample_investigation_plots/"

res <- read_tsv("7_outputs_DESeq_ConditionsByCluster/ALL.combined.2021-02-17.tsv")


# Adding location, cell type, and origin columns 
res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
res <- res %>% filter(!is.na(pvalue))

clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
names(clust2Names)<-c(0:23)
res$Cell_type<-clust2Names[res$Cell_type]

#res_genes<-res[res$gene_name %in% c("CALM1","CALM2","PRKG1", "PDE5A",  "PLN", "GUCY1A1", "GUCY1B1", "LMOD1", "JUN","IGFBP5"),]
res_genes<-res[res$gene_name %in% c("CALM1","PPP1R16A","PRKG1", "PDE5A",  "PLN", "GUCY1A1", "GUCY1B1", "LMOD1", "JUN","IGFBP5"),]

unique(res_genes$gene_name)

##,size=(Cell_type=="Smooth muscle cells-1")

p <- res_genes %>%
  ggplot(aes(x=gene_name, y = log2FoldChange, color=Cell_type, alpha = (Cell_type=="Smooth muscle cells-1"))) +
  geom_errorbar(aes(ymax = log2FoldChange + 1.96*lfcSE, ymin = log2FoldChange - 1.96*lfcSE), size=1.0, width=0, position=position_dodge(width=0.7)) +
  geom_point(aes(size=padj,color=Cell_type), position=position_dodge(width=0.7)) +
  ##scale_color_manual(guide = guide_legend(reverse = TRUE) ) +
  scale_size("q-values", trans="log10", range=c(7, 1),limits=c(1E-10,1), breaks=c(1E-12,1E-6,0.001,0.01,0.1)) +
  geom_hline(yintercept=0,lty=2) + 
  scale_color_manual(values=cluster.Colors) +
  xlab("Gene") + ylab(expression(log[2](Fold~Change))) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0,hjust=0,vjust=0.5),strip.background=element_rect(fill="white",color="white")) +
  coord_flip()+
  xlab(NULL)
#facet_grid(gene_name ~ Cell_type) 
##facet_grid(~gene_name ) 
#facet_grid(gene_name ~  ,scales = "free_y",space="free") 

ggsave(paste0(outFolder,"forestPlot_metaplot_v2.pdf"),p,width=8,height=12)


