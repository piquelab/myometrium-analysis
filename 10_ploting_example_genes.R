
##################################################################
### Plotting the top 3 down-regulated and up-regulated genes ###
##  Plotting DE genes per cell types
##################################################################

library(Seurat)
library(Matrix)
library(tidyverse)
library(tidyverse)
library(dplyr)
library(stringr)


outFolder <- paste0("./8_example_genes_Plots/")
system(paste0("mkdir -p ",outFolder))


##loading all DE genes FDR < 0.1
res <- read_tsv("7_outputs_DESeq_ConditionsByCluster/SIG.combined.2021-02-17.tsv")


# Adding location, cell type, and origin columns 
res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
res <- res %>% filter(!is.na(pvalue))
res<-res[order(abs(res$log2FoldChange), abs(res$baseMean),decreasing = TRUE),]


clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
names(clust2Names)<-c(0:23)
res$Cell_type<-clust2Names[res$Cell_type]


################################################
# Violin pot
################################################

selected_cell_types<-c("Stromal-1","Endothelial-1","Monocyte","Stromal-2","Macrophage-1","Macrophage-2","Decidual","LED","Smooth muscle cells-1","Macrophage-3")
res_examplegenes<-res  %>% filter (Cell_type %in% selected_cell_types)

res_examplegenes<-res_examplegenes %>% group_by(Cell_type) %>% arrange(log2FoldChange) 

top5 <- res_examplegenes %>% group_by(Cell_type) %>% top_n(n = 5, wt = log2FoldChange)
bottom5<-res_examplegenes %>% group_by(Cell_type) %>% top_n(n = 5, wt = -1*log2FoldChange)
res_violin<-rbind(top5,bottom5)


## Load the data seruat object and the list of DE genes.
sc <- read_rds("6_harmony_cellClass_plots_res0.8_final/SeuratObject.rds")

#meta data
md <- sc@meta.data



a2 <- FetchData(sc,c("UMAP_1","UMAP_2","cluster_name","Location","Condition","Origin","Pregnancy_ID","Condition")) 

myscale = 1/colSums(sc@assays$RNA@counts)*1000000


sc$cluster_name<-clust2Names[sc$seurat_clusters]
##rec <- map_dfr(1:nrow(res),function(ii){
rec <- map_dfr(1:(dim(res_violin)[1]),function(ii){
  a2 <- FetchData(sc,c("UMAP_1","UMAP_2","Location","cluster_name","Condition","Origin","Pregnancy_ID")) 
  a2$gene_name=res_violin$gene_name[ii]
  a2$cname=paste(res_violin$gene_name[ii],res_violin$cname[ii],sep="_")
  a2$Expression=sc@assays$RNA@counts[res_violin$kbid[ii],]*myscale
  a2$log2FoldChange<-res_violin$log2FoldChange[ii]
  a2$baseMean<-res_violin$baseMean[ii]
  aux <- a2 %>% filter(Origin==res_violin$Origin[ii],cluster_name==res_violin$Cell_type[ii])
  ## a2$Expression=sc2@assays$RNA@scale.data[selGenes$kbid[ii],]
  aux
})
dim(rec)
# splited violin pot
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                          draw_group = function(self, data, ..., draw_quantiles = NULL) {
                            data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                            grp <- data[1, "group"]
                            newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                            newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                            newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                            
                            if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                              stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                        1))
                              quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                              aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                              aesthetics$alpha <- rep(1, nrow(quantiles))
                              both <- cbind(quantiles, aesthetics)
                              quantile_grob <- GeomPath$draw_panel(both, ...)
                              ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                            }
                            else {
                              ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                            }
                          })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}



rec<-rec %>% group_by(cluster_name) %>% arrange(log2FoldChange,baseMean, .by_group = TRUE) 
rec$cluster_name <- factor(rec$cluster_name,levels=unique(rec$cluster_name))


#######################################################
#violin plot split
#######################################################
pdf(paste0(outFolder,"violinplot_genes_per_celltypes_immune_cells_v2.pdf"),width=60,height=10)
#pdf(paste0(outFolder,"violinplot_genes_per_celltypes_non-immune_cells_v2.pdf"),width=30,height=10)

#gg <-

rec %>% ggplot(aes(x=reorder(gene_name,(log2FoldChange)),y=log10(Expression),fill=Condition)) + 
  #geom_violin()+
  #c(0.25, 0.5, 0.75)
  #geom_violin(draw_quantiles = c(0.5),col="white",position = position_dodge(width = 0.9)) +
  geom_split_violin()+
  #geom_boxplot(width=.1)+
  #scale_color_manual(values=c("TNL"="#333399","TIL"="#A50021"))
  #scale_fill_manual(values=c("Control"="#3333A3","w/COVID-19"="tomato"))+
  #scale_fill_manual(values=c("Control"="#3333A3","w/COVID-19"="#A50021"))+
  scale_fill_manual(values=c("TNL"="#333399","TIL"="#A50021"))+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.title.x=element_blank())+
  theme_bw()+
  #geom_jitter()+
  #geom_jitter(shape=16, position=position_jitter(width = NULL, height = NULL))+
  theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +  #,text = element_text(size=30)
  facet_grid(~cluster_name, scale="free",space="free")+
 scale_y_continuous(limits = c(1, 5))
dev.off()




################################################
#Forest Plot 
################################################

res_select<-res %>% group_by(Cell_type) %>% filter (padj <=0.01 & abs(log2FoldChange)>=2)
DE_number_per_celltype<-tapply(res_select$padj ,factor(res_select$Cell_type),length )

res_select$DE_number_per_celltype<-DE_number_per_celltype[res_select$Cell_type]

res_select<-res_select[order(res_select$DE_number_per_celltype,decreasing = TRUE),]


res_select$Cell_type <- factor(res_select$Cell_type,levels=unique(res_select$Cell_type))


## forestPlot 
#coloring
res_select$padj[res_select$padj<1E-6]<- 1E-6

p<-res_select %>% 
  ggplot(aes(x=reorder(gene_name,(log2FoldChange),color="black"),y=log2FoldChange),color="black") +
  geom_errorbar(aes(ymax = log2FoldChange + 1.96*lfcSE, ymin = log2FoldChange - 1.96*lfcSE),size=1.2,color="black",alpha=1,width=0,position=position_dodge(width=0.2)) +
  geom_point(aes(size=padj),position=position_dodge(width=0.2),color="#111111") +
      ##scale_color_manual(guide = guide_legend(reverse = TRUE) ) +
    scale_size("q-values", trans="log10", range=c(7, 1),limits=c(1E-10,1), breaks=c(1E-12,1E-6,0.001,0.01,0.1)) +
    ## scale_alpha_manual(values=c(0.3, 1.0),guide=FALSE) +
  geom_hline(yintercept=0,lty=2) + 
  xlab("Gene") + ylab(expression(log[2](Fold~Change))) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0,hjust=0,vjust=0.5),strip.background=element_rect(fill="white",color="white")) +
  coord_flip() + 
  #theme(strip.background =element_rect(fill=cluster.Colors))+
  facet_grid(Cell_type ~ . ,scales = "free_y",space="free") 

g <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-', g$layout$name))
fills <- c("Stromal-1"="#DF7D99"  ,"Monocyte"="#2BA3D3","Endothelial-1"="#FFC000","Stromal-2"="#ED4315","Macrophage-3"="#72C6C8"  ,"Macrophage-1"="#4E65A6",
           "Decidual"="#D14357", "Macrophage-2"="#838EDF",
           "Endothelial-2"="#E4652C" ,"LED"="#D5438E" )
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col<-fills[k]
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$lwd=4
  k <- k+1
}
ggsave(paste0(outFolder,"forestPlot_DEgenes_cname_v5.pdf"),g,width=10,height=30)
