library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Matrix)
library(tidyverse)
library(Seurat)


##########################################################################
# Whether the cell-cell communication is enhanced or not
# The interaction between which cell types is significantly changed
# How the major sources and targets change from one condition to another
##########################################################################





clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
names(clust2Names)<-c(0:23)
# cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D","#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48","#AA5485","#4E747A","#C59A89","#C9C76F")   
# names(cluster.Colors)<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte", "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast", "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")


cell.type.annotation<-read.delim("data.colors.mouse.human.txt")
cluster.Colors<-cell.type.annotation$cluster.Colors_human2
names(cluster.Colors)<-cell.type.annotation$cluster_human
cluster.Colors<-cluster.Colors[cluster.Colors!=""]



# res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
# res <- res %>% filter(!is.na(pvalue))
# res<-res %>% filter(padj<0.1)
# res$Cell_type<-clust2Names[res$Cell_type]


threshold<-50
#threshold<-100
#threshold<-200
outFolder<-paste0("./10_CellChat_comparison_conditions_after_filter_",threshold,"_plots/")
#outFolder="./10_CellChat_comparison_conditions_after_filter_200_plots/"

system(paste0("mkdir -p ", outFolder))


future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)


####################################
# Load data
####################################


conditions<-c("TNL" ,"TIL")


##################################################################
# circle plots for paper 
##################################################################





##################################################
#comparison between ecoli and pbs plots
######################################################

#threshold<-100
#threshold<-200
outFolder<-paste0("./10_CellChat_comparison_conditions_after_filter_",threshold,"_plots/")

  cellchat_TNL<-cellchat_TNL_original<-read_rds(paste0("./10_CellChat_analysis_default_after_filter_",threshold,"/","cellchat_TNL","_","2022-02-11",".rds"))
  cellchat_TIL<-cellchat_TIL_original<-read_rds(paste0("./10_CellChat_analysis_default_after_filter_",threshold,"/","cellchat_TIL","_","2022-02-11",".rds"))
 
  
  # pathways_TNL<-cellchat_TNL_original@netP$pathways
  # controls<-cbind(pathways_TNL, rep("control",length(pathways_TNL)))
  # controls<-as.data.frame(controls)
  # colnames(controls)<-c("pathway","control")
  # 
  # 
  # pathways_TIL<-cellchat_TIL_original@netP$pathways
  # Ecolis<-cbind(pathways_TIL, rep("TIL",length(pathways_TIL)))
  # Ecolis<-as.data.frame(Ecolis)
  # colnames(Ecolis)<-c("pathway","TIL")
  # 
  # allpathways<-Ecolis %>% full_join(controls)
  # 
  # allpathways<-allpathways %>% arrange(desc(!is.na(control),(!is.na(TIL))))
  # 
  #write.csv(allpathways,file=paste0(outFolder,"pathways.csv"))
  
  # shared_pathways<-intersect(pathways_TNL,pathways_TIL)
  # pathways_TNL_only<-pathways_TNL[!pathways_TNL %in% pathways_TIL]
  # pathways_TIL_only<-pathways_TIL[! pathways_TIL %in% pathways_TNL]
  # 
  
  
  
  # there is an additional population  30_B cell specific to cellchat_TNL compared to Ecoli
  # we lift up Ecoli by lifting up the cell groups to the same cell labels as cellchat_TNL 
  
  TNL_cells<-rownames(cellchat_TNL@net$prob)
  TIL_cells<-rownames(cellchat_TIL@net$prob)
  
  TIL_cells_only<-TIL_cells [!TIL_cells %in%TNL_cells ]
  TNL_cells_only<-TNL_cells [! TNL_cells  %in%TIL_cells ]
  
  if (length(TIL_cells_only)==0)
  {
    group.new = levels(cellchat_TNL@idents)
    cellchat_TIL <- liftCellChat(cellchat_TIL, group.new)
  }else if (length(TNL_cells_only)==0){
    group.new = levels(cellchat_TIL@idents)
    cellchat_TNL <- liftCellChat(cellchat_TNL, group.new)} else 
  {
    group.new = union(levels( cellchat_TIL@idents), levels(cellchat_TNL@idents))
    cellchat_TNL <- liftCellChat(cellchat_TNL, group.new)
    cellchat_TIL <- liftCellChat(cellchat_TIL, group.new)
    }
  
  
  # now merge
  object.list <- list(TNL = cellchat_TNL, TIL = cellchat_TIL)
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))
  
  color.use.TIL=cluster.Colors[rownames(cellchat@netP$TIL$prob)]
  color.use.TNL=cluster.Colors[rownames(cellchat@netP$TNL$prob)]
  coloruses<-list(color.use.TNL,color.use.TIL)
  
  overalcoloruses<-cluster.Colors[unique(c(rownames(cellchat@netP$TNL$prob),rownames(cellchat@netP$TIL$prob)))]
  
  
#   
  pdf(paste0(outFolder,"informationflow.pdf"),width=10,height=16)
  #gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
  gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use=c("#333399","#A50021"),font.size=12)
  gg2
  dev.off()
#   
#   
#   

  
  pdf(paste0(outFolder,"diffInteraction.pdf"),width=24,height=24)
  #par(mfrow = c(1,2), xpd=TRUE)
  #gg1<-netVisual_diffInteraction(cellchat, weight.scale = T,top=0.25,color.use=cluster.Colors[rownames(cellchat@netP$Ecoli$prob)])
  gg2<-netVisual_diffInteraction(cellchat, arrow.width = 2,arrow.size=1, vertex.weight = 15, vertex.size.max=8,weight.scale = T, measure = "weight",top=0.25,color.use=cluster.Colors[rownames(cellchat@netP$TIL$prob)],vertex.label.cex = 1,edge.width.max=20)
  #gg1 + gg2
  gg2
  dev.off()
  
  
  pdf(paste0(outFolder,"diffInteraction_nolabel.pdf"),width=24,height=24)
  #par(mfrow = c(1,2), xpd=TRUE)
  #gg1<-netVisual_diffInteraction(cellchat, weight.scale = T,top=0.25,color.use=cluster.Colors[rownames(cellchat@netP$Ecoli$prob)])
  gg2<-netVisual_diffInteraction(cellchat, arrow.width = 2,arrow.size=1, vertex.weight = 15, vertex.size.max=8,weight.scale = T, measure = "weight",top=0.25,color.use=cluster.Colors[rownames(cellchat@netP$TIL$prob)],vertex.label.cex = 0.00001,edge.width.max=20)
  #gg1 + gg2
  gg2
  dev.off()
  
  
  
  weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
  pdf(paste0(outFolder,"interactions_strength_conditions_nolabel.pdf"),width=26,height=15)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_circle(object.list[[i]]@net$weight, arrow.width = 2,arrow.size=1, vertex.weight = 15, vertex.size.max=8, vertex.label.cex=0.000001,top=0.25,weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction strength - ", names(object.list)[i]),color.use = coloruses[[i]])
  }
  dev.off()
  
  
  weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
  pdf(paste0(outFolder,"interactions_strength_conditions.pdf"),width=26,height=15)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_circle(object.list[[i]]@net$weight,arrow.width = 2,arrow.size=1, vertex.weight = 15, vertex.size.max=8,top=0.25,weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction strength - ", names(object.list)[i]),color.use = coloruses[[i]])
  }
  dev.off()
  
#   
  pdf(paste0(outFolder,"diffInteraction_heatmap.pdf"),width=6,height=6)
  #gg1 <- netVisual_heatmap(cellchat,color.use=cluster.Colors[rownames(cellchat@netP$Ecoli$prob)])
  #> Do heatmap based on a merged object
  gg2 <- netVisual_heatmap(cellchat, measure = "weight",color.use=cluster.Colors[rownames(cellchat@netP$TIL$prob)])
  #> Do heatmap based on a merged object
  #gg1 + gg2
  gg2
  dev.off()

  color.use.TIL=cluster.Colors[rownames(cellchat@netP$TIL$prob)]
  color.use.TNL=cluster.Colors[rownames(cellchat@netP$TNL$prob)]
  coloruses<-list(color.use.TNL,color.use.TIL)


  
  
#############################################################################
# scatter plot  (arrow plot) 
#############################################################################

  clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
  names(clust2Names)<-c(0:23)
  cell.type.annotation<-read.delim("data.colors.mouse.human.txt")
  cluster.Colors<-cell.type.annotation$cluster.Colors_human2
  names(cluster.Colors)<-cell.type.annotation$cluster_human
  cluster.Colors<-cluster.Colors[cluster.Colors!=""]
  

  
pdf(paste0(outFolder,"outgoing_incoming_conditions.pdf"),width=19,height=10)
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  object.listi <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP")

  gg[[i]] <- netAnalysis_signalingRole_scatter(object.listi, title = names(object.list)[i], weight.MinMax = weight.MinMax,color.use = coloruses[[i]])
}
 p<-patchwork::wrap_plots(plots = gg)
 plot(p)
 dev.off()

 
label.size = 4
dot.size = c(2, 6)
dot.alpha = 0.6
font.size.title = 13
font.size = 15
df1<-gg[[1]]$data
df1$group<-names(object.list)[1]
df2<-gg[[2]]$data
df2$group<-names(object.list)[2]
newdf<-rbind(df1,df2)
newdf$color_custome<-cluster.Colors[newdf$labels]

xlabel = "Outgoing interaction strength"
ylabel = "Incoming interaction strength"

overalcoloruses<-overalcoloruses[newdf$labels]


#############################################################################
# scatter plot with lables 
#############################################################################

require(grid)



newdf<-newdf %>% arrange(labels)
# newdf$color_custome <- factor(newdf$color_custome,levels=unique(newdf$color_custome))
newdf$labels <- factor(newdf$labels,levels=unique(newdf$labels))
# newdf$group <- factor(newdf$group,levels=unique(newdf$group))


gg <- ggplot(data  = newdf, aes(x, y),show.legend = F) + geom_point(aes(size = Count, colour = labels, fill = labels),show.legend = F)
gg <- gg + CellChat_theme_opts() + theme(text = element_text(size = font.size), legend.key.height = grid::unit(0.15, "in")) + labs(title = "TIL + TNL", x = xlabel, y = ylabel) + theme(plot.title = element_text(size = font.size.title, face = "plain")) + theme(axis.line.x = element_line(size = 0.25),                                                                                                                                                                               axis.line.y = element_line(size = 0.25))
#gg <- gg + scale_fill_manual(values = ggplot2::alpha(overalcoloruses, alpha = dot.alpha), drop = FALSE) + guides(fill = FALSE)
gg <- gg + scale_size_continuous(limits = weight.MinMax, range = dot.size)
gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels, colour = labels), size = label.size, show.legend = F, segment.size = 0.2, segment.alpha = 0.5)

gg<-gg+geom_path(aes(colour=labels,size = 1.5),arrow = arrow(type = "closed",length=unit(0.2, "inches")))
#gg <- gg + scale_fill_manual(values = ggplot2::alpha(overalcoloruses, alpha = dot.alpha), drop = FALSE) + guides(fill = FALSE)
gg <- gg + scale_colour_manual(values =overalcoloruses, drop = FALSE) + guides(colour = FALSE)+ guides(fill = FALSE)

pdf(paste0(outFolder,"outgoing_incoming_conditions_both.pdf"),width=12,height=8)
gg+ theme(legend.position = "none")
dev.off()


#############################################################################
# scatter plot without lables 
#############################################################################

require(grid)
#newdf <- data.frame(x = outgoing.cells, y = incoming.cells, labels = names(incoming.cells), Count = num.link)
#newdf$labels <- factor(newdf$labels, levels = names(incoming.cells))
gg <- ggplot(data = newdf, aes(x, y)) + geom_point(aes(size = Count, colour = labels, fill = labels),show.legend = FALSE)
gg <- gg + CellChat_theme_opts() + theme(text = element_text(size = font.size), legend.key.height = grid::unit(0.15, "in")) + labs(title = "Ecoli + control", x = xlabel, y = ylabel) + theme(plot.title = element_text(size = font.size.title, face = "plain")) + theme(axis.line.x = element_line(size = 0.25),axis.line.y = element_line(size = 0.25))
#gg <- gg + scale_fill_manual(values = ggplot2::alpha(overalcoloruses, alpha = dot.alpha), drop = FALSE) + guides(fill = FALSE)
gg <- gg + scale_size_continuous(limits = weight.MinMax, range = dot.size)+ guides(fill=FALSE, color=FALSE)
gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels, colour = labels), size = 0.0001, show.legend = F, segment.size = 0.2, segment.alpha = 0.5)
gg<-gg+geom_path(aes(colour=labels,size = 1.5),arrow = arrow(type = "closed",length=unit(0.2, "inches")))
#gg<-gg+geom_path(arrow = arrow(length=unit(0.15,"cm"), ends="last", type = "closed"))+
gg <- gg + scale_colour_manual(values =overalcoloruses, drop = FALSE) + guides(colour = FALSE)+ guides(fill = FALSE)
#gg + theme(legend.position = "none")
pdf(paste0(outFolder,"outgoing_incoming_conditions_nolabel_both.pdf"),width=12,height=8)
gg+ theme(legend.position = "none")

dev.off()










