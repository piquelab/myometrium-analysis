library(NMF)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Matrix)
library(tidyverse)
library(Seurat)
library(ggalluvial)

################################################################################
# Cell-cell communication using CellChat

# Outgoing and incoming patterns (alluvial plots)
################################################################################


outFolder="./10_CellChat_analysis_customized/"
pathways.show<-c("contraction","COLLAGEN","LAMININ","MIF","GALECTIN","CCL","IL1","COMPLEMENT","PDGF","TNF","THY1","TGFb","IFN-II","IL6","IL2")
nPatterns = 3

# filter out weak links (contribution score)

#cutoff<-0.5
cutoff<-0.7



load(paste0(outFolder,"cellchat.RData"))
groupSize <- as.numeric(table(cellchat@idents))

cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)


# cell type colors
cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D","#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48","#AA5485","#4E747A","#C59A89","#C9C76F")   
names(cluster.Colors)<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte", "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast", "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")

############################################################################
# outgoing data
############################################################################
outgoing_signaling<-cellchat@netP$pattern$outgoing$pattern$signaling
outgoing_signaling<-outgoing_signaling %>% filter(Contribution>cutoff)

colnames(outgoing_signaling)[3]<-"Contribution_outgoing_signaling"
colnames(outgoing_signaling)[1]<-paste0(colnames(outgoing_signaling)[1],"_outgoing")

outgoing_cell<-cellchat@netP$pattern$outgoing$pattern$cell
outgoing_cell<-outgoing_cell %>% filter(Contribution>cutoff)

colnames(outgoing_cell)[3]<-"Contribution_outgoing_cell"
colnames(outgoing_cell)[2]<-paste0(colnames(outgoing_cell)[2],"_outgoing")

outgoing_data<-outgoing_signaling %>% inner_join(outgoing_cell)
outgoing_data<-outgoing_data %>% filter(Signaling %in% pathways.show)


############################################################################
# incoming data
############################################################################
incoming_signaling<-cellchat@netP$pattern$incoming$pattern$signaling
incoming_signaling<-incoming_signaling %>% filter(Contribution>cutoff)
#incoming_signaling$Contribution[incoming_signaling$Contribution < cutoff] <- 0

colnames(incoming_signaling)[3]<-"Contribution_incoming_signaling"
colnames(incoming_signaling)[1]<-paste0(colnames(incoming_signaling)[1],"_incoming")

incoming_cell<-cellchat@netP$pattern$incoming$pattern$cell
incoming_cell<-incoming_cell %>% filter(Contribution>cutoff)
#incoming_cell$Contribution[incoming_cell$Contribution < cutoff] <- 0

colnames(incoming_cell)[3]<-"Contribution_incoming_cell"
colnames(incoming_cell)[2]<-paste0(colnames(incoming_cell)[2],"_incoming")

incoming_data<-incoming_signaling %>% inner_join(incoming_cell)

incoming_data<-incoming_data %>% filter(Signaling %in% pathways.show)

outgoing_data<-outgoing_data %>% mutate(CellGroup_outgoing=CellGroup)#,Signaling,Contribution-outgoing_cell,Contribution_outgoing_signaling)
incoming_data<-incoming_data %>% mutate(CellGroup_incoming=CellGroup)#,Signaling,Contribution_incoming_signal,Contribution_incoming_cell)
signalins_incoming<-unique(incoming_data$Signaling)
outgoing_incoming<-unique(outgoing_data$Signaling)
overlap_signalings<-intersect(signalins_incoming,outgoing_incoming)


############################################################################
# outgoing alluvial plot
############################################################################
outgoing_data<- outgoing_data %>% filter(Signaling %in% overlap_signalings)
outgoing_data$Signaling<-as.character(outgoing_data$Signaling)
outgoing_data<-outgoing_data[order(outgoing_data$Signaling,decreasing = FALSE),]
outgoing_data$Signaling <- factor(outgoing_data$Signaling,levels=unique(outgoing_data$Signaling))

fname=paste0(outFolder,"outgoing_ggalluvial.pdf");
pdf(fname,width=20,height=28)
ggplot(data = outgoing_data,
       aes(axis1 = CellGroup_outgoing, axis2 = Signaling, y = Contribution_outgoing_signaling)) +
  geom_alluvium(aes(fill = Signaling)) +
  geom_stratum() +
  geom_text(stat = "stratum",aes(label = after_stat(stratum)),size=5) +
  scale_x_discrete(limits = c("Outgoing cell", "Pathway"),
                   expand = c(0.3, 0.1)) +
  theme_bw()
dev.off()


# fname=paste0(outFolder,"outgoing_ggalluvial_cuttoff.0.7.pdf");
# pdf(fname,width=20,height=28)
# ggplot(data = outgoing_data,
#        aes(axis1 = CellGroup_outgoing, axis2 = Signaling, y = Contribution_outgoing_signaling,fill=CellGroup_outgoing)) +
#   geom_alluvium(aes(fill=CellGroup_outgoing)) +#aes(fill = Signaling)
#   geom_stratum() +
#   geom_flow()+
#   scale_fill_manual("Cell type",values=cluster.Colors) +
#   geom_text(stat = "stratum",aes(label = after_stat(stratum)),size=5) +
#   scale_x_discrete(limits = c("Outgoing cell", "Pathway"),
#                    expand = c(0.3, 0.1)) +
#   theme_bw()
# dev.off()



############################################################################
# incoming alluvial plot
############################################################################
incoming_data<- incoming_data %>% filter(Signaling %in% overlap_signalings)
incoming_data$Signaling<-as.character(incoming_data$Signaling)
incoming_data<-incoming_data[order(incoming_data$Signaling,decreasing = FALSE),]
incoming_data$Signaling <- factor(incoming_data$Signaling,levels=unique(incoming_data$Signaling))


fname=paste0(outFolder,"incoming_ggalluvial.pdf");
pdf(fname,width=20,height=28)
ggplot(data = incoming_data,
       aes(axis1 = Signaling, axis2 = CellGroup_incoming, y = Contribution_incoming_signaling)) +
  geom_alluvium(aes(fill = Signaling)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),size=5) +
  scale_x_discrete(limits = c("Pathway", "Incoming cell"),
                   expand = c(0.15, 0.05)) +
  theme_bw()
dev.off()




#######################################################
# cellchat table 7
#######################################################

total_signalings<-unique(c(incoming_data$Signaling, outgoing_data$Signaling))

incoming_data_df<-incoming_data %>% select(Signaling,CellGroup_incoming)
outgoing_data_df<-outgoing_data  %>% select(Signaling,CellGroup_outgoing)
data<-incoming_data_df %>% full_join (outgoing_data_df)

senders<-sapply(unique(as.character(data$Signaling)), function(x)
  {
  sender<-data %>% filter (Signaling ==x) %>% select (CellGroup_outgoing) %>% unlist %>% unique
  sender<-as.character(sender)
  sender<-paste(sender,collapse = ", ")
  sender
  
  })


receivers<-sapply(unique(as.character(data$Signaling)), function(x)
  {
  receiever<-data %>% filter (Signaling ==x) %>% select (CellGroup_incoming) %>% unlist %>% unique
  receiever<-as.character(receiever)
  receiever<-paste(receiever,collapse = ", ")
  receiever
  })
# 






groupSize <- as.numeric(table(cellchat@idents))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

################################
# pathway / ligand/receptor from pathways that are labor-associated DEGs / celltype with DEGs
################################

res <- read_tsv("7_outputs_DESeq_ConditionsByCluster_bath_library/SIG.combined.2021-02-17.tsv")
res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
res <- res %>% filter(!is.na(pvalue))
res<-res %>% filter(padj<0.1)
clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
names(clust2Names)<-c(0:23)
res$Cell_type<-clust2Names[res$Cell_type]

cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D","#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48","#AA5485","#4E747A","#C59A89","#C9C76F")   
names(cluster.Colors)<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte", "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast", "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")



pathways_genes_celltypes<-lapply(unique(as.character(data$Signaling)),function(x){
  
  pairLR <- extractEnrichedLR(cellchat, signaling = x, geneLR.return = FALSE)
  LR_genes<-unique(unlist(str_split(pairLR$interaction_name,"_")))
  
  
  celltypes<-res %>% filter(gene_name %in% LR_genes)%>% dplyr::select(Cell_type)%>% unlist %>% unique()
  
  LR_genes <-res %>% filter (gene_name %in% LR_genes) %>% dplyr::select(gene_name)%>% unlist %>% unique()
  LR_genes<-paste(LR_genes,collapse = ", ")
  
  celltypes<-paste(celltypes,collapse = ", ")
  if (celltypes=="") celltypes<-"NA"
  return(c(x,LR_genes,celltypes))
})


pathways_genes_celltypes_db<-do.call(rbind,pathways_genes_celltypes)
colnames(pathways_genes_celltypes_db)<-c("Pathway","Genes in pathway","Cell type")
pathways_genes_celltypes_db<-as.data.frame(pathways_genes_celltypes_db)


final<-cbind(pathways_genes_celltypes_db,senders,receivers)
write.csv(final,file=paste0(outFolder,"cellchat_table.csv"))


# fname=paste0(outFolder,"incoming_ggalluvial.cuttof.0.7.pdf");
# pdf(fname,width=20,height=28)
# ggplot(data = incoming_data,
#        aes(axis1 = Signaling, axis2 = CellGroup_incoming, y = Contribution_incoming_signaling)) +
#   geom_alluvium() +#aes(fill = Signaling)
#   geom_stratum(aes(fill=CellGroup_incoming)) +
#   geom_flow(aes(fill=CellGroup_incoming))+
#   scale_fill_manual("Cell type",values=cluster.Colors) +
#   geom_text(stat = "stratum",
#             aes(label = after_stat(stratum)),size=5) +
#   scale_x_discrete(limits = c("Pathway", "Incoming cell"),
#                    expand = c(0.15, 0.05)) +
#   theme_bw()
# dev.off()

############################################################################
# both incoming and outgoing alluvial plot
############################################################################
# outgoing_data<-outgoing_data[,-4]
# incoming_data<-incoming_data[,-4]
# 
# outgoing_incoming_data<-outgoing_data %>% inner_join(incoming_data)
# 
# outgoing_incoming_data$Contribution_outgoing_signaling<-1
# outgoing_incoming_data$Contribution_incoming_signaling<-1
# #
# 
# 
# library(ggalluvial)
# 
# fname=paste0(outFolder,"incoming_outgoing_ggalluvial.pdf");
# pdf(fname,width=20,height=28)
# 
# ggplot(data = outgoing_incoming_data,
#        aes(axis1 = CellGroup_outgoing,   # First variable on the X-axis
#            axis2 = Signaling, # Second variable on the X-axis
#            axis3 = CellGroup_incoming,   # Third variable on the X-axis
#            y = Contribution_outgoing_signaling)) +  # what about incoming
#   geom_alluvium(aes(fill = Signaling)) +
#   geom_stratum() +
#   geom_text(stat = "stratum",
#             aes(label = after_stat(stratum)),size=5) +
#   scale_x_discrete(limits = c("CellGroup_outgoing", "CellGroup_incoming"),
#                    expand = c(0.15, 0.05)) +
#   theme_bw()
# dev.off()