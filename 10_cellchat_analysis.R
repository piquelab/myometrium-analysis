library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Matrix)
library(tidyverse)
library(Seurat)
library(ggalluvial)


################################################################################
# Cell-cell communication using CellChat
################################################################################



outFolder="./10_CellChat_analysis_customized/"
system(paste0("mkdir -p ", outFolder))


future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)


########################################
# ligand-receptor for contraction 
########################################
lr_contraction<-read_tsv("LR-contraction.txt")
rownames(lr_contraction)<-lr_contraction$interaction_name
write.csv(lr_contraction,file=paste0(outFolder,"lr_contraction.csv"))

##################################################################
# customized CellChatDB.human
##################################################################
old<-CellChatDB.human$interaction
customized<-rbind(lr_contraction,old)
customized<-as.data.frame(customized)
CellChatDB.human$interaction<-customized



##################################################################
# load single cell data (labor-associated DEGs)
##################################################################

res <- read_tsv("7_outputs_DESeq_ConditionsByCluster_bath_library/ALL.combined.2021-08-30.tsv")

res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
res <- res %>% filter(!is.na(pvalue))
res<-res %>% filter(padj<0.1)
clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
names(clust2Names)<-c(0:23)
res$Cell_type<-clust2Names[res$Cell_type]


##################################################################
# cluster colors
##################################################################

cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D","#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48","#AA5485","#4E747A","#C59A89","#C9C76F")   
names(cluster.Colors)<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte", "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast", "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")



# CellChat requires two user inputs: 
# one is the gene expression data of cells, 
# and the other is either user assigned cell labels (i.e., label-based mode) 
# or a low-dimensional representation of the single-cell data

####################################
# Load single cell data
####################################


# Here we load a scRNA-seq data matrix and its associated cell meta data

#sc <- read_rds("4_harmony_cellClass_soupx_doubletfinder_chrM/sc.NormByLibrary.cellclassify_newfilter-res0.5.2021-06-28.rds")
sc <- read_rds("6_harmony_cellClass_plots_res0.8_final/SeuratObject.rds")
sc2<-read_rds("4_harmony/sc.NormByLibrary.Harmony.StringentFiltering.res0.8.rds")

sc2 <- NormalizeData(sc2, verbose=TRUE)
data.input<-sc2@assays$RNA@data
meta<-sc@meta.data
data.input<-data.input[,rownames(meta)]

meta$labels<-clust2Names[meta$seurat_clusters]


anno <- read_rds("3_MergeDemux_Output/anno.rds")

gene_symbol <- anno$gene_name
names(gene_symbol) <- anno$kbid

rownames(data.input)<-gene_symbol[rownames(data.input)]


####################################
# Create a CellChat object
####################################
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")


####################################
# Add cell information into meta slot of the object (Optional)
####################################

cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

####################################
# Set the ligand-receptor interaction database
####################################

CellChatDB <- CellChatDB.human


# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)


# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB


CellChatDB.use <- CellChatDB 


# set the used database in the object
cellchat@DB <- CellChatDB.use

########################################################################
# Preprocessing the expression data for cell-cell communication analysis
########################################################################

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database



cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

# Part II: Inference of cell-cell communication network


# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)



########################################################################
# Infer the cell-cell communication at a signaling pathway level
########################################################################

# NB: The inferred intercellular communication network of each ligand-receptor pair and each signaling pathway is stored in the slot ‘net’ and ‘netP’, respectively.
cellchat <- computeCommunProbPathway(cellchat)


########################################################################
# Calculate the aggregated cell-cell communication network
########################################################################
cellchat <- aggregateNet(cellchat)

save(cellchat,file=paste0(outFolder,"cellchat.RData"))
#######################################################################
# Visualization
#######################################################################

outFolder="./10_CellChat_analysis_customized/"
load(paste0(outFolder,"cellchat.RData"))
groupSize <- as.numeric(table(cellchat@idents))

# pdf(paste0(outFolder,"interaction_top10.pdf"),width=10,height=10)
# netVisual_circle(cellchat@net$weight,remove.isolate=TRUE,color.use=cluster.Colors[rownames(cellchat@net$weight)], vertex.weight = groupSize, weight.scale = T, title.name = "Top 10 % interactions (myometrium)",vertex.label.cex = 1,arrow.width = 6)
# dev.off()
# 
# 
# weight<-cellchat@net$weight
# weight.vector<-as.vector(weight)
# weight.vector<-weight.vector[order(weight.vector,decreasing = TRUE)]
# q<-weight.vector[round(length(weight.vector)*.05)]
# 
# mat<-weight
# mat[mat<q]<-0
# isolated<-rownames(mat) [intersect(which(rowSums(mat)==0)  , which(colSums(mat)==0))]
# 
# mat<-mat[which(!rownames(mat) %in% isolated),which(!rownames(mat) %in% isolated)]
# 
# pdf(paste0(outFolder,"interaction_top10percent.pdf"),width=10,height=10)
# netVisual_circle(mat,color.use=cluster.Colors[rownames(mat)], vertex.weight = groupSize[which(rownames(mat) %in% rownames(weight))], weight.scale = T, title.name = "Top 5 % interactions (myometrium)",vertex.label.cex = 1,arrow.width = 6)
# dev.off()


#######################################################
# signaling roles heatmap plots
#######################################################

pathways<-cellchat@netP$pathways
system(paste0("mkdir -p ", outFolder,"Pathways/"))
sapply(pathways,function(pathways.show){
  #pathways.show <- c("CXCL") 
  # Hierarchy plot
  # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  pdf(paste0(outFolder,"Pathways/",pathways.show,".pdf"),width=25,height=6)
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 25, height = 6, font.size = 10)
  
  dev.off()
})


#######################################################

# Top 10 % interactions per pathway (p <0.05)
# Circle plots
#######################################################


pathways<-cellchat@netP$pathways
system(paste0("mkdir -p ", outFolder,"Pathways_circleplots_top10/"))


sapply(pathways,function(pathways.show){
  #pathways.show <- c("CXCL") 
  # Hierarchy plot
  # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
  #cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  pdf(paste0(outFolder,"Pathways_circleplots_top10/",pathways.show,".pdf"),width=25,height=25)
  par(mfrow=c(1,1))
  #arrow(length = unit(.02, "inches"),type = "closed",angle = 40)
  netVisual_aggregate(cellchat, signaling = pathways.show[1], layout = "circle",color.use=cluster.Colors[rownames(cellchat@net$weight)],vertex.label.cex=2,top=0.10,arrow.width = 15)
  dev.off()
})



# for manuscript, hiding the labels 
pathways<-cellchat@netP$pathways
system(paste0("mkdir -p ", outFolder,"Pathways_circleplots_white_labels_top10/"))
sapply(pathways,function(pathways.show){
  pdf(paste0(outFolder,"Pathways_circleplots_white_labels_top10/",pathways.show,".pdf"),width=25,height=25)
  par(mfrow=c(1,1))
  #arrow(length = unit(.02, "inches"),type = "closed",angle = 40)
  netVisual_aggregate(cellchat, signaling = pathways.show[1], layout = "circle",color.use=cluster.Colors[rownames(cellchat@net$weight)],vertex.label.cex=0.0001,top=0.10,arrow.width = 15,vertex.label.color="white")
  dev.off()
})




# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways


# cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# pdf(paste0(outFolder,"outgoing_incoming.pdf"),width=50,height=20)
# ht1 + ht2
# dev.off()
# 
# pdf(paste0(outFolder,"outgoing.pdf"),width=35,height=70)
# ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",width=35,height=70,font.size=15,font.size.title=30)
# ht1
# dev.off()
# 
# pdf(paste0(outFolder,"incoming.pdf"),width=35,height=70)
# ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",width=35,height=70,font.size=15,font.size.title=30)
# ht2
# dev.off()
# 
# 
# #### 
# cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# pdf(paste0(outFolder,"outgoing_interfron_collagen.pdf"),width=35,height=6)
# ht1 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("COLLAGEN","IFN-II"),color.use=cluster.Colors[rownames(cellchat@net$weight)],pattern = "outgoing",width=35,height=6,font.size=15,font.size.title=30)
# ht1
# dev.off()
# 
# pdf(paste0(outFolder,"incoming_interfron_collagen.pdf"),width=35,height=6)
# ht2 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("COLLAGEN","IFN-II"),color.use=cluster.Colors[rownames(cellchat@net$weight)],pattern = "incoming",width=35,height=6,font.size=15,font.size.title=30)
# ht2
# dev.off()
# 
# 
# cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# pdf(paste0(outFolder,"outgoing_interfron_collagen_v2.pdf"),width=35,height=6)
# ht1 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("COLLAGEN","IFN-II","IL6","IL2","IL1"),color.use=cluster.Colors[rownames(cellchat@net$weight)],pattern = "outgoing",width=35,height=6,font.size=15,font.size.title=30)
# ht1
# dev.off()
# 
# pdf(paste0(outFolder,"incoming_interfron_collagen_v2.pdf"),width=35,height=6)
# ht2 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("COLLAGEN","IFN-II","IL6","IL2","IL1"),color.use=cluster.Colors[rownames(cellchat@net$weight)],pattern = "incoming",width=35,height=6,font.size=15,font.size.title=30)
# ht2
# dev.off()
# 
# 
# cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# pdf(paste0(outFolder,"allsingnals_interfron_collagen_v2.pdf"),width=35,height=6)
# ht1 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("COLLAGEN","IFN-II","IL6","IL2","IL1"),color.use=cluster.Colors[rownames(cellchat@net$weight)],pattern = "all",width=35,height=6,font.size=15,font.size.title=30)
# ht1
# dev.off()
# 
# pdf(paste0(outFolder,"circleplot_interfron_collagen_v2.pdf"),width=25,height=25)
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat, signaling =c("COLLAGEN","IFN-II","IL6","IL2","IL1"), layout = "circle",color.use=cluster.Colors[rownames(cellchat@net$weight)],vertex.label.cex=2,top=0.1,arrow.width = 8)
# dev.off()

###################################################################################
# ligand receptor information per pathways
###################################################################################


pathways<-cellchat@netP$pathways
system(paste0("mkdir -p ", outFolder,"Pathways_LR/"))
sapply(pathways,function(pathways.show){
  pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
  LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
  # Hierarchy plot
  vertex.receiver = seq(1,10) # a numeric vector
  pdf(paste0(outFolder,"Pathways_LR/",pathways.show,".pdf"),width=25,height=25)
  par(mfrow=c(1,1))
  netVisual_individual(cellchat, signaling = pathways.show,  color.use=cluster.Colors[rownames(cellchat@net$weight)],pairLR.use = LR.show, vertex.receiver = vertex.receiver,layout="hierarchy")
  
  dev.off()
})


#pathways.show<-"COLLAGEN"
#pathways.show<-"IFN-II"
#pathways.show<-"IL6"
#pathways.show<-"IL2"
#pathways.show<-"IL1"
#pathways.show<-c("COLLAGEN","IFN-II","IL6","IL1")



pathways.show<-cellchat@netP$pathways
LR_genes_list<-lapply(pathways.show,function(x){
  
  pairLR <- extractEnrichedLR(cellchat, signaling = x, geneLR.return = FALSE)
  LR_genes<-unique(unlist(str_split(pairLR$interaction_name,"_")))
  return(LR_genes)
})


res_list<-lapply(pathways.show,function(x){
  
  pairLR <- extractEnrichedLR(cellchat, signaling = x, geneLR.return = FALSE)
  LR_genes<-unique(unlist(str_split(pairLR$interaction_name,"_")))
  
  res_LR_genes<-res %>% filter(gene_name %in% LR_genes)
  res_LR_genes$pathway_name<-x
  write.csv(res_LR_genes,file=paste0(outFolder,x,"_","res_LR_genes.csv"))
  
  return(res_LR_genes)
})

res_all<-do.call(rbind,res_list)
write.csv(res_all,file=paste0(outFolder,"res_LR_pathways.csv"))




#################################################################################

# Making data frame, including pathways, ligand/receptors genes that are involved in those pathways and labor-associated DEGs , cell types with those DEGs, and their rowLogSumExps()


#################################################################################
pathways<-cellchat@netP$pathways
res <- read_tsv("7_outputs_DESeq_ConditionsByCluster_bath_library/SIG.combined.2021-02-17.tsv")
res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
res <- res %>% filter(!is.na(pvalue))
res<-res %>% filter(padj<0.1)
clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
names(clust2Names)<-c(0:23)
res$Cell_type<-clust2Names[res$Cell_type]

cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D","#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48","#AA5485","#4E747A","#C59A89","#C9C76F")   
names(cluster.Colors)<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte", "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast", "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")



load(paste0(outFolder,"cellchat.RData"))
groupSize <- as.numeric(table(cellchat@idents))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

################################
# pathway / ligand/receptor from pathways that are labor-associated DEGs / celltype with DEGs
################################
pathways_genes_celltypes<-lapply(pathways,function(x){
  
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

write.csv(pathways_genes_celltypes_db,file = paste0(outFolder,"pathways_genes_celltypes_db.csv"))


################################
# cells with sending role (from cellchat)
################################
senders<-lapply(pathways,function(x){
  
  scores<-cellchat@netP$centr[[x]]$outdeg
  minx<-min(scores)
  maxx<-max(scores)
  scores_scaled<-sapply( scores, function(x){(x-minx)/(maxx-minx)})
  senders<-names(scores_scaled)[which(scores_scaled>=0.75)]
  roles<-rep("sender",length(senders))
  path<-rep(x,length(senders))
  return(cbind(senders,roles,path))
})
senders_df<-do.call(rbind,senders)



################################
# cells with receiving role (from cellchat)
################################

receiver<-lapply(pathways,function(x){
  
  scores<-cellchat@netP$centr[[x]]$indeg
  minx<-min(scores)
  maxx<-max(scores)
  scores_scaled<-sapply( scores, function(x){(x-minx)/(maxx-minx)})
  senders<-names(scores_scaled)[which(scores_scaled>=0.75)]
  roles<-rep("receiver",length(senders))
  path<-rep(x,length(senders))
  return(cbind(senders,roles,path))
})
receiver_df<-do.call(rbind,receiver)


################################
# merging sender/receiever cells
################################



receiver_senders_df<-rbind(senders_df,receiver_df)
colnames(receiver_senders_df)<-c("celltype","role","pathway")
receiver_senders_df<-receiver_senders_df[!duplicated(receiver_senders_df), ]
#receiver_senders_df<-receiver_senders_df[order(receiver_senders_df$pathway,decreasing = FALSE),]
receiver_senders_df<-as.data.frame(receiver_senders_df)
receiver_senders_df<-receiver_senders_df %>% arrange(pathway)

#write.csv(receiver_senders_df,file=paste0(outFolder,"receiver_senders_df.csv"))
write.csv(receiver_senders_df,file=paste0(outFolder,"receiver_senders_0.75_df.csv"))

sender_receiver_columns<-sapply(1: nrow(pathways_genes_celltypes_db), function(x){
  path<-pathways_genes_celltypes_db$Pathway[x]
  celltypes<-pathways_genes_celltypes_db$`Cell type`[x]
  celltypes<-unlist(strsplit(celltypes,", ",fixed=TRUE))
  # senders<- receiver_senders_df %>% filter (pathway==path & celltype %in% celltypes & role=="sender") %>% select(celltype)%>% unlist %>% unique
  # receivers<-receiver_senders_df %>% filter (pathway==path & celltype %in% celltypes & role=="receiver") %>% select(celltype)%>% unlist %>% unique
  senders<- receiver_senders_df %>% filter (pathway==path &  role=="sender") %>% select(celltype)%>% unlist %>% unique
  receivers<-receiver_senders_df %>% filter (pathway==path & role=="receiver") %>% select(celltype)%>% unlist %>% unique
  
  senders<-paste(senders,collapse = ", ")
  receivers<-paste(receivers,collapse = ", ")
  if (senders=="") senders<-"NA"
  if (receivers=="") receivers<-"NA"
  return(c(senders,receivers))
})



################################
# make the final data frame 
################################



sender_receiver_columns<-t(sender_receiver_columns)
pathways_genes_celltypes_db<-cbind(pathways_genes_celltypes_db,sender_receiver_columns)
colnames(pathways_genes_celltypes_db)[c(4,5)]<-c("senders","receivers")
pathways_genes_celltypes_roles_db<-pathways_genes_celltypes_db
colnames(pathways_genes_celltypes_roles_db)[2]<-"Genes in pathway that are labor-associated DEGs"
colnames(pathways_genes_celltypes_roles_db)[3]<-"Cell type with DEGs (genes from pathway)"
colnames(pathways_genes_celltypes_roles_db)[c(4,5)]<-c("senders-cellchat"      ,    "receivers-cellchat"     )
#write.csv(pathways_genes_celltypes_roles_db,file = paste0(outFolder,"pathways_genes_celltypes_roles_db.csv"))
write.csv(pathways_genes_celltypes_roles_db,file = paste0(outFolder,"pathways_DEgenes_celltypes_roles_db.csv"))




###########################################################################################
# Finding outgoing and incoming patterns 
###########################################################################################
library(NMF)
library(ggalluvial)


outFolder="./10_CellChat_analysis_customized/"
pathways.show<-c("contraction","COLLAGEN","LAMININ","MIF","GALECTIN","CCL","IL1","COMPLEMENT","PDGF","TNF","THY1","TGFb","IFN-II","IL6","IL2")

load(paste0(outFolder,"cellchat.RData"))
groupSize <- as.numeric(table(cellchat@idents))

# selected pathways
pathways.show<-c("contraction","COLLAGEN","LAMININ","MIF","GALECTIN","CCL","IL1","COMPLEMENT","PDGF","TNF","THY1","TGFb","IFN-II","IL6","IL2")


# finding the best number of patterns based on cophenetic and silouette metrics

# choozing best number of outgoing patterns
# pdf(paste0(outFolder,"n_pattern_outgoing.pdf"),width=25,height=10)
# selectK(cellchat, pattern = "outgoing")
# dev.off()

# nPatterns=3 is chosen based on previous plot 
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

# outgoing patterns
pdf(paste0(outFolder,"outgoing_selected_riverplot_pattern3.pdf"),width=30,height=20)
netAnalysis_river(cellchat, pattern = "outgoing",color.use=cluster.Colors[rownames(cellchat@net$weight)],signaling=pathways.show,font.size = 5,main.title = NULL)
dev.off()


# choozing best number of incoming patterns
# pdf(paste0(outFolder,"n_pattern_incoming.pdf"),width=25,height=10)
# selectK(cellchat, pattern = "incoming")
# dev.off()

# nPatterns=3 is chosen based on previous plot 
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)


pdf(paste0(outFolder,"incoming_selected_riverplot_pattern3.pdf"),width=30,height=20)
netAnalysis_river(cellchat, pattern = "incoming",color.use=cluster.Colors[rownames(cellchat@net$weight)],signaling=pathways.show,font.size = 5,main.title = NULL)
dev.off()



# pdf(paste0(outFolder,"outgoing_riverplot_pattern3.pdf"),width=30,height=15)
# netAnalysis_river(cellchat, pattern = "outgoing",color.use=cluster.Colors[rownames(cellchat@net$weight)],signaling=c("COLLAGEN","IFN-II","IL6","IL1"),font.size = 5,main.title = NULL)
# dev.off()
# 
# pdf(paste0(outFolder,"incoming_riverplot_pattern3.pdf"),width=30,height=15)
# netAnalysis_river(cellchat, pattern = "incoming",color.use=cluster.Colors[rownames(cellchat@net$weight)],signaling=c("COLLAGEN","IFN-II","IL6","IL1"),font.size = 5,main.title = NULL)
# dev.off()

# pdf(paste0(outFolder,"outgoing_all_riverplot_pattern3.pdf"),width=30,height=50)
# netAnalysis_river(cellchat, pattern = "outgoing",color.use=cluster.Colors[rownames(cellchat@net$weight)],font.size = 5,main.title = NULL)
# dev.off()



# pdf(paste0(outFolder,"incoming_all_riverplot_pattern3.pdf"),width=30,height=50)
# netAnalysis_river(cellchat, pattern = "incoming",color.use=cluster.Colors[rownames(cellchat@net$weight)],font.size = 5,main.title = NULL)
# dev.off()


# Threshold for contribution 
#cutoff<-0.5
cutoff<-0.7



outgoing_signaling<-cellchat@netP$pattern$outgoing$pattern$signaling
outgoing_signaling<-outgoing_signaling %>% filter(Contribution>cutoff)
#outgoing_signaling$Contribution[outgoing_signaling$Contribution < cutoff] <- 0


colnames(outgoing_signaling)[3]<-"Contribution_outgoing_signaling"
colnames(outgoing_signaling)[1]<-paste0(colnames(outgoing_signaling)[1],"_outgoing")

outgoing_cell<-cellchat@netP$pattern$outgoing$pattern$cell
outgoing_cell<-outgoing_cell %>% filter(Contribution>cutoff)
#outgoing_cell$Contribution[outgoing_cell$Contribution < cutoff] <- 0
colnames(outgoing_cell)[3]<-"Contribution_outgoing_cell"
colnames(outgoing_cell)[2]<-paste0(colnames(outgoing_cell)[2],"_outgoing")
outgoing_data<-outgoing_signaling %>% inner_join(outgoing_cell)
outgoing_data<-outgoing_data %>% filter(Signaling %in% pathways.show)




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



# Joining incoming and outgoing 
incoming_data<-incoming_signaling %>% inner_join(incoming_cell)
incoming_data<-incoming_data %>% filter(Signaling %in% pathways.show)


outgoing_data<-outgoing_data %>% mutate(CellGroup_outgoing=CellGroup)#,Signaling,Contribution-outgoing_cell,Contribution_outgoing_signaling)
incoming_data<-incoming_data %>% mutate(CellGroup_incoming=CellGroup)#,Signaling,Contribution_incoming_signal,Contribution_incoming_cell)

signalins_incoming<-unique(incoming_data$Signaling)
outgoing_incoming<-unique(outgoing_data$Signaling)
overlap_signalings<-intersect(signalins_incoming,outgoing_incoming)



########### 
# ploting patterns 
# refer to 10_alluvial_plot.R (the same code)

library(ggalluvial)
outgoing_data<- outgoing_data %>% filter(Signaling %in% overlap_signalings)
outgoing_data$Signaling<-as.character(outgoing_data$Signaling)
outgoing_data<-outgoing_data[order(outgoing_data$Signaling,decreasing = FALSE),]
outgoing_data$Signaling <- factor(outgoing_data$Signaling,levels=unique(outgoing_data$Signaling))


fname=paste0(outFolder,"outgoing_ggalluvial_cuttoff.0.7.pdf");
pdf(fname,width=20,height=28)
ggplot(data = outgoing_data,
       aes(axis1 = CellGroup_outgoing, axis2 = Signaling, y = Contribution_outgoing_signaling,fill=CellGroup_outgoing)) +
  geom_alluvium(aes(fill=CellGroup_outgoing)) +#aes(fill = Signaling)
  geom_stratum() +
  geom_flow()+
  scale_fill_manual("Cell type",values=cluster.Colors) +
  geom_text(stat = "stratum",aes(label = after_stat(stratum)),size=5) +
  scale_x_discrete(limits = c("Outgoing cell", "Pathway"),
                   expand = c(0.3, 0.1)) +
  theme_bw()
dev.off()



incoming_data<- incoming_data %>% filter(Signaling %in% overlap_signalings)
incoming_data$Signaling<-as.character(incoming_data$Signaling)
incoming_data<-incoming_data[order(incoming_data$Signaling,decreasing = FALSE),]
incoming_data$Signaling <- factor(incoming_data$Signaling,levels=unique(incoming_data$Signaling))



fname=paste0(outFolder,"incoming_ggalluvial.cuttof.0.7.pdf");
pdf(fname,width=20,height=28)
ggplot(data = incoming_data,
       aes(axis1 = Signaling, axis2 = CellGroup_incoming, y = Contribution_incoming_signaling)) +
  geom_alluvium() +#aes(fill = Signaling)
  geom_stratum(aes(fill=CellGroup_incoming)) +
  geom_flow(aes(fill=CellGroup_incoming))+
  scale_fill_manual("Cell type",values=cluster.Colors) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),size=5) +
  scale_x_discrete(limits = c("Pathway", "Incoming cell"),
                   expand = c(0.15, 0.05)) +
  theme_bw()
dev.off()



outgoing_data<-outgoing_data[,-4]
incoming_data<-incoming_data[,-4]

outgoing_incoming_data<-outgoing_data %>% inner_join(incoming_data)

outgoing_incoming_data$Contribution_outgoing_signaling<-1
outgoing_incoming_data$Contribution_incoming_signaling<-1
#


library(ggalluvial)

fname=paste0(outFolder,"incoming_outgoing_ggalluvial.pdf");
pdf(fname,width=20,height=28)

ggplot(data = outgoing_incoming_data,
       aes(axis1 = CellGroup_outgoing,   # First variable on the X-axis
           axis2 = Signaling, # Second variable on the X-axis
           axis3 = CellGroup_incoming,   # Third variable on the X-axis
           y = Contribution_outgoing_signaling)) +  # what about incoming
  geom_alluvium(aes(fill = Signaling)) +
  geom_stratum() +
  scale_fill_manual("Cell type",values=cluster.Colors) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),size=5) +
  scale_x_discrete(limits = c("CellGroup_outgoing", "CellGroup_incoming"),
                   expand = c(0.15, 0.05)) +
  theme_bw()
dev.off()



