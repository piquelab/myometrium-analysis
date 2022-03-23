library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Matrix)
library(tidyverse)
library(Seurat)


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



#######################
# Location- specific analysis
#######################

res <- read_tsv("7_outputs_DESeq_ConditionsByCluster_bath_library/ALL.combined.2021-08-30.tsv")

res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
res <- res %>% filter(!is.na(pvalue))
res<-res %>% filter(padj<0.1)
clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
names(clust2Names)<-c(0:23)
res$Cell_type<-clust2Names[res$Cell_type]

cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D","#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48","#AA5485","#4E747A","#C59A89","#C9C76F")   
names(cluster.Colors)<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte", "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast", "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")






future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)

# CellChat requires two user inputs: 
# one is the gene expression data of cells, 
# and the other is either user assigned cell labels (i.e., label-based mode) 
# or a low-dimensional representation of the single-cell data

####################################
# Load data
####################################




# Here we load a scRNA-seq data matrix and its associated cell meta data


sc <- read_rds("6_harmony_cellClass_plots_res0.8_final/SeuratObject.rds")
sc2<-read_rds("4_harmony/sc.NormByLibrary.Harmony.StringentFiltering.res0.8.rds")

sc2 <- NormalizeData(sc2, verbose=TRUE)


cmd <- paste0("zcat ",
              "/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_31/gencode.v31.annotation.gff3.gz",
              " | awk '$3~/gene/'",
              " | sed 's/ID=.*gene_id=//;s/;gene_type=/\\t/;s/;gene_name=/\\t/;s/;.*//'")
cat(cmd,"\n")

## Check 0 or 1-based coordinates. 

aux <- data.table::fread(cmd=cmd) %>% mutate(TSS=ifelse(V7=="+",V4,V5)) %>%
  select(Chr=V1,Min=V4,Max=V5,kbid=V9,TSS,Strand=V7,Type=V10,gene_name=V11) 
##anno <- tibble(kbid=rownames(sc)) %>% mutate(ensgene=gsub("\\..*","",kbid)) %>% left_join(grch38)

anno <- tibble(kbid=rownames(sc2),rs=rowSums(sc2@assays$RNA@data)) %>% filter(rs>0) %>% left_join(aux) %>%
  filter(!is.na(Chr))

sc2 <- sc2[anno$kbid,]



locations<-unique(sc2$Location)

#locations<-"Myometrium"
#locations<-"Myometrium"
conditions<-unique(sc2$Condition)
#threshold<-100
threshold<-50

mclapply(conditions, function(xcondition)

  {
   
  
  outFolder=paste0("./10_CellChat_analysis_default_after_filter_",threshold,"/")
  system(paste0("mkdir -p ", outFolder))
  

   sc2_filter<-subset(sc2, Condition==xcondition) 
    
    
    metadata<-sc@meta.data
    mapping<-metadata$seurat_clusters
    names(mapping)<-rownames(metadata)
    
    
    sc2_filter$barcode<-colnames(sc2_filter)
    sc2_filter$seurat_clusters<-mapping[sc2_filter$barcode]
    sc2_filter$celltype<-clust2Names[sc2_filter$seurat_clusters]
    cluster_filter<-names(which(table(sc2_filter$seurat_clusters)>threshold))
    
    

    metadata<-metadata %>% filter (Condition==xcondition) 
    
    
    metadata<-metadata %>% filter( seurat_clusters %in% cluster_filter )
    sc2_filter<-subset(sc2_filter,barcode %in% rownames(metadata) & seurat_clusters %in% cluster_filter) 
    
    metadata$labels<-clust2Names[metadata$seurat_clusters]
    data<-sc2_filter@assays$RNA@data
    
    convertensymbol<-anno$gene_name
    names(convertensymbol)<-anno$kbid
    rownames(data)<-convertensymbol[rownames(data)]
    
    
    
    
    cellchat <- createCellChat(object = data, meta = metadata, group.by = "labels")
    cellchat <- addMeta(cellchat, meta = metadata)
    cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
    levels(cellchat@idents) # show factor levels of the cell labels
    groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
    
    
    ####################################
    # Set the ligand-receptor interaction database
    ####################################
    
    CellChatDB <- CellChatDB.human
    
    # pdf(paste0(outFolder,"showDatabaseCategory.pdf"))
    # showDatabaseCategory(CellChatDB)
    # dev.off()
    
    # Show the structure of the database
    dplyr::glimpse(CellChatDB$interaction)
  
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
    
    write_rds(cellchat,file=paste0(outFolder,"cellchat_",xcondition,"_",Sys.Date(),".rds"))
    
},mc.cores=6)
  










