library(Seurat)
library(Matrix)
library(tidyverse)
library(diem)


#############################################################
# Downstream analysis of the kallisto output
#############################################################


rm(list=ls())
outdir <- "./2_kb_Output/"
system("mkdir -p 2_kb_Output/")


################################################
# generate folders containing h5ad data produced by kallisto 
################################################
basefolder <- "/wsu/home/groups/prbgenomics/labor_myo/kallisto/bus/"
expNames <- dir(basefolder,"^Labor*")

expNames

subfolder <- ""
folders <- paste0(basefolder,expNames,subfolder)
ind <- file.info(folders)$isdir
ind[is.na(ind)]<- FALSE
folders <- folders[ind]
expNames <- expNames[ind]
names(folders) <- expNames
folders


################################################
# read h5ad data into seurat
################################################


readKallisto  <- function (run, prefixFile="spliced/s",expPrefix=NULL) 
{
  if (!dir.exists(paths = run)) {
    stop("Directory provided does not exist")
  }
  barcode.loc <- file.path(run, paste0(prefixFile,".barcodes.txt"))
  gene.loc <- file.path(run, paste0(prefixFile,".genes.txt"))
  ##    features.loc <- file.path(run, "features.tsv.gz")
  matrix.loc <- file.path(run, paste0(prefixFile,".mtx"))
  if (!file.exists(barcode.loc)) {
    stop("Barcode file missing")
  }
  if (!file.exists(gene.loc)) {
    stop("Gene name or features file missing")
  }
  if (!file.exists(matrix.loc)) {
    stop("Expression matrix file missing")
  }
  
  # Read in the count matrix that was output by `kb`.
  data <- readMM(file = matrix.loc)
  
  cell.names <- readLines(barcode.loc)
  if(is.null(expPrefix)){
    rownames(data) <- cell.names
  }
  else{
    rownames(data) <- paste0(expPrefix,"_",cell.names)
  }    
  feature.names <- readLines(gene.loc)
  colnames(data) <- feature.names
  
  # The matrix read has cells as rows
  t(data)
}

###
adata <- sapply(expNames,function(ii){
  ##
  expPrefix = ii;
  cat("#Loading ",paste0(folders[ii], "/counts_unfiltered/spliced"), " ...")    
  sFull <- readKallisto(folders[ii], prefixFile="/counts_unfiltered/spliced" , expPrefix)  #"/spliced/s"
  cat(dim(sFull),"\n")
  cat("#Loading ",paste0(folders[ii], "/counts_unfiltered/unspliced"), " ...")
  uFull <- readKallisto(folders[ii], prefixFile="/counts_unfiltered/unspliced", expPrefix) #"/unspliced/u"
  cat(dim(uFull),"\n")
  ##
  scs <- colSums(sFull)
  #rownames(sFull) <- paste0("S-",rownames(sFull))
  ucs <- colSums(uFull)
  #rownames(uFull) <- paste0("U-",rownames(uFull))

  sel <- intersect(colnames(sFull),colnames(uFull))
  sel <- intersect(colnames(sFull)[scs>0],colnames(uFull)[ucs>0])
  count0 <- rbind(sFull[,sel],uFull[,sel])
  
  ## May need to rename rownames or split...
  sce <- CreateSeuratObject(count0)
  #sce <- create_SCE(count0)
  cat(dim(sce),"\n")
  ## Remove debris...
  #sce <- diem(sce,top_n = 16000)
  #sc <- convert_to_seurat(sce)
  
  #cat(dim(sce),"\n")
  
  sce
})

write_rds(adata,paste0(outdir,"adata.rds"))
sc <- merge(adata[[1]],adata[-1], project="kbMyometrium")
write_rds(sc,paste0(outdir,"kb_Seurat.obj.rds"))



