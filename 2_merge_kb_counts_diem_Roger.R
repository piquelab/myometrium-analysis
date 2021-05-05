## options(repos = c(CRAN = "http://cran.rstudio.com"))
##install.packages("tidyverse")
## library("rhdf5")

library(Seurat)

library(Matrix)

library(tidyverse)

library(diem)

library(SeuratDisk)

## Load the datasets
##basefolder <- "/nfs/rprdata/scilab/cc1-Placenta/counts/ALL-recount-6.17.19-preRNA/"
##basefolder <- "/nfs/rprdata/scilab/endometrium/kallisto2/bus/"
##expNames <- dir(basefolder,"^s*")

#basefolder <- "/nfs/rprdata/scilab/labor2/kallisto.covid19/bus/"
#expNames <- dir(basefolder,"^C19C*")

basefolder <- "/wsu/home/groups/prbgenomics/parturition/kallisto/bus/"
expNames <- dir(basefolder)

expNames

##subfolder <- "/outs/filtered_gene_bc_matrices/hg19/"
##subfolder <- "/outs/raw_gene_bc_matrices/hg19/"
##subfolder <- "/outs/"
subfolder <- ""
folders <- paste0(basefolder,expNames,subfolder)
ind <- file.info(folders)$isdir
ind[is.na(ind)]<- FALSE
folders <- folders[ind]
expNames <- expNames[ind]
names(folders) <- expNames
folders

##aa2 <- ReadH5AD(file = "../kallisto2/bus/HPL20289_R_C_1/counts_filtered/adata.h5ad",verbose=TRUE)

ii<-expNames[1]

adata <- map(expNames,function(ii){
    ##
    
    expPrefix = ii;
    folders[ii]
    fName= paste0(folders[ii],"/counts_unfiltered/adata.h5ad")
    cat("#Loading ",fName,"... \n")
    ##adata = ReadH5AD(file = fName,verbose=TRUE)

    fName2= paste0(folders[ii],"/counts_unfiltered/adata.h5seurat")
##    Convert(fName,fName2,overwrite = TRUE)
    adata <- LoadH5Seurat(fName2)
    
    cat(dim(adata),"\n")
    ## May need to rename rownames or split...
    sce <- create_SCE(rbind(adata@assays$spliced@data,adata@assays$unspliced@data))
    cat(dim(sce),"\n")
    
    ## Remove debris...  And this seems to have changed ... or consider switching to SoupX. 
    sce <- diem(sce)
    ##
    sc = convert_to_seurat(sce)
    cat("#Final: ",dim(sc),"\n")
    adata <- adata[,colnames(sc)]
    ## I could keep the background debris model, or also return cells not debris
    ##  adata <- SCTransform(adata, verbose = TRUE)
    adata
})


names(adata) <- expNames


system("mkdir -p 2_kb_diem_Output/")
write_rds(adata,"2_kb_diem_Output/kb_diem_Seurat.list.rds")

sparse.size <- object.size(adata)
sparse.size


sc <- merge(adata[[1]], y = adata[-1], project = "KallistoDiem",add.cell.ids=expNames)

dim(sc)

system("mkdir -p 2_kb_diem_Output/")
write_rds(sc,"2_kb_diem_Output/kb_diem_Seurat.obj.rds")

##############################
