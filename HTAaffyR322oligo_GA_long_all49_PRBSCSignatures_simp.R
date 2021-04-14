library(lme4)
library(tidyverse)
#load expression data


cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
                  "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48",
                  "#AA5485","#4E747A","#C59A89","#C9C76F")   
names(cluster.Colors)<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
                         "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast",
                         "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
names(cluster.Colors)<-paste0(c(0:23),"_",names(cluster.Colors))

outFolder="./13_Longitudinal_study/"
system(paste0("mkdir -p ", outFolder))



load("CellularTransc260Longit.RData")




# read myometrium single cell sig (cell marker gene set)
SCGeneSets2 = read_tsv("5_harmonyClustersDGE/ClusterDEG.tsv")
colnames(SCGeneSets2)[which(colnames(SCGeneSets2)=="avg_log2FC")]<-"avg_logFC"
colnames(SCGeneSets2)[which(colnames(SCGeneSets2)=="gene")]<-"ENSG"
colnames(SCGeneSets2)[which(colnames(SCGeneSets2)=="symbol")]<-"gene"
SCGeneSets<-SCGeneSets2






clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(c(0:23),"_",clust2Names)
names(clust2Names)<-c(0:23)
SCGeneSets$cluster<-clust2Names[as.character(SCGeneSets$cluster)]
SCGeneSet<-SCGeneSets

### based on DE genes 
res <- read_tsv("7_outputs_DESeq_ConditionsByCluster/SIG.combined.2021-02-17.tsv")
res<-res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
res <- res %>% filter(!is.na(pvalue))
res$Cell_type<-clust2Names[res$Cell_type]
res<-res %>% select("pvalue" , "padj","gene_name" ,"Cell_type" ,"log2FoldChange","kbid")
colnames(res)<-c("p_val","p_val_adj","gene","cluster","avg_logFC","ENSG")

#upregulated signatures 
res_up<-res %>% filter(p_val_adj <= 0.01 & avg_logFC>=0.5)
res_up <- res_up %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)



res_down<-res %>% filter(p_val_adj <= 0.01 & avg_logFC<=-0.5)
res_down<-res_down %>% group_by(cluster) %>% top_n(n = 20, wt = -1*avg_logFC)


res_combined<-res %>% filter(p_val_adj <= 0.01 )
res_combined <- res_combined %>% group_by(cluster) %>% top_n(n = 20, wt = abs(avg_logFC))


#
DE_ref<-read.delim("supp_adi.txt")
DE_ref<-DE_ref %>%filter(adj.P.Val<0.1)

res<-res %>%filter(p_val_adj<0.1)

gene_intersect<-DE_ref$SYMBOL[which( DE_ref$SYMBOL %in% res$gene )]

res_intersect<-res%>%filter(gene %in% gene_intersect)
system(paste0("mkdir -p ", paste0(outFolder,"intersect/")))
write.csv(res_intersect,file=paste0(outFolder,"intersect/res_intersect_ref_Adi.csv"))
write.csv(gene_intersect,file=paste0(outFolder,"intersect/gene_intersect_intersect_ref_Adi.csv"))


#read gene sets
#SCGeneSets=read.table("MarkersAll.txt",sep=" ",header=TRUE,stringsAsFactors = FALSE)

signature_plots(SCGeneSets=res_up,outFolder="./13_Longitudinal_study/up_regulated_sig/")
signature_plots(SCGeneSets=res_down,outFolder="./13_Longitudinal_study/down_regulated_sig/")
signature_plots(SCGeneSets=res_combined,outFolder="./13_Longitudinal_study/combined_sig/")
signature_plots(SCGeneSets=SCGeneSet,outFolder="./13_Longitudinal_study/")

signature_plots<-function(SCGeneSets,outFolder="./13_Longitudinal_study/up_regulated_sig/")
{

  system(paste0("mkdir -p ", outFolder))
  pile=NULL
  for(a in unique(SCGeneSets$cluster)){
    tmp=SCGeneSets[SCGeneSets$cluster==a,]
    tmp=tmp[order(tmp$avg_logFC,decreasing=TRUE),]
    pile=rbind(pile,tmp[1:min(c(20,dim(tmp)[1])),c("gene","cluster")])
  }
  names(pile)<-c("Symbol","Type")
  genes=pile
  #genes=rbind(genes,data.frame(Symbol=genes[genes$Type%in%c("Tcells-activated","Tcells-resting"),"Symbol"],Type="Tcell"))
  #genes=rbind(genes,data.frame(Symbol=genes[genes$Type%in%c("5_CD4_T-cell" ,"7_CD8_T-cell"),"Symbol"],Type="Tcell"))
  
  
  celsig=list()
  for(cs in unique(genes$Type)){
    celsig[[cs]]<-intersect(as.character(genes$Symbol[genes$Type==cs]),aT1$SYMBOL)
  }
  
  myb=b[rownames(aT1),]
  
  esetM=apply(myb[,ano$GA>=37],1,mean)
  esetS=apply(myb[,ano$GA>=37],1,sd)
  eset2=apply(myb,2,function(x){(x-esetM)/esetS})
  #eset2=apply(myb,2,function(x){(x-esetM)})
  eset2=myb
  
  
  SYMBOLS<-aT1$SYMBOL
  
  toad=NULL
  names_celsig<-c()
  for (cs in names(celsig)){
    
   if (length(which(SYMBOLS%in%celsig[[cs]]))>1)
   {
     toad=rbind(toad,apply(as.matrix(eset2[SYMBOLS%in%celsig[[cs]],]),2,mean))
     names_celsig<-c(names_celsig,cs)
   }
  
   }
  rownames(toad)<-names_celsig #names(celsig)
  
  tg=data.frame(ID=rownames(toad),ID2=rownames(toad))
  
  f3=function(xx,tg,nm="_Sigs.pdf"){
    
    #pdf(paste("./",nm,sep=""))
    pdf(paste(outFolder,nm,sep=""))
    par(mfrow=c(3,3))
    for(an in tg$ID){
      ano2=ano
      ano2$Y=(xx[an,])
      
      
      
      mod1=lmer(Y~x+x2+x3+(1|Main_Index),data = ano2,control=lmerControl(optimizer="bobyqa"),REML=FALSE)  
      mod2=lmer(Y~1+(1|Main_Index),data = ano2,control=lmerControl(optimizer="bobyqa"),REML=FALSE)  
      
      p=anova(mod1,mod2)$"Pr(>Chisq)"[2]
      pred=expand.grid(GA=seq(10,40,by=0.1),Main_Index="newpoint")
      pred$x=pred$GA/40
      pred$x2=pred$x^2
      pred$x3=pred$x^3
      pred$Y=predict(mod1,pred,allow.new.levels=TRUE)
      FC=round(2^(max(pred$Y)-min(pred$Y)),1)
      
      #white 
      par(bg = "white")
      plot(0,0,xlim=c(9,40),cex.lab=1.2,cex.axis=1.5,ylim=range(ano2$Y),
           main=paste(an," p=",round(p,3),"Fold Change=",FC),ylab="Average log2 Intensity",
           xlab="Gestational age (weeks)",cex.main=0.7)  
      for(s in unique(ano2$Main_Index)){
        tmp=ano2[ano2$Main_Index==s,]
        tmp=tmp[order(tmp$GA),]
        points(Y~GA,tmp,type="l",col="darkgrey",lty=1)
        
        points(Y~GA,tmp,col="darkgrey")
        
      }
      points(Y~GA,data=pred[pred$GA<=40&pred$GA>=10,],lwd=3,type="l",col="blue",lty=1)
      
    }
    dev.off()
    
  }
  
  
  #f3(xx=toad,tg=tg[c(13,11,7),],nm="Sigs_sel.pdf")
  f3(xx=toad,tg=tg,nm="Sigs_sel.pdf")
  
}


sig_cells_up<-signature_one_plot(SCGeneSets=res_up,outFolder="./13_Longitudinal_study/up_regulated_sig/")
sig_cells_down<-signature_one_plot(SCGeneSets=res_down,outFolder="./13_Longitudinal_study/down_regulated_sig/")
sig_cells_combined<-signature_one_plot(SCGeneSets=res_combined,outFolder="./13_Longitudinal_study/combined_sig/")
sig_cells_markers<-signature_one_plot(SCGeneSets=SCGeneSet,outFolder="./13_Longitudinal_study/")
#############################
# one plot

plot(1,1)
#legend("topleft",legend=sig_cells_up, fill=cluster.Colors[sig_cells_up],bty = "n",title="",cex=1)
#legend("topleft",legend=sig_cells_down, fill=cluster.Colors[sig_cells_down],bty = "n",title="",cex=1)
sig_cells_markers<-sig_cells_markers[which(!sig_cells_markers %in% c("0_Stromal-1","2_Macrophage-1","6_Decidual","8_LED","14_Macrophage-3" ,"20_Smooth muscle cells-3","12_Smooth muscle cells-1","13_Myofibroblast"))]
sig_cells_up<-sig_cells_up[which(!sig_cells_up %in% c("12_Smooth muscle cells-1"))]
sig_cells_down<-sig_cells_down[which(!sig_cells_down %in% c("3_Endothelial-1","6_Decidual","8_LED" ,"14_Macrophage-3"))]

legend("topleft",legend=sig_cells_up, fill=cluster.Colors[sig_cells_up],bty = "n",title="",cex=1)

signature_one_plot<-function(SCGeneSets,outFolder="./13_Longitudinal_study/up_regulated_sig/")
{

  system(paste0("mkdir -p ", outFolder))
  pile=NULL
  for(a in unique(SCGeneSets$cluster)){
    tmp=SCGeneSets[SCGeneSets$cluster==a,]
    tmp=tmp[order(tmp$avg_logFC,decreasing=TRUE),]
    pile=rbind(pile,tmp[1:min(c(20,dim(tmp)[1])),c("gene","cluster")])
  }
  names(pile)<-c("Symbol","Type")
  genes=pile
  #genes=rbind(genes,data.frame(Symbol=genes[genes$Type%in%c("Tcells-activated","Tcells-resting"),"Symbol"],Type="Tcell"))
  #genes=rbind(genes,data.frame(Symbol=genes[genes$Type%in%c("5_CD4_T-cell" ,"7_CD8_T-cell"),"Symbol"],Type="Tcell"))
  
  
  celsig=list()
  for(cs in unique(genes$Type)){
    celsig[[cs]]<-intersect(as.character(genes$Symbol[genes$Type==cs]),aT1$SYMBOL)
  }
  
  myb=b[rownames(aT1),]
  
  esetM=apply(myb[,ano$GA>=37],1,mean)
  esetS=apply(myb[,ano$GA>=37],1,sd)
  eset2=apply(myb,2,function(x){(x-esetM)/esetS})
  #eset2=apply(myb,2,function(x){(x-esetM)})
  eset2=myb
  
  
  SYMBOLS<-aT1$SYMBOL
  
  toad=NULL
  names_celsig<-c()
  for (cs in names(celsig)){
    
    if (length(which(SYMBOLS%in%celsig[[cs]]))>1)
    {
      toad=rbind(toad,apply(as.matrix(eset2[SYMBOLS%in%celsig[[cs]],]),2,mean))
      names_celsig<-c(names_celsig,cs)
    }
    
  }
  rownames(toad)<-names_celsig #names(celsig)
  
  tg=data.frame(ID=rownames(toad),ID2=rownames(toad))
  
  f3=function(xx,tg,nm="_Sigs.pdf"){
    
    pdf(paste(outFolder,nm,sep=""),width=4,height=5)
    
    flag=0
    #pdf(paste("./",nm,sep=""))
    
    #par(mar=c(4, 14 ,4.1 ,2.1))
    
    #par(mar=c(5, 4, 4, 8))
    
    sig_celltypes<-c()
    
    for(an in tg$ID){
      ano2=ano
      ano2$Y=(xx[an,])
      
      
      
      mod1=lmer(Y~x+x2+x3+(1|Main_Index),data = ano2,control=lmerControl(optimizer="bobyqa"),REML=FALSE)  
      mod2=lmer(Y~1+(1|Main_Index),data = ano2,control=lmerControl(optimizer="bobyqa"),REML=FALSE)  
      
      p=anova(mod1,mod2)$"Pr(>Chisq)"[2]
      if(p<=0.05)
    
      sig_celltypes<-c(sig_celltypes,an) 
    
      }
    #par(mfrow=c(3,3))
    #sig_celltypes<-c()
    for(an in tg$ID){
      ano2=ano
      ano2$Y=(xx[an,])
      
      
      
      mod1=lmer(Y~x+x2+x3+(1|Main_Index),data = ano2,control=lmerControl(optimizer="bobyqa"),REML=FALSE)  
      mod2=lmer(Y~1+(1|Main_Index),data = ano2,control=lmerControl(optimizer="bobyqa"),REML=FALSE)  
      
      p=anova(mod1,mod2)$"Pr(>Chisq)"[2]
      pred=expand.grid(GA=seq(10,40,by=0.1),Main_Index="newpoint")
      pred$x=pred$GA/40
      pred$x2=pred$x^2
      pred$x3=pred$x^3
      pred$Y=predict(mod1,pred,allow.new.levels=TRUE)
      FC=round(2^(max(pred$Y)-min(pred$Y)),1)
      
     
      #white 
      par(bg = "white")
      if(p<=0.05)
      {
        if (flag==0)
        {
          
          plot(0,0,xlim=c(9,40),cex.lab=1.2,cex.axis=1.5,ylim=range(xx[sig_celltypes,]),
               main="",ylab="Average log2 Intensity",
               xlab="Gestational age (weeks)",cex.main=0.7)#,labels=FALSE) #labels=FALSE #main=paste(an," p=",round(p,3),"Fold Change=",FC)
          
          #axis(side = 2, at = seq(0, range(xx)[2], by = 0.1))
          #axis(side = 1, at = seq(0, 40, by = 5))
          flag=1
          }
          
        
        #points(Y~GA,data=pred[pred$GA<=40&pred$GA>=10,],lwd=3,type="l",col="blue",lty=1)
        points(Y~GA,data=pred[pred$GA<=40&pred$GA>=10,],lwd=3,type="l",col=cluster.Colors[an],lty=1)
       
       
      }
    }
    
    
    #legend(x=5, y=3,legend=sig_celltypes, fill=cluster.Colors[sig_celltypes],bty = "n",title="",cex=1)
    #par(mar=c(5, 4, 4, 2) + 0.1)
    
    dev.off()
    
  }
  
  
  #f3(xx=toad,tg=tg[c(13,11,7),],nm="Sigs_sel.pdf")
  f3(xx=toad,tg=tg,nm="all_significant_sigs.pdf")
  return(names_celsig)
}

