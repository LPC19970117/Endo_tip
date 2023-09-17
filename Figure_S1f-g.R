rm(list = ls())
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(hdf5r)
library(pheatmap)


MHCII <- readxl::read_xlsx("MHCII.xlsx")
#转换成list
sce_VC= SetIdent(sce_VC,value="subcelltype2") 
#4.提取指定单细胞亚群
gene <- as.list(MHCII)
sce_VC <- AddModuleScore(
  object = sce_VC,
  features = gene,
  ctrl = 100, #默认值是100
  name = 'MHCII'
)
pdf(paste0("./",sam.name,"/MHCII",".pdf"),width = 10,height = 7)
VlnPlot(sce_VC,features = 'MHCII1',pt.size = 0
)+geom_boxplot(width=0.2,col="black",fill="white")
dev.off()
write.table(sce_VC@meta.data,"sce_VC@meta.data_Angiogenesis",sep="\t",quote = F)




#4.提取指定单细胞亚群sce_BVEC
sce_BVEC<-subset(sce_VC,subcelltype2 %in% c("Tip","Immature","Artery","Capillary","Vein"))

gene <- as.list(MHCII)
sce_BVEC<- AddModuleScore(
  object = sce_BVEC,
  features = gene,
  ctrl = 100, #默认值是100
  name = 'MHCII'
)
pdf(paste0("./",sam.name,"/sce_BVEC_MHCII",".pdf"),width = 10,height = 7)
VlnPlot(sce_BVEC,group.by="subcelltype2",features = 'MHCII1',pt.size = 0
)+geom_boxplot(width=0.2,col="black",fill="white")
dev.off()
write.table(sce_BVEC@meta.data,"sce_BVEC@meta.data_MHCII",sep="\t",quote = F)

###MHCI
MHCI <- readxl::read_xlsx("MHCI.xlsx")
gene <- as.list(MHCI)
sce_BVEC<- AddModuleScore(
  object = sce_BVEC,
  features = gene,
  ctrl = 100, #默认值是100
  name = 'MHCI'
)
pdf(paste0("./",sam.name,"/sce_BVEC_MHCI",".pdf"),width = 10,height = 7)
VlnPlot(sce_BVEC,group.by="subcelltype2",features = 'MHCI1',pt.size = 0
)+geom_boxplot(width=0.2,col="black",fill="white")
dev.off()
write.table(sce_BVEC@meta.data,"sce_BVEC@meta.data_MHCI",sep="\t",quote = F)

