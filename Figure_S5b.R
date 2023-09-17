#rm(list = ls())
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(hdf5r)
#### input data
load(file = 'sce1.1_united.RData')

##output data file
sam.name <- "Major2"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}

features=c("VEGFA","HIF1A","SPP1","KDR","ESM1")
#提取特定celltype亚群
sce_Myeloid<-subset(sce,celltype %in% c("Myeloid"))
sce_Endo<-subset(sce,celltype %in% c("Endothelial"))
sce_Epi<-subset(sce,celltype %in% c("Epithelial"))

#Myeloid
sce_Myeloid= SetIdent(sce_Myeloid,value="sample") 
heatmap_AveE <- AverageExpression(sce_Myeloid, assays = "RNA", features = features,verbose = TRUE) %>% .$RNA
#输出Myeloid heatmap_AveE数据
write.table(heatmap_AveE,file=paste0("./",sam.name,"/","sce_Myeloid_heatmap_AveE.txt"),sep="\t",quote = F)

#Endo
sce_Endo= SetIdent(sce_Endo,value="sample") 
heatmap_AveE <- AverageExpression(sce_Endo, assays = "RNA", features = features,verbose = TRUE) %>% .$RNA
#输出Endo heatmap_AveE数据
write.table(heatmap_AveE,file=paste0("./",sam.name,"/","sce_Endo_heatmap_AveE.txt"),sep="\t",quote = F)

#Epi
sce_Epi= SetIdent(sce_Epi,value="sample") 
heatmap_AveE <- AverageExpression(sce_Epi, assays = "RNA", features = features,verbose = TRUE) %>% .$RNA
#输出Epi heatmap_AveE数据
write.table(heatmap_AveE,file=paste0("./",sam.name,"/","sce_Epi_heatmap_AveE.txt"),sep="\t",quote = F)
#Endothelial_NormalvsTumor基因的差异
sce= SetIdent(sce,value="celltype") 
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
features <- c(
              "VEGFA","VEGFB","VEGFC","PGF",
              "PDGFA","PDGFB","ANGPT2",
              "APLN","APLNR",
              "FLT1","KDR","FLT4"#VEGFR1/2/3
)
pdf(paste0("./",sam.name,"/Endo_celltype_GF2",".pdf"),width =3.5,height = 3)
DotPlot(sce,group.by = "celltype", features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+RotatedAxis()
dev.off()
##celltype UMAP图
pdf(paste0("./",sam.name,"/celltype-DimPlot_umap_2",max(dim.use),".pdf"),width = 5,height = 4)
DimPlot(sce, reduction = "umap",group.by = 'celltype',label = T,raster=FALSE)
dev.off()
