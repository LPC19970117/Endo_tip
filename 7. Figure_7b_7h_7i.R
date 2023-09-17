library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(hdf5r)
#### input data
load(file = 'Endothelial.RData')

##output data file
sam.name <- "Endo_01"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}

## Mast_cell亚群重新聚类分组
sce_Endo <- NormalizeData(sce_Endo, normalization.method = "LogNormalize", scale.factor = 10000) 
sce_Endo <- FindVariableFeatures(sce_Endo, selection.method = 'vst', nfeatures = 2000)
sce_Endo<- ScaleData(sce_Endo, vars.to.regress = "percent.mt")
sce_Endo<- RunPCA(sce_Endo, features = VariableFeatures(object = sce_Endo)) 
ElbowPlot(sce_Endo, reduction = "pca",ndims = 50)

dim.use<-1:50
sce_Endo <- FindNeighbors(sce_Endo, dims = 1:50)
sce_Endo <- FindClusters(sce_Endo, resolution = 0.8 )
sce_Endo <- RunUMAP(sce_Endo, dims = 1:50)

##sce_Endo umap图
pdf(paste0("./",sam.name,"/sce_Endo_umap",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(sce_Endo, reduction = 'umap',label=T)
dev.off()
##sce_Endo umap图(split.by = "tissue")
pdf(paste0("./",sam.name,"/sce_Endo_umap_split_tissue",max(dim.use),"PC.pdf"),width = 10,height = 4)
DimPlot(sce_Endo, split.by = "tissue",reduction = 'umap',label=T)
dev.off()
##sce_Endo umap图(group.by = "tissue")
pdf(paste0("./",sam.name,"/sce_Endo_umap_group",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(sce_Endo, group.by = "group",reduction = 'umap',label=T)
dev.off()
##sce_Endo umap图(group.by = "tissue")
pdf(paste0("./",sam.name,"/sce_Endo_umap_tissue",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(sce_Endo, group.by = "tissue",reduction = 'umap',label=T)
dev.off()
#细胞类群相似树
sce <- BuildClusterTree(
  sce_Endo,
  dims = dim.use,
  reorder = F,
  reorder.numeric = F)
pdf(paste0("./",sam.name,"/Endo-ClusterTree_",max(dim.use),"PC.pdf"),width = 10,height = 8)
PlotClusterTree(sce)
dev.off()

##marker_gene_umap
pdf(paste0("./",sam.name,"/umap_marker_sce_Endo.pdf"),width = 15,height = 15)
FeaturePlot(sce_Endo, features = c("RGS5","ACTA2",#mural
                                 "GJA4","GJA5","FBLN5","TSPAN2",#artery
                                 "CA4","CD36","BTNL9",#capillary
                                 "ACKR1","SELP","CPE",#Vein
                                 "LYVE1","CCL21","PROX1",#lymphatic
                                 "ESM1","PGF","COL4A1",#Tip-like
                                 "MKI67","CD3D","CD79A","EPCAM1"
),
cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,
reduction = "umap")
dev.off()

##major_marker_gene点状热图 aginogenesis
library(RColorBrewer)
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
features <- c("GJA5","GJA4","FBLN5","SERPINE2","TSPAN2",
              "BTNL9","CD320","CD36","CA4",
              "FABP5","PLVAP","RPLP1","RPS19",
              "ESM1","PGF","NID2","VCAN",
              "SELP","ACKR1","CPE",
              "PROX1","CCL21",
              "LYVE1","TNC","PRRX1","RGS5",
              "ACTA2","CD74","HLA-DRA",
              "EPCAM","CD3D",#VEGFR1/2/3
              "MKI67","CD79A","LYZ","RGCC","FBLN2")
pdf(paste0("./",sam.name,"/niche_marker_DotPlot_CRC1_Endo_marker.pdf"),width =7,height = 7 )
DotPlot(sce_Endo,group.by = "seurat_clusters",features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+RotatedAxis()
dev.off()

#### 2.1 计算mast cluster marker gene all
#sce_Endo= SetIdent(sce_Endo,value="subcelltype2") 
all.markers <- FindAllMarkers(sce_Endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/","sce_Endo",max(dim.use),"PC.txt"),sep="\t",quote = F)

#Endothelial_NormalvsTumor基因的差异
sce_Endo= SetIdent(sce_Endo,value="tissue") 
markers <- FindMarkers(sce_Endo,  ident.1="Normal", ident.2="Tumor",
                       assay = 'RNA',slot = 'counts',
                       logfc.threshold =0.25,min.pct = 0.25 )
write.table(markers,file=paste0("./",sam.name,"/","sce_Endo_NvsT.txt"),sep="\t",quote = F)
#Endothelial_NormalvsTumor基因的差异
sce_Endo= SetIdent(sce_Endo,value="group") 
markers <- FindMarkers(sce_Endo,  ident.1="A", ident.2="B",
                       assay = 'RNA',slot = 'counts',
                       logfc.threshold =0,min.pct = 0 )
write.table(markers,file=paste0("./",sam.name,"/","sce_Endo_AvsB_all.txt"),sep="\t",quote = F)


##major_marker_gene点状热图 aginogenesis
library(RColorBrewer)
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
features <- c("VEGFA","VEGFB","VEGFC","VEGFD","VEGFE","PGF",
              "HIF1A","EGF","HGF",
              "PDGFA","PDGFB","PDGFC","PDGFD",
              "ANGPT1","ANGPT2","ANGPT4","TNF",
              "FGF2","FGF12","TGFB1",
              "APLN","APLNR",
              "CD34","PROM1","CD274","PTN",
              "FLT1","KDR","FLT4",#VEGFR1/2/3
              "FGFR1","FGFR2","FGFR3","FGFR4","TIE1","TEK",#="TIE2",
              "PDGFRA","PDGFRB","EPHA2","EPHB1","EPHB4",
              
              "NOTCH1","ID1","ID2","ID3","HES1","JAG1","DLL4","COL4A1","COL4A2","LAMB1","LAMA4",
              "FBLN2","PXDN","CD93",
              "CXCR2","CXCR4","CXCL12","CCL2",
              "MMP2","MMP9","MMP15")
pdf(paste0("./",sam.name,"/Endo_vegf.pdf"),width =3.5,height = 8)
DotPlot(sce_Endo,group.by = "group",features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+RotatedAxis()
dev.off()

##major_marker_gene点状热图
library(RColorBrewer)
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
features <- c("HLA-A","HLA-B","HLA-C","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DQB1",
              "HLA-DQB2","HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-E","HLA-F","HLA-G")
pdf(paste0("./",sam.name,"/niche_marker_DotPlot_CRC1_MHC_Endo.pdf"),width =4,height = 4.5)
DotPlot(sce_Endo,group.by = "group",features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+RotatedAxis()
dev.off()

##major_marker_gene点状热图
library(RColorBrewer)
#Endothelial_NormalvsTumor基因的差异
sce_Endo= SetIdent(sce_Endo,value="group") 
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
features <- c("VEGFA","VEGFB","VEGFC","PGF",
              "PDGFA","PDGFB","ANGPT2",
              "APLN","APLNR",
              "FLT1","KDR","FLT4"#VEGFR1/2/3
)
pdf(paste0("./",sam.name,"/niche_marker_DotPlot_CRC2_2",".pdf"),width =4,height = 3)
DotPlot(sce_Endo, features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+RotatedAxis()
dev.off()
##major_marker_gene_umap
pdf(paste0("./",sam.name,"/sce_VC_ESM1.pdf"),width = 5,height = 4)
FeaturePlot(sce_Endo, features = "ESM1",
            cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,
            reduction = "umap")
dev.off()
##major_marker_gene_umap
pdf(paste0("./",sam.name,"/sce_ESM1_VlnPlot..pdf"),width = 5,height = 4)
VlnPlot(sce_Endo,group.by = "group", features = "ESM1",pt.size = 0
)
dev.off()

#4.提取Endothelial亚组
sce_EndoABC<-subset(sce_Endo,group %in% c("A","B","C"))
table(sce_EndoABC@meta.data$celltype)
#4.存储Endothelial数据
save(sce_EndoABC,file=paste0("./",sam.name,"/","EndoABC.RData"))



##major_marker_gene点状热图
library(RColorBrewer)
#Endothelial_NormalvsTumor基因的差异
sce_EndoABC= SetIdent(sce_EndoABC,value="group") 
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
features <- c("ESM1","NID2","COL4A1","COL4A2","VWA1","LAMA4","MCAM",
              "VEGFA","VEGFB","VEGFC","PGF",
              "PDGFA","PDGFB","ANGPT2",
              "APLN","APLNR",
              "FLT1","KDR","FLT4"#VEGFR1/2/3
)
pdf(paste0("./",sam.name,"/EndoABC_DotPlot_GF",".pdf"),width =3.1,height = 3.8)
DotPlot(sce_EndoABC,group.by = "group", features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+RotatedAxis()
dev.off()

##major_marker_gene点状热图
library(RColorBrewer)
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
features <- c("HLA-A","HLA-B","HLA-C","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DQB1",
              "HLA-DQB2","HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-E","HLA-F","HLA-G")
pdf(paste0("./",sam.name,"/EndoABC_MHC_Endo.pdf"),width =3.3,height = 4.2)
DotPlot(sce_EndoABC, group.by = "group",features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+RotatedAxis()
dev.off()

