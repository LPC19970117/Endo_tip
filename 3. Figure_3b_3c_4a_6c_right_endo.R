rm(list = ls())
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(hdf5r)
#### input data
load(file = 'Endothelial.RData')
##output data file
sam.name <- "BVEC_热图"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}
## Endo_cell亚群重新聚类分组
sce_Endo <- NormalizeData(sce_Endo, normalization.method = "LogNormalize", scale.factor = 10000) 
sce_Endo <- FindVariableFeatures(sce_Endo, selection.method = 'vst', nfeatures = 2000)
sce_Endo<- ScaleData(sce_Endo, vars.to.regress = "percent.mt")
sce_Endo<- RunPCA(sce_Endo, features = VariableFeatures(object = sce_Endo)) 
ElbowPlot(sce_Endo, reduction = "pca",ndims = 50)

dim.use<-1:50
sce_Endo <- FindNeighbors(sce_Endo, dims = 1:50)
sce_Endo <- FindClusters(sce_Endo, resolution = 0.8 )
sce_Endo <- RunUMAP(sce_Endo, dims = 1:50)

cluster5celltype <-c("0"="Tip",
                     "1"="Vein",
                     "2"="Immature",
                     "3"="Capillary",
                     "4"="Vein",
                     "5"="Capillary",
                     "6"="Mural",
                     "7"="Artery",
                     "8"="Artery",
                     "9"="Other_EC",
                     "10"="Lymphatic",
                     "11"="Other_EC",
                     "12"="Capillary",
                     "13"="Other_EC",
                     "14"="Other_EC",
                     "15"="Mural",
                     "16"="Mural")
sce_Endo[['subcelltype2']] = unname(cluster5celltype[sce_Endo@meta.data$seurat_clusters])

cluster5celltype <-c("0"="EC07_ESM1",
                     "1"="EC09_ACKR1",
                     "2"="EC06_PRCP",
                     "3"="EC05_CA4",
                     "4"="EC08_VCAN",
                     "5"="EC04_CD36",
                     "6"="EC13_RGS5",
                     "7"="EC01_GJA5",
                     "8"="EC02_SERPINE2",
                     "9"="EC17_LowQ",
                     "10"="EC10_PROX1",
                     "11"="EC14_MKI67",
                     "12"="EC03_BTNL9",
                     "13"="EC15_CDH1",
                     "14"="EC16_CD3D",
                     "15"="EC12_PRRX1",
                     "16"="EC11_TNC")
sce_Endo[['ECtype']] = unname(cluster5celltype[sce_Endo@meta.data$seurat_clusters])
table(sce_Endo@meta.data$ECtype,sce_Endo@meta.data$tissue)
#提取endo特定亚群
sce_VC<-subset(sce_Endo,subcelltype2 %in% c("Tip","Immature","Artery","Capillary","Vein","Lymphatic","Mural"))

## Mast_cell亚群重新聚类分组
sce_VC <- NormalizeData(sce_VC, normalization.method = "LogNormalize", scale.factor = 10000) 
sce_VC <- FindVariableFeatures(sce_VC, selection.method = 'vst', nfeatures = 2000)
sce_VC<- ScaleData(sce_VC, vars.to.regress = "percent.mt")
sce_VC<- RunPCA(sce_VC, features = VariableFeatures(object = sce_VC)) 

sce_VC <- FindNeighbors(sce_VC, dims = 1:50)
sce_VC <- FindClusters(sce_VC, resolution = 0.8 )
sce_VC <- RunUMAP(sce_VC, dims = 1:50)


##marker_gene_umap  Figure 3c right
pdf(paste0("./",sam.name,"/umap_sce_VC_ESM1_01.pdf"),width = 8,height = 7)
FeaturePlot(sce_VC, features = c("ESM1"),
            cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,
            reduction = "umap",raster=FALSE)
dev.off()

##major_marker_gene点状热图 Figure 4a right
library(RColorBrewer)
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
features <- c("VEGFA","VEGFB","VEGFC","VEGFD","PGF",
              "PDGFA","PDGFB",#"PDGFC","PDGFD",
              "ANGPT2",#"FGF12","EPHA2","EPHB1","EPHB4",
              "APLN","APLNR",
              "FLT1","KDR","FLT4",#VEGFR1/2/3
              #"FGFR1",#"FGFR2","FGFR3","FGFR4",
              "TIE1","TEK"#="TIE2",
)
pdf(paste0("./",sam.name,"/DotPlot_CRC1_sce_VC.pdf"),width =3.8,height = 3.8)
DotPlot(sce_VC,group.by = "subcelltype2",
        features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+RotatedAxis()
dev.off()


##major_marker_gene点状热图 Figure 3b right
library(RColorBrewer)
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
features <- c("ESM1")
pdf(paste0("./",sam.name,"/DotPlot_sce_VC_ESM1.pdf"),width =3.5,height = 3)
DotPlot(sce_VC,group.by = "subcelltype2",features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+RotatedAxis()
dev.off()


##major_marker_gene点状热图 Figure 6c right
library(RColorBrewer)
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
features <- c("NRP1","NRP2")
pdf(paste0("./",sam.name,"/DotPlot_sce_VC_NRP1.pdf"),width =4,height = 1.2)
DotPlot(sce_VC,group.by = "subcelltype2",features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+RotatedAxis()
dev.off()
