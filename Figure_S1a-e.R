library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(hdf5r)
#### input data
load(file = 'sce_Endo.RData')
##output data file
sam.name <- "5_cohort_Endo3._2.0"
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
pdf(paste0("./",sam.name,"/sce_Endo_umap_group_tissue",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(sce_Endo, group.by = "tissue",reduction = 'umap',label=T)
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
              "ACTA2","CXCL9","CXCL10","CXCL11",
              "EPCAM","CD3D",#VEGFR1/2/3
              "MKI67","CD79A","LYZ","RGCC","FBLN2")
pdf(paste0("./",sam.name,"/niche_marker_DotPlot_CRC1_Endo_marker.pdf"),width =7,height = 7 )
DotPlot(sce_Endo,group.by = "seurat_clusters",features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+RotatedAxis()
dev.off()

table(sce_Endo@meta.data$seurat_clusters,sce_Endo@meta.data$tissue)
table(sce_Endo@meta.data$sample,sce_Endo@meta.data$tissue)
table(sce_Endo@meta.data$sample,sce_Endo@meta.data$seurat_clusters)
write.table(table(sce_Endo@meta.data$seurat_clusters,sce_Endo@meta.data$tissue),
            file=paste0("./",sam.name,"/sce_Endo_clusters_tissue.txt"),sep="\t",quote = F)
write.table(table(sce_Endo@meta.data$sample,sce_Endo@meta.data$tissue),
            file=paste0("./",sam.name,"/sce_Endo_sample_tissue.txt"),sep="\t",quote = F)
write.table(table(sce_Endo@meta.data$sample,sce_Endo@meta.data$seurat_clusters),
            file=paste0("./",sam.name,"/sce_Endo_sample_clusters.txt"),sep="\t",quote = F)


cluster5celltype <-c("0"="Vein",
                     "1"="Tip",
                     "2"="Vein",
                     "3"="Capillary",
                     "4"="Tip",
                     "5"="Artery",
                     "6"="Capillary",
                     "7"="Lymphatic",
                     "8"="Vein",
                     "9"="Vein",
                     "10"="Vein",
                     "11"="Tip",
                     "12"="Other_EC",
                     "13"="Capillary",
                     "14"="Mural",
                     "15"="Other_EC",
                     "16"="Mural",
                     "17"="Other_EC",
                     "18"="Vein",
                     "19"="Other_EC",
                     "20"="Other_EC",
                     "21"="Other_EC")
sce_Endo[['subcelltype2']] = unname(cluster5celltype[sce_Endo@meta.data$seurat_clusters])


##sce_Endo umap图(group.by = "tissue")
pdf(paste0("./",sam.name,"/sce_Endo_umap_group_subcelltype2",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(sce_Endo, group.by = "subcelltype2",reduction = 'umap',label=T)
dev.off()

#提取endo特定亚群
sce_VC<-subset(sce_Endo,subcelltype2 %in% c("Artery","Capillary","Tip","Vein","Lymphatic","Mural"))

## Mast_cell亚群重新聚类分组
sce_VC <- NormalizeData(sce_VC, normalization.method = "LogNormalize", scale.factor = 10000) 
sce_VC <- FindVariableFeatures(sce_VC, selection.method = 'vst', nfeatures = 2000)
sce_VC<- ScaleData(sce_VC, vars.to.regress = "percent.mt")
sce_VC<- RunPCA(sce_VC, features = VariableFeatures(object = sce_VC)) 

sce_VC <- FindNeighbors(sce_VC, dims = 1:50)
sce_VC <- FindClusters(sce_VC, resolution = 0.8 )
sce_VC <- RunUMAP(sce_VC, dims = 1:50)

##sce_VC umap图
sce_VC= SetIdent(sce_VC,value="subcelltype2") 
pdf(paste0("./",sam.name,"/sce_VC_umap_subcelltype2.pdf"),width = 5.5,height = 4)
DimPlot(sce_VC, reduction = 'umap',label=F)
dev.off()
##sce_VC umap图
sce_VC= SetIdent(sce_VC,value="seurat_clusters") 
pdf(paste0("./",sam.name,"/sce_VC_umap_seurat_clusters.pdf"),width = 5,height = 4)
DimPlot(sce_VC, reduction = 'umap',label=F)
dev.off()
##sce_VC umap图
sce_VC= SetIdent(sce_VC,value="tissue") 
pdf(paste0("./",sam.name,"/sce_VC_umap_tissue.pdf"),width = 5.5,height = 4)
DimPlot(sce_VC, reduction = 'umap',label=F)
dev.off()
##sce_VC umap图
sce_VC= SetIdent(sce_VC,value="ECtype") 
pdf(paste0("./",sam.name,"/sce_VC_umap_ECtype.pdf"),width = 15,height = 14)
DimPlot(sce_VC, reduction = 'umap',label=T)
dev.off()

##sce_VC umap图
sce_Endo= SetIdent(sce_Endo,value="seurat_clusters") 
cluster5celltype <-c("0"="EC08_IGFBP5",
                     "1"="EC05_VWA1",
                     "2"="EC11_ACKR1",
                     "3"="EC03_CD36",
                     "4"="EC06_ESM1",
                     "5"="EC01_GJA5",
                     "6"="EC02_BTNL9",
                     "7"="EC14_LYVE1",
                     "8"="EC12_MADCAM1",
                     "9"="EC10_ADIRF",
                     "10"="EC09_SELP",
                     "11"="EC07_PGF",
                     "12"="EC19_CD3D",
                     "13"="EC04_CA4",
                     "14"="EC15_ACTA2",
                     "15"="EC17_MKI67",
                     "16"="EC16_RGS5",
                     "17"="EC20_LYZ",
                     "18"="EC13_CPE",
                     "19"="EC21_MZB1",
                     "20"="EC18_CXCL9",
                     "21"="EC22_MS4A1")
sce_Endo[['ECtype']] = unname(cluster5celltype[sce_Endo@meta.data$seurat_clusters])
table(sce_Endo@meta.data$ECtype,sce_Endo@meta.data$tissue)
##major_marker_gene点状热图 
library(RColorBrewer)
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
features <- c("GJA5","GJA4","FBLN5",#"SERPINE2","TSPAN2",
              "BTNL9",#"CD320",
              "CD36","FABP5","CA4",
              #"PLVAP",#"RPLP1","RPS19",
             #"VWA1",
              "ESM1","NID2","PGF","VCAN",
              "SELP","ADIRF","ACKR1","MADCAM1","CPE",
              "PROX1","CCL21",
              "LYVE1",#"TNC","PRRX1",
              "RGS5",
              "ACTA2","MKI67","CXCL9",#"CXCL10","CXCL11",
              "CD3D",#VEGFR1/2/3
              "LYZ","MZB1","MS4A1")
pdf(paste0("./",sam.name,"/ECtype_DotPlot_CRC1_Endo_marker.pdf"),width =7,height = 7 )
DotPlot(sce_Endo,group.by = "ECtype",features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+RotatedAxis()
dev.off()

#提取endo特定亚群
sce_VC<-subset(sce_Endo,subcelltype2 %in% c("Artery","Capillary","Tip","Vein","Lymphatic","Mural"))

##major_marker_gene点状热图 
library(RColorBrewer)
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
features <- c("GJA5","GJA4","FBLN5",#"SERPINE2","TSPAN2",
              "BTNL9",#"CD320",
              "CD36","FABP5","CA4",
              #"PLVAP",#"RPLP1","RPS19",
              #"VWA1",
              "ESM1","NID2","PGF","VCAN",
              "SELP","ADIRF","ACKR1","MADCAM1","CPE",
              "PROX1","CCL21",
              "LYVE1",#"TNC","PRRX1",
              "RGS5",
              "ACTA2")
pdf(paste0("./",sam.name,"/ECtype_DotPlot_sce_VC_marker.pdf"),width =5.8,height = 5.2 )
DotPlot(sce_VC,group.by = "ECtype",features = features)+coord_flip()+
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
pdf(paste0("./",sam.name,"/HLA_DotPlot_CRC1_MHC_BVEC.pdf"),width =6,height = 4.5)
DotPlot(sce_Endo,group.by = "subcelltype2",features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+RotatedAxis()
dev.off()

##major_marker_gene_umap
pdf(paste0("./",sam.name,"/sce_Endo_RGCC_umap.pdf"),width = 5,height = 4)
FeaturePlot(sce_Endo, features = "RGCC",
            cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,
            reduction = "umap")
dev.off()

##major_marker_gene_umap
pdf(paste0("./",sam.name,"/sce_Endo_FBLN2_umap..pdf"),width = 5,height = 4)
FeaturePlot(sce_Endo, features = "FBLN2",
            cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,
            reduction = "umap")
dev.off()

##major_marker_gene_umap
pdf(paste0("./",sam.name,"/sce_VC_RGCC_umap.pdf"),width = 5,height = 4)
FeaturePlot(sce_VC, features = "RGCC",
            cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,
            reduction = "umap")
dev.off()

##major_marker_gene_umap
pdf(paste0("./",sam.name,"/sce_VC_FBLN2_umap..pdf"),width = 5,height = 4)
FeaturePlot(sce_VC, features = "FBLN2",
            cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,
            reduction = "umap")
dev.off()

##major_marker_gene_umap
pdf(paste0("./",sam.name,"/sce_VC_VEGFA_umap..pdf"),width = 5,height = 4)
FeaturePlot(sce_VC, features = "VEGFA",
            cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,
            reduction = "umap")
dev.off()



#Endothelial_NormalvsTumor基因的差异
sce_VC= SetIdent(sce_VC,value="subcelltype2") 
markers <- FindMarkers(sce_VC,  ident.1="Capillary", ident.2="Tip",
                       assay = 'RNA',slot = 'counts',
                       logfc.threshold =0,min.pct = 0 )
write.table(markers,file=paste0("./",sam.name,"/","sce_VC_Capillary_vs_Tip_All.txt"),sep="\t",quote = F)

#### 2.1 计算mast cluster marker gene all
sce_VC= SetIdent(sce_VC,value="subcelltype2") 
all.markers <- FindAllMarkers(sce_VC, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.table(all.markers,file=paste0("./",sam.name,"/","sce_VC_subcelltype2_all",max(dim.use),"PC.txt"),sep="\t",quote = F)

##major_marker_gene点状热图 aginogenesis
library(RColorBrewer)
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
features <- c("VEGFA","VEGFB","VEGFC","PGF",
              "PDGFA","PDGFB","ANGPT2",
              "APLN","APLNR",
              "FLT1","KDR","FLT4"#VEGFR1/2/3
)
pdf(paste0("./",sam.name,"/niche_marker_DotPlot_CRC1_VEGFR_BVEC.pdf"),width =3.7,height = 3)
DotPlot(sce_VC,group.by = "subcelltype2",features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+RotatedAxis()
dev.off()

#ESM1
features <- c("ESM1","CA4","PGF"
)
pdf(paste0("./",sam.name,"/ESM1.pdf"),width =3.7,height = 3)
DotPlot(sce_VC,group.by = "subcelltype2",features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+RotatedAxis()
dev.off()



#提取endo特定亚群sce_VBEC
sce_BVEC<-subset(sce_Endo,subcelltype2 %in% c("Artery","Capillary","Tip","Vein"))
##major_marker_gene_umap
pdf(paste0("./",sam.name,"/VBEC_RGCC_VlnPlot..pdf"),width = 5,height = 4)
VlnPlot(sce_VBEC,group.by = "subcelltype2", features = "RGCC",pt.size = 0
)
dev.off()

##major_marker_gene_umap
pdf(paste0("./",sam.name,"/sce_BVEC_FBLN2_VlnPlot..pdf"),width = 5,height = 4)
VlnPlot(sce_VBEC,group.by = "subcelltype2", features = "FBLN2",pt.size = 0
)
dev.off()

##major_marker_gene_umap
pdf(paste0("./",sam.name,"/VC_RGCC_VlnPlot..pdf"),width = 5,height = 4)
VlnPlot(sce_VC,group.by = "subcelltype2", features = "RGCC",pt.size = 0
)
dev.off()

##major_marker_gene_umap
pdf(paste0("./",sam.name,"/sce_BC_FBLN2_VlnPlot..pdf"),width = 5,height = 4)
VlnPlot(sce_VC,group.by = "subcelltype2", features = "FBLN2",pt.size = 0
)
dev.off()


##major_marker_gene点状热图
library(RColorBrewer)
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
features <- c("HLA-A","HLA-B","HLA-C","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DQB1",
              "HLA-DQB2","HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-E","HLA-F","HLA-G")
pdf(paste0("./",sam.name,"/HLA_DotPlot_CRC1_MHC_BVEC.pdf"),width =3.5,height = 4.5)
DotPlot(sce_BVEC,group.by = "subcelltype2",features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+RotatedAxis()
dev.off()

