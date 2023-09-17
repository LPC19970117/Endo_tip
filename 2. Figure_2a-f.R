library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(hdf5r)
#### input data
load(file = 'Endothelial.RData')

##output data file
sam.name <- "GSE178341_Endo3_ratio_dot"
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

#细胞类群相似树
sce <- BuildClusterTree(
  sce_Endo,
  dims = dim.use,
  reorder = F,
  reorder.numeric = F)
pdf(paste0("./",sam.name,"/Endo-ClusterTree_",max(dim.use),"PC.pdf"),width = 10,height = 8)
PlotClusterTree(sce)
dev.off()

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

#### 2.1 计算mast cluster marker gene all
sce_Endo= SetIdent(sce_Endo,value="ECtype") 
all.markers <- FindAllMarkers(sce_Endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/","sce_Endo_ECtype.txt"),sep="\t",quote = F)

##Endo ECtype umap
pdf(paste0("./",sam.name,"/Endo_ECtype_umap_",max(dim.use),".pdf"),width = 8,height = 7)
DimPlot(sce_Endo, reduction = "umap",group.by = 'ECtype',label = T)
dev.off()

##Endo subcelltype2 umap
pdf(paste0("./",sam.name,"/Endo_subcelltype2_umap_",max(dim.use),".pdf"),width = 8,height = 7)
DimPlot(sce_Endo, reduction = "umap",group.by = 'subcelltype2',label = T)
dev.off()

#### 2.1 计算mast cluster marker gene all
sce_Endo= SetIdent(sce_Endo,value="subcelltype2") 
all.markers <- FindAllMarkers(sce_Endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/","sce_Endo_subcelltype2_all",max(dim.use),"PC.txt"),sep="\t",quote = F)

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

##sce_VC umap图
sce_VC= SetIdent(sce_VC,value="subcelltype2") 
pdf(paste0("./",sam.name,"/sce_VC_umap",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(sce_VC, reduction = 'umap',label=F)
dev.off()

#将sce_VC数据写到文件中一边后续分析
save(sce_VC,file=paste0("./",sam.name,"/","sce_VC.RData"))

##marker_gene_umap
pdf(paste0("./",sam.name,"/umap_marker_geneVEC.",max(dim.use),"PC.pdf"),width = 15,height = 15)
FeaturePlot(sce_VC, features = c("RGS5","ACTA2",#mural
                                  "GJA4","GJA5","FN1","TSPAN2",#artery
                                  "CA4","RGCC","BTNL9",#capillary
                                  "ACKR1","PLVAP","NR2F2",#Vein
                                  "SELP",#STALK
                                  "LYVE1",#lymphatic
                                  "MKI67",#proliferation
                                  "MT-CO2",#MT
                                  "VWA1","VWF",#Vein
                                  "CXCR4","LXN","PGF","COL4A1"#Tip-like
),
cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,
reduction = "umap")
dev.off()

##marker_gene_umap
pdf(paste0("./",sam.name,"/umap_marker_gene_dotVEC.",max(dim.use),"PC.pdf"),width = 8,height = 7)
DotPlot(sce_VC, features = c("RGS5","ACTA2",#mural
                              "GJA4","GJA5","FN1","TSPAN2",#artery
                              "CA4","RGCC","BTNL9",#capillary
                              "ACKR1","PLVAP","NR2F2",#Vein
                              "LYVE1",#lymphatic
                              "MKI67",#proliferation
                              "MT-CO2",#MT
                              "VWA1","VWF",#Vein
                              "CXCR4","LXN","PGF","COL4A1"#Tip-like
))+RotatedAxis()
dev.off()

##marker_gene_umap
sce_Endo= SetIdent(sce_Endo,value="ECtype") 
pdf(paste0("./",sam.name,"/umap_marker_gene_dot2_ECtype.",max(dim.use),"PC.pdf"),width = 14,height = 6)
DotPlot(sce_Endo,group.by = "ECtype",
        features =
          c("GJA5","FBLN5", "TSPAN2", "SOX17", #Artery:
            "PGF", "ESM1", "NID2", "LXN", "CXCR4", "ADM", "COL4A1", "COL4A2",#Tip:
            "HSPG2", "VWA1", "APLNR", #Immature:
            "CA4","CD36","BTNL9","PLVAP","KDR", # "RGCC",#"CD24", #Capillary:
            "MFSD2A", "RGCC",  #Microvascular:
            "ACKR1","SELP","VWF","VCAM1", "DUSP23",  #Vein: 
            "PROX1", "PDPN", "LYVE1",    #Lymphatic:
            "RGS5", #Mural:
            "ACTA2", #vSMC:
            "MKI67",#Proliferation:
            "CD3D", #Doublet-T: 
            "CDH1" #MC_CDH1:
          ))+RotatedAxis()
dev.off()
#### 2.1 计算mast cluster marker gene all
sce_VC= SetIdent(sce_VC,value="subcelltype2") 
all.markers <- FindAllMarkers(sce_VC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/","sce_VC_subcelltype2",max(dim.use),"PC.txt"),sep="\t",quote = F)
#### 2.1 计算mast cluster marker gene all
sce_VC= SetIdent(sce_VC,value="subcelltype2") 
all.markers <- FindAllMarkers(sce_VC, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.table(all.markers,file=paste0("./",sam.name,"/","sce_VC_subcelltype2_all",max(dim.use),"PC.txt"),sep="\t",quote = F)

#### 2.1 计算mast cluster marker gene all
sce_VC= SetIdent(sce_VC,value="ECtype") 
all.markers <- FindAllMarkers(sce_VC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/","sce_VC_ECtype",max(dim.use),"PC.txt"),sep="\t",quote = F)

##marker_gene_umap
sce_VC= SetIdent(sce_VC,value="ECtype") 
pdf(paste0("./",sam.name,"/umap_marker_gene_dot2_ECtype.",max(dim.use),"PC.pdf"),width = 10,height = 5)
DotPlot(sce_VC,group.by = "ECtype",
        features =
          c("GJA5","GJA4","FBLN5", "SERPINE2","TSPAN2", #Artery:
            "BTNL9", "CD36", "CA4", "FABP5","PLVAP",#cap
            "RPLP1", "RPS19", #Immature:
            "COL4A1","COL4A2","ESM1","PGF", # TIP
             # "MFSD2A", "RGCC",  #Microvascular:
            "ACKR1","VCAN","SELP","VWF", "CPE",  #Vein: 
            "CCL21","PROX1", "LYVE1", #Lymphatic:
            "RGS5", #Mural:
            "ACTA2" #vSMC:
          ))+RotatedAxis()
dev.off()

##major_marker_gene点状热图
library(RColorBrewer)
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
features <- c("GJA5","GJA4","FBLN5", "SERPINE2","TSPAN2", #Artery:
              "BTNL9", "CD320","CD36","CA4","FABP5","PLVAP",#cap
              "RPLP1", "RPS19",#"PRCP","LGALS1", #Immature:
             # "COL4A1","COL4A2",
              "ESM1","PGF","NID2", # TIP
              #"RGCC",  #Microvascular:
              "VCAN","SELP","ACKR1", "CPE",  #Vein: 
              "PROX1","CCL21", "LYVE1", #Lymphatic:
              "TNC","PRRX1","RGS5", #Mural:
              "ACTA2"#vSMC:
)
pdf(paste0("./",sam.name,"/niche_marker_DotPlot_CRC1",".pdf"),width =5.5,height = 6.5)
DotPlot(sce_VC,group.by = "ECtype",features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+RotatedAxis()
dev.off()

##sce_VC umap图(split.by = "tissue")
pdf(paste0("./",sam.name,"/sce_VC_umap_split_tissue",max(dim.use),"PC.pdf"),width = 10,height = 4)
DimPlot(sce_VC, split.by = "tissue",reduction = 'umap',label=F)
dev.off()
##sce_VC umap图(group.by = "tissue")
pdf(paste0("./",sam.name,"/sce_VC_umap_group_tissue",max(dim.use),"PC.pdf"),width = 5,height = 4.5)
DimPlot(sce_VC, group.by = "tissue",reduction = 'umap',label=F)
dev.off()

##sce_VC umap图(group.by = "seurat_clusters")
pdf(paste0("./",sam.name,"/sce_VC_umap_seurat_clusters_tissue",max(dim.use),"PC.pdf"),width = 5,height = 4.5)
DimPlot(sce_VC, group.by = "seurat_clusters",reduction = 'umap',label=F)
dev.off()
##sce_VC umap图(group.by = "ECtype")
pdf(paste0("./",sam.name,"/sce_VC_umap_ECtype_tissue",max(dim.use),"PC.pdf"),width = 5,height = 4.5)
DimPlot(sce_VC, group.by = "ECtype",reduction = 'umap',label=F)
dev.off()
#细胞类群相似树
sce <- BuildClusterTree(
  sce_VC,
  dims = dim.use,
  reorder = F,
  reorder.numeric = F)
pdf(paste0("./",sam.name,"/Endo-ClusterTree_",max(dim.use),"PC.pdf"),width = 10,height = 8)
PlotClusterTree(sce)
dev.off()

#输出sample_tissue
write.csv(table(sce_VC@meta.data$ECtype,sce_VC@meta.data$tissue),file="ECtype_tissue_VC.csv")
table(sce_VC@meta.data$ECtype,sce_VC@meta.data$tissue)

#提取endo特定亚群
sce_Endo_Mural<-subset(sce_Endo,subcelltype2 %in% c("Mural"))
#sce_Endo_Mural_NormalvsTumor基因的差异
#Endothelial_NormalvsTumor基因的差异
sce_Endo_Mural= SetIdent(sce_Endo_Mural,value="tissue") 
markers <- FindMarkers(sce_Endo_Mural,  ident.1="N", ident.2="T",
                       assay = 'RNA',slot = 'counts',
                       logfc.threshold =0,min.pct = 0 )
write.table(markers,file=paste0("./",sam.name,"/","sce_Endo_Mural_NvsT_All.txt"),sep="\t",quote = F)



##major_marker_gene点状热图
library(RColorBrewer)
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
features <- c("HLA-A","HLA-B","HLA-C","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DQB1",
          "HLA-DQB2","HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-E","HLA-F","HLA-G")
pdf(paste0("./",sam.name,"/niche_marker_DotPlot_CRC1_MHC",".pdf"),width =4,height = 4)
DotPlot(sce_VC,group.by = "subcelltype2",features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+RotatedAxis()
dev.off()


sce_BVEC<-subset(sce_VC,subcelltype2 %in% c("Tip","Immature","Artery","Capillary","Vein"))

##major_marker_gene点状热图
library(RColorBrewer)
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
features <- c("HLA-A","HLA-B","HLA-C","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DQB1",
              "HLA-DQB2","HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-E","HLA-F","HLA-G")
pdf(paste0("./",sam.name,"/niche_marker_DotPlot_CRC1_MHC_BVEC.pdf"),width =4,height = 4.5)
DotPlot(sce_BVEC,group.by = "subcelltype2",features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+RotatedAxis()
dev.off()

##major_marker_gene_umap
pdf(paste0("./",sam.name,"/sce_VC_RGCC_umap.pdf"),width = 5,height = 4)
FeaturePlot(sce_BVEC, features = "RGCC",
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
pdf(paste0("./",sam.name,"/sce_VC_RGCC_VlnPlot..pdf"),width = 5,height = 4)
VlnPlot(sce_VC,group.by = "subcelltype2", features = "RGCC",pt.size = 0
)
dev.off()

##major_marker_gene_umap
pdf(paste0("./",sam.name,"/sce_BVEC_FBLN2_VlnPlot..pdf"),width = 5,height = 4)
VlnPlot(sce_BVEC,group.by = "subcelltype2", features = "FBLN2",pt.size = 0
)
dev.off()

##major_marker_gene_umap
pdf(paste0("./",sam.name,"/sce_VC_maker_umap..pdf"),width =7,height = 9)
FeaturePlot(sce_VC, features = c("GJA5","CA4","ESM1","ACKR1","PROX1","RGS5"),
            cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,
            reduction = "umap")
dev.off()







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
pdf(paste0("./",sam.name,"/niche_marker_DotPlot_CRC1_VEGFR_BVEC.pdf"),width =3.5,height = 8)
DotPlot(sce_BVEC,group.by = "subcelltype2",features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+RotatedAxis()
dev.off()


#Endothelial_NormalvsTumor基因的差异
sce_VC= SetIdent(sce_VC,value="subcelltype2") 
markers <- FindMarkers(sce_VC,  ident.1="Capillary", ident.2="Tip",
                       assay = 'RNA',slot = 'counts',
                       logfc.threshold =0,min.pct = 0 )
write.table(markers,file=paste0("./",sam.name,"/","sce_VC_Capillary_vs_Tip_All.txt"),sep="\t",quote = F)

##输出样本临床和细胞数，用于后续比较分析
table(sce_Endo@meta.data$ECtype,sce_Endo@meta.data$tissue)
table(sce_Endo@meta.data$sample,sce_Endo@meta.data$ECtype)
table(sce_Endo@meta.data$sample,sce_Endo@meta.data$tissue)
write.csv(table(sce_Endo@meta.data$ECtype,sce_Endo@meta.data$tissue),file="sce_Endo_ECtype_tissue.csv")
write.csv(table(sce_Endo@meta.data$sample,sce_Endo@meta.data$ECtype),file="sce_Endo_sample_ECtype.csv")
write.csv(table(sce_Endo@meta.data$sample,sce_Endo@meta.data$tissue),file="sce_Endo_sample_tissue.csv")

write.table(table(sce_Endo@meta.data$seurat_clusters,sce_Endo@meta.data$tissue),
            file=paste0("./",sam.name,"/sce_Endo_clusters_tissue.txt"),sep="\t",quote = F)
write.table(table(sce_Endo@meta.data$sample,sce_Endo@meta.data$tissue),
            file=paste0("./",sam.name,"/sce_Endo_sample_tissue.txt"),sep="\t",quote = F)
write.table(table(sce_Endo@meta.data$sample,sce_Endo@meta.data$ECtype),
            file=paste0("./",sam.name,"/sce_Endo_sample_ECtype.txt"),sep="\t",quote = F)

