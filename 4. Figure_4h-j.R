
##载入SpaceRanger数据
rm(list = ls())
#载入R包；
library(Seurat)
library(ggplot2)
library(patchwork)
library(data.table)
library(dplyr)
library(hdf5r)
# Error in library(hdf5r) : 不存在叫‘hdf5r’这个名字的程辑包
# 缺什么就安装什么
# BiocManager::install("hdf5r")
#library(hdf5r)

################################################################################
CRC <- Load10X_Spatial(data.dir ="./", #本地文件目录下运行
                       filename = "filtered_feature_bc_matrix.h5",#读取h5表达谱文件
                       slice ="CRC")

#这里文件夹的名字可以修改，但最好只用英文字母
sam.name <- "Endo"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}
# An object of class Seurat 
# 17943 features across 4371 samples within 1 assay 
# Active assay: Spatial (17943 features, 0 variable features)

p1 <- VlnPlot(CRC, 
              features ="nCount_Spatial",
              pt.size = 0,
              cols ="red")
pdf(paste0("./",sam.name,"/nCount.pdf"),width = 8,height = 4.5)
SpatialFeaturePlot(CRC, features ="nCount_Spatial")
dev.off()
p1 | p2

##################################################
## spatial天生就是Seurat对象
# 数据标准化，使用SCTransform方法进行标准化

CRC <- SCTransform(CRC, assay = "Spatial", verbose = FALSE)
CRC <- RunPCA(CRC, assay = "SCT", verbose = FALSE) 
plot1 <- DimPlot(CRC, reduction = "pca", group.by="orig.ident")
plot2 <- ElbowPlot(CRC, ndims=20, reduction="pca") 
plot1+plot2
#选取主成分
pc.num=1:15
## 细胞聚类
# 图聚类分群
CRC <- FindNeighbors(CRC, reduction = "pca", dims = pc.num)
CRC <- FindClusters(CRC, verbose = FALSE)
# UMAP降维可视化
CRC <- RunUMAP(CRC, reduction = "pca", dims =  pc.num)
p1 <- DimPlot(CRC, reduction = "umap", label = TRUE)
# 使用SpatialDimPlot函数进行可视化
p2 <- SpatialDimPlot(CRC, label = TRUE, label.size = 3)
p1 + p2
pdf(paste0("./",sam.name,"/cluster.pdf"),width = 8,height = 4.5)
p1 + p2
dev.off()
save(CRC,file = "CRC.Rds")
View(CRC@meta.data)
################################################################################
pdf(paste0("./",sam.name,"/MC_5markers.pdf"),width = 15,height = 14)
SpatialFeaturePlot(CRC, features =c("VEGFA","VEGFB","VEGFC","VEGFD","KDR","ESM1","PGF","CA4","CD36","PECAM1","PECAM1","EPCAM","COL1A2"),
)
dev.off()


#以下重命名髓系亚群
#方法二：使用unname函数配合向量：
Region_type <-c("0"="Tumor",
                      "1"="Tumor",
                      "2"="Stromal",
                      "3"="Tumor",
                      "4"="Normal",
                      "5"="Normal",
                      "6"="Normal",
                      "7"="Tumor",
                      "8"="Tumor",
                      "9"="Normal"
                      
)
CRC[['Region']] = unname(Region_type[CRC@meta.data$seurat_clusters])
View(CRC@meta.data) 

pdf(paste0("./",sam.name,"/Region_type_",".pdf"),width = 8,height = 7)
DimPlot(CRC, reduction = "umap",group.by = 'Region',label = T)
dev.off()

pdf(paste0("./",sam.name,"/niche_region_type_ST",".pdf"),width = 8,height = 7)
SpatialDimPlot(CRC, group.by = 'Region',label = TRUE, label.size = 3)
dev.off()

CRC= SetIdent(CRC,value="Region") 
pdf(paste0("./",sam.name,"/region_MC_5markers_Dot.pdf"),width = 6,height = 3)
DotPlot(CRC, features =c("VEGFA","KDR","ESM1","PGF","CA4","CD36","PECAM1","KIT","KITLG"))+ RotatedAxis()
dev.off()

features =c("VEGFA","VEGFB","VEGFC","VEGFD","PGF","FLT1","KDR","FLT4","ESM1","PECAM1","DCN","CD68","LYZ","CEACAM5","EPCAM")

#### 2.1 计算marker基因
CRC= SetIdent(CRC,value="Region") 
all.markers <- FindAllMarkers(CRC, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.table(all.markers,file=paste0("./",sam.name,"/Region","PC.txt"),sep="\t",quote = F)

# gene filter 
CRC_RM <- CRC[!rownames(CRC) %in% GenesRM,]
all.markers <- FindAllMarkers(CRC_RM , only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(all.markers,file=paste0("./",sam.name,"/CRC_RM_Region","PC.txt"),sep="\t",quote = F)

write.table(CRC@assays$Spatial@counts,file=paste0("./",sam.name,"/CRC@assays$Spatial@counts","PC.txt"),sep="\t",quote = F)

# Read in the txt file
dat <- read.table('CRC@assays$Spatial@countsPC.txt', header = TRUE)

# Extract first row and row with related genes (for Figure 4i)
dat1 <- dat[c(which(dat$ID == "PECAM1" | dat$ID == "PGF"|dat$ID == "KDR" | dat$ID == "VEGFA")),]

# Write the output as a txt file
write.table(dat1, file = "output.txt", sep="\t", row.names = FALSE)
