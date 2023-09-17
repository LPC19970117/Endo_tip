#每次运行都需要先加载相关的包才能调用后面的分析脚本
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)

#这里文件夹的名字可以修改，但最好只用英文字母
sam.name <- "GSE205506united1"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}

#### 2. 读入原始表达数据 ####
#以下两种方式二选一
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample111/")
sce111 <- CreateSeuratObject(experiment.data)
sce111@meta.data$sample="111"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample112/")
sce112 <- CreateSeuratObject(experiment.data)
sce112@meta.data$sample="112"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample121/")
sce121 <- CreateSeuratObject(experiment.data)
sce121@meta.data$sample="121"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample122/")
sce122 <- CreateSeuratObject(experiment.data)
sce122@meta.data$sample="122"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample141/")
sce141 <- CreateSeuratObject(experiment.data)
sce141@meta.data$sample="141"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample151/")
sce151 <- CreateSeuratObject(experiment.data)
sce151@meta.data$sample="151"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample152/")
sce152 <- CreateSeuratObject(experiment.data)
sce152@meta.data$sample="152"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample171/")
sce171 <- CreateSeuratObject(experiment.data)
sce171@meta.data$sample="171"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample172/")
sce172 <- CreateSeuratObject(experiment.data)
sce172@meta.data$sample="172"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample181/")
sce181 <- CreateSeuratObject(experiment.data)
sce181@meta.data$sample="181"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample182/")
sce182 <- CreateSeuratObject(experiment.data)
sce182@meta.data$sample="182"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample191/")
sce191 <- CreateSeuratObject(experiment.data)
sce191@meta.data$sample="191"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample192/")
sce192 <- CreateSeuratObject(experiment.data)
sce192@meta.data$sample="192"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample211/")
sce211 <- CreateSeuratObject(experiment.data)
sce211@meta.data$sample="211"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample213/")
sce213 <- CreateSeuratObject(experiment.data)
sce213@meta.data$sample="213"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample233/")
sce233 <- CreateSeuratObject(experiment.data)
sce233@meta.data$sample="233"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample241/")
sce241 <- CreateSeuratObject(experiment.data)
sce241@meta.data$sample="241"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample243/")
sce243 <- CreateSeuratObject(experiment.data)
sce243@meta.data$sample="243"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample251/")
sce251 <- CreateSeuratObject(experiment.data)
sce251@meta.data$sample="251"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample252/")
sce252 <- CreateSeuratObject(experiment.data)
sce252@meta.data$sample="252"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample253/")
sce253 <- CreateSeuratObject(experiment.data)
sce253@meta.data$sample="253"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample261/")
sce261 <- CreateSeuratObject(experiment.data)
sce261@meta.data$sample="261"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample271/")
sce271 <- CreateSeuratObject(experiment.data)
sce271@meta.data$sample="271"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample272/")
sce272 <- CreateSeuratObject(experiment.data)
sce272@meta.data$sample="272"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample273/")
sce273 <- CreateSeuratObject(experiment.data)
sce273@meta.data$sample="273"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample281/")
sce281 <- CreateSeuratObject(experiment.data)
sce281@meta.data$sample="281"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample282/")
sce282 <- CreateSeuratObject(experiment.data)
sce282@meta.data$sample="282"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample283/")
sce283 <- CreateSeuratObject(experiment.data)
sce283@meta.data$sample="283"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample291/")
sce291 <- CreateSeuratObject(experiment.data)
sce291@meta.data$sample="291"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample292/")
sce292 <- CreateSeuratObject(experiment.data)
sce292@meta.data$sample="292"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample301/")
sce301 <- CreateSeuratObject(experiment.data)
sce301@meta.data$sample="301"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample302/")
sce302 <- CreateSeuratObject(experiment.data)
sce302@meta.data$sample="302"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample303/")
sce303 <- CreateSeuratObject(experiment.data)
sce303@meta.data$sample="303"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample311/")
sce311<- CreateSeuratObject(experiment.data)
sce311@meta.data$sample="311"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample312/")
sce312 <- CreateSeuratObject(experiment.data)
sce312@meta.data$sample="312"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample313/")
sce313 <- CreateSeuratObject(experiment.data)
sce313@meta.data$sample="313"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample321/")
sce321 <- CreateSeuratObject(experiment.data)
sce321@meta.data$sample="321"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample322/")
sce322 <- CreateSeuratObject(experiment.data)
sce322@meta.data$sample="322"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample323/")
sce323 <- CreateSeuratObject(experiment.data)
sce323@meta.data$sample="323"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)
#如果下再的是10X格式的数据，则用下面代码读入
experiment.data <- Read10X("./sample333/")
sce333 <- CreateSeuratObject(experiment.data)
sce333@meta.data$sample="333"
#write.table(sce@meta.data,file=paste0("./",sam.name,"/meta.data.txt"),sep="\t",quote = F)

#合并数据
sce1 <- merge(sce111, y=c( sce112, sce121, sce122, sce141,sce151, sce152, sce171, sce172, sce181, sce182,
                           sce191, sce192, sce211, sce213,sce233, sce241, sce243, sce251, sce252, sce253,
                           sce261, sce271, sce272, sce273,sce281, sce282, sce283, sce291, sce292, sce301,
                           sce302, sce303, sce311, sce312,sce313, sce321, sce322, sce323, sce333))
#4.存储数据
save(sce1,file=paste0("./",sam.name,"/sce1_united.RData"))

sce<-sce1
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
pdf(paste0("./",sam.name,"/QC-VlnPlot.pdf"),width = 8,height = 4.5)
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)
dev.off()

plot1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pdf(paste0("./",sam.name,"/QC-FeatureScatter.pdf"),width = 8,height = 4.5)
CombinePlots(plots = list(plot1, plot2),legend = "none")
dev.off()
rm(plot1,plot2)

#### 5. 筛选cell ####
cat("Before filter :",nrow(sce@meta.data),"cells\n")
#sce <- subset(sce, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20)
#cat("After filter :",nrow(sce@meta.data),"cells\n")

### Normalization
#############################################################################################################
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)  # Checkpoint

### Feature selection
#############################################################################################################
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sce), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sce)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf(file = paste0(sam.name,"/Norm-feature_variable_plot.pdf"),width = 8,height = 5)
CombinePlots(plots = list(plot1, plot2),legend = "none")
dev.off()
rm(plot1,plot2)

#均一化（需要一点时间）
#这是2000高变基因的归一化
sce <- ScaleData(sce,
                 vars.to.regress = c("percent.mt"))
#############################################################################################################
## Reduction
#############################################################################################################
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
pdf(paste0("./",sam.name,"/PCA-ElbowPlot.pdf"),width = 5,height = 4)
ElbowPlot(sce, reduction = "pca",ndims = 50)
dev.off()

### Cluster
dim.use <- 1:50
# ############################################################################################################
sce <- RunUMAP(sce, dims = dim.use) # umap
sce <- FindNeighbors(sce, dims = 1:50)  # louvain cluster, graph based
sce <- FindClusters(sce, resolution = 0.8)

##seurat_clusters_DimPlot
pdf(paste0("./",sam.name,"/CellCluster-DimPlot_umap_",max(dim.use),".pdf"),width = 10,height = 9)
DimPlot(sce, reduction = "umap",group.by = 'seurat_clusters',label = T)
dev.off()

##major_marker_gene_umap
pdf(paste0("./",sam.name,"/umap.",max(dim.use),"PC.pdf"),width = 15,height = 14)
FeaturePlot(sce, features = c("CD3D","KLRF1","MS4A1","MZB1","LYZ","TPSAB1","EPCAM","DCN","TNFRSF17","COL1A1","COL1A2","PLVAP","LILRA4","FCGR3B"),
            cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,
            reduction = "umap")
dev.off()

##major_marker_gene_umap2
pdf(paste0("./",sam.name,"/GSE178341_major.",max(dim.use),"PC.pdf"),width = 15,height = 20)
FeaturePlot(sce, features <- c("TFF3","AGR2","EPCAM",#Epithelial
                               "CD3E","CD3D","CCL5",#T
                               "PLVAP","VWF","PECAM1",#
                               "JCHAIN","MZB1","IGHA1",#plasma
                               "COL1A1","COL1A2","DCN",#Fibroblast
                               "LYZ","CD68","C1QC",#Myeloid
                               "MS4A1","CD79A","CD19",#B
                               "TPSB2","TPSAB1","CPA3"#Mast
),
cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,
reduction = "umap")
dev.off()

## major_marker_gene 气泡图
features <- c("TFF3","AGR2","EPCAM",#Epithelial
              "CD3E","CD3D","CCL5",#T
              "PLVAP","VWF","PECAM1",#
              "JCHAIN","MZB1","IGHA1",#plasma
              "COL1A1","COL1A2","DCN",#Fibroblast
              "LYZ","CD68","C1QC",#Myeloid
              "MS4A1","CD79A","CD19",#B
              "TPSB2","TPSAB1","CPA3"#Mast
)
pdf(paste0("./",sam.name,"/DotPlot_umap.",max(dim.use),"PC.pdf"),width = 10,height = 10)
DotPlot(sce, features = unique(features)) + RotatedAxis()
dev.off()

##细胞类群相似树
sce <- BuildClusterTree(
  sce,
  dims = dim.use,
  reorder = F,
  reorder.numeric = F)
pdf(paste0("./",sam.name,"/CellCluster-ClusterTree_",max(dim.use),"PC.pdf"),width = 10,height = 8)
PlotClusterTree(sce)
dev.off()

#### 2.1 计算major_marker gene
#sce= SetIdent(sce,value="celltype") 
all.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/","main_marker_genes",max(dim.use),"PC.txt"),sep="\t",quote = F)

table(sce@meta.data$seurat_clusters)

##在meta.data中重命名seurat_clusters为celltype
#方法二：使用unname函数配合向量：
cluster2celltype <- c( "0"="Epithelial",
                       "1"="B", 
                       "2"="Epithelial", 
                       "3"= "T", 
                       "4"= "T", 
                       "5"= "Endothelial",
                       "6"= "Epithelial", 
                       "7"= "T", 
                       "8"= "Myeloid",
                       "9"= "Epithelial",
                       "10"="Epithelial",
                       "11"="Plasma", 
                       "12"="Epithelial", 
                       "13"= "Epithelial", 
                       "14"= "Epithelial", 
                       "15"= "Epithelial",
                       "16"= "Epithelial", 
                       "17"= "Fibroblast", 
                       "18"= "Epithelial",
                       "19"= "Epithelial",
                       "20"="Epithelial",
                       "21"="Epithelial", 
                       "22"="B", 
                       "23"= "Epithelial",
                       "24"= "Fibroblast", 
                       "25"= "Epithelial",
                       "26"= "Epithelial", 
                       "27"= "Epithelial", 
                       "28"= "Epithelial",
                       "29"= "Epithelial",
                       "30"="Epithelial",
                       "31"="T", 
                       "32"="T", 
                       "33"= "Epithelial", 
                       "34"= "Epithelial", 
                       "35"= "Endothelial",
                       "36"= "Epithelial", 
                       "37"= "Epithelial", 
                       "38"= "Epithelial",
                       "39"= "Epithelial",
                       "40"="T",
                       "41"="B", 
                       "42"="Mast", 
                       "43"= "Plasma", 
                       "44"= "Epithelial", 
                       "45"= "Endothelial"
)
sce[['celltype']] = unname(cluster2celltype[sce@meta.data$seurat_clusters])
##celltype UMAP图
pdf(paste0("./",sam.name,"/celltype-DimPlot_umap_",max(dim.use),".pdf"),width = 5,height = 4)
DimPlot(sce, reduction = "umap",group.by = 'celltype',label = T)
dev.off()

table(sce@meta.data$celltype)
#4.存储数据
save(sce,file=paste0("./",sam.name,"/sce1.1_united.RData"))
#### 2.1 计算major_marker gene
sce= SetIdent(sce,value="celltype") 
all.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/","main_celltype.txt"),sep="\t",quote = F)
#输出meta.data临床数据，在excel中修改补充
write.table(sce@meta.data,file=paste0("./",sam.name,"/","sce@meta.data",max(dim.use),"PC.txt"),sep="\t",quote = F)
#补充meta.data临床数据(txt)
ldf=read.table("patient.txt",header = T)[,1]
sce@meta.data$patient <- ldf
ldf=read.table("tissue.txt",header = T)[,1]
sce@meta.data$tissue <- ldf
ldf=read.table("response.txt",header = T)[,1]
sce@meta.data$response <- ldf
ldf=read.table("treatment.txt",header = T)[,1]
sce@meta.data$treatment<- ldf
ldf=read.table("group.txt",header = T)[,1]
sce@meta.data$group <- ldf
rm(ldf)

##临床参数的UMAP图
pdf(paste0("./",sam.name,"/umap_patient.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="patient",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/umap_tissue.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="tissue",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/umap_group.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="group",reduction='umap')
dev.off()

#4.存储数据
save(sce,file=paste0("./",sam.name,"/sce1.1_united.RData"))

#4.提取Endothelial亚组
sce_Endo<-subset(sce,celltype %in% c("Endothelial"))
table(sce_Endo@meta.data$celltype)
#4.存储Endothelial数据
save(sce_Endo,file=paste0("./",sam.name,"/","Endothelial.RData"))
