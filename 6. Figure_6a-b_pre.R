rm(list=ls())
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(hdf5r)
#### input data
load(file = 'sce1.RData')

##output data file
sam.name <- "Major_抽样500_CellphoneDB"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}

#补充meta.data临床数据(txt)_subcelltype2
ldf=read.table("subcelltype2.txt",header = T)[,1]
sce@meta.data$subcelltype2 <- ldf
#输出meta.data临床数据，在excel中修改补充
write.table(sce@meta.data,file=paste0("./",sam.name,"/","sce@meta.data.txt"),sep="\t",quote = F)
##major_marker_gene点状热图
library(RColorBrewer)
sce= SetIdent(sce,value="celltype")
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
features <- c("KRT8","PHGR1","EPCAM",#Epithelial
              "CD7","CD3D","CCL5",#T
              "COL1A1","COL1A2","DCN",#Fibroblast
              "CD68","TYROBP","LYZ",#Myeloid
              "PLVAP","VWF","CLDN5",#Endo
              "MS4A1","CD79A","CD79B",#B
              "TPSB2","TPSAB1","CPA3",#Mast
              "JCHAIN","MZB1","IGHA2"#plasma
)
pdf(paste0("./",sam.name,"/niche_marker_DotPlot_CRC1",".pdf"),width =4.2,height = 5.2)
DotPlot(sce, features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+RotatedAxis()
dev.off()
# 每个细胞亚群抽50 
# 每个细胞亚群抽500 
sce= SetIdent(sce,value="subcelltype2")
allCells=names(Idents(sce))
allType = levels(Idents(sce))
choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(sce)== x ]
  cg=sample(cgCells,500)
  cg
}))
cg_sce = sce[,allCells %in% choose_Cells]
cg_sce
as.data.frame(table(Idents(cg_sce)))

#将数据写到文件中一边后续分析
#提取endo特定亚群
sce_500<-subset(cg_sce,subcelltype2 %in% c("Epithelial","T","Plasma","Fibroblast","Myeloid","B",
                                           "Mast","Endo_Tip","Endo_Immature","Endo_Artery","Endo_Capillary","Endo_Vein"))
as.data.frame(table(Idents(sce_500)))
save(sce_500,file=paste0("./",sam.name,"/","sce_500.RData"))
