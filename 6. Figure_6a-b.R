#终端部分################################################################################################################################################################
# 创建名为cellphonedb的虚拟环境
#conda create -n cellphonedb python=3.7 
# 激活虚拟环境
#conda activate cellphonedb 
# 在虚拟环境中下载软件
#pip install cellphonedb
# 假如网络不行，就加上  -i http://mirrors.aliyun.com/pypi/simple/   
# 很简单就安装成功， 试试看运行它 获取帮助信息
#cellphonedb --help
########################################################################################################################################################################

#R部分由于细胞类型已经注释好，接下来准备cellphonedb的文件：
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(hdf5r)
load(file = "sce_500.RData")
sce<-sce_500
#表达谱文件cellphonedb_count.txt和细胞分组注释文件cellphonedb_meta.txt。
write.table(as.matrix(sce@assays$RNA@data), 'cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(sce@meta.data), 500@meta.data[,'subcelltype2', drop=F])  #cell_type改成subcelltype2，之后在txt中改回
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA
write.table(meta_data, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)

#终端########################################################################################################################################################################
#先定位 cd XXXXXXX
#cellphonedb method statistical_analysis  cellphonedb_meta.txt  cellphonedb_count.txt      --counts-data=gene_name

#$cellphonedb plot dot_plot$cellphonedb plot heatmap_plot cellphonedb_meta.txt
#如果我们count的基因是基因名格式，需要添加参数--counts-data=gene_name，如果行名为ensemble名称的话，可以不添加这个参数，使用默认值即可。

#下面没卵用
#conda activate cellphonedb
# 必须要保证当前路径下面有前面的步骤输出的out文件夹哦 
#cellphonedb plot dot_plot 
#cellphonedb plot heatmap_plot cellphonedb_meta.txt 
########################################################################################################################################################################

# 做完后为了跟前面的区分，我把 out文件夹，修改名字为 out-by-symbol 文件夹啦 

#R 热图 可视化
library(tidyverse)
library(RColorBrewer)
library(scales)

pvalues=read.table("./out/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
pvalues=pvalues[,12:dim(pvalues)[2]] #此时不关注前11列
statdf=as.data.frame(colSums(pvalues < 0.05)) #统计在某一种细胞pair的情况之下，显著的受配体pair的数目；阈值可以自己选
colnames(statdf)=c("number")

#排在前面的分子定义为indexa；排在后面的分子定义为indexb
statdf$indexb=str_replace(rownames(statdf),"^.*\\.","")
statdf$indexa=str_replace(rownames(statdf),"\\..*$","")
#设置合适的细胞类型的顺序
rankname=sort(unique(statdf$indexa)) 
#转成因子类型，画图时，图形将按照预先设置的顺序排列
statdf$indexa=factor(statdf$indexa,levels = rankname)
statdf$indexb=factor(statdf$indexb,levels = rankname)

statdf%>%ggplot(aes(x=indexa,y=indexb,fill=number,display_numbers = T))+geom_tile(color="white")+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,165))+
  scale_x_discrete("cluster 1 produces molecule 1")+
  scale_y_discrete("cluster 2 produces molecule 2")+
  theme_minimal()+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = NULL, angle = 45),
    panel.grid = element_blank(),display_numbers = T
  )+RotatedAxis()
ggsave(filename = "interaction.num.1.pdf",device = "pdf",width = 12,height = 10,units = c("cm"))

#对称热图
library(tidyverse)
library(RColorBrewer)
library(scales)

pvalues=read.table("./out/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
pvalues=pvalues[,12:dim(pvalues)[2]]
statdf=as.data.frame(colSums(pvalues < 0.1))
colnames(statdf)=c("number")
statdf$indexb=str_replace(rownames(statdf),"^.*\\.","")
statdf$indexa=str_replace(rownames(statdf),"\\..*$","")
statdf$total_number=0

for (i in 1:dim(statdf)[1]) {
  tmp_indexb=statdf[i,"indexb"]
  tmp_indexa=statdf[i,"indexa"]
  if (tmp_indexa == tmp_indexb) {
    statdf[i,"total_number"] = statdf[i,"number"]
  } else {
    statdf[i,"total_number"] = statdf[statdf$indexb==tmp_indexb & statdf$indexa==tmp_indexa,"number"]+
      statdf[statdf$indexa==tmp_indexb & statdf$indexb==tmp_indexa,"number"]
  }
}
statdf$indexa
rankname=sort(unique(statdf$indexa)) 
statdf$indexa=factor(statdf$indexa,levels = rankname)
statdf$indexb=factor(statdf$indexb,levels = rankname)
table(statdf$total_number)
statdf$total_number
statdf%>%ggplot(aes(x=indexa,y=indexb,fill=total_number,display_numbers = T))+geom_tile(color="white")+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,270))+#limits=c(-5,100))能显著改变颜色
  scale_x_discrete("cluster 1")+display_numbers = T
  scale_y_discrete("cluster 2")+
  theme_minimal()+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = NULL, angle = 45),
    panel.grid = element_blank(60)
  )+RotatedAxis()
ggsave(filename = "interaction.num.3——.pdf",device = "pdf",width = 15,height = 10,units = c("cm"))

write.table(statdf$total_number,file=paste0("./statdf.txt"),sep="\t",quote = F)

#气泡图
source("CCC.bubble.R")
CCC(
  pfile="./out/pvalues.txt",
  mfile="./out/means.txt",
  #neg_log10_th= -log10(0.01),
  means_exp_log2_th=1,
  neg_log10_th2=3,
  #means_exp_log2_th2=c(-4,6),
  #notused.cell=c("B"),
  #used.cell=c("Mcell"),
  #cell.pair=c("Mcell.Scell","Mcell.NKcell","Mcell.Tcell","Scell.Mcell","NKcell.Mcell","Tcell.Mcell"),#这里是自定义的顺序，若是可选细胞对的子集，则只展示子集，若有交集则只展示交集；空值情况下，会根据可选细胞对自动排序
  gene.pair=c("VEGFA_KDR",#"VEGFA_FLT1",
              "NRP1_VEGFA","NRP2_VEGFA")#作用同上
)
ggsave(filename = "VEGFA_.pdf",device = "pdf",width =20,height = 12,units = "cm")
#再细化
CCC(
  pfile="./out/pvalues.txt",
  mfile="./out/means.txt",
  cell.pair=c("Mcell.Scell","Mcell.NKcell","Mcell.Tcell","Scell.Mcell","NKcell.Mcell","Tcell.Mcell"),#这里是自定义的顺序，若是可选细胞对的子集，则只展示子集，若有交集则只展示交集；空值情况下，会根据可选细胞对自动排序
  gene.pair=c("VEGFA_KDR","VEGFA_FLT1","NRP1_VEGFA","NRP2_VEGFA")#作用同上
)
ggsave(filename = "interaction.detail.2.pdf",device = "pdf",width =16,height = 10,units = "cm")




####
#安装
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
devtools::install_github('zktuong/ktplots', dependencies = TRUE)
library(ktplots)
#作图需要的文件，means.txt， pvalues.txt, 单细胞seurat对象
load("sce_Mast_ar.RData")
pvals <- read.delim(paste0(sce,"pvalues.txt"), check.names = FALSE)
means <- read.delim(paste0(sce,"means.txt"), check.names = FALSE)




 pvals <- read.delim("pvalues.txt", check.names = FALSE)
 means <- read.delim("means.txt", check.names = FALSE)
pdf(paste0("./","/cellphoneDB_Endo_Tip.",".pdf"),width = 10,height = 20)
 plot_cpdb(cell_type1 = "Endo_Tip", cell_type2 = "Myeloid", scdata = sce_500,
           idents = 'subcelltype2', means = means, pvals = pvals, 
           #gene.family = 'costimulatory',
           #split.by = 'subcelltype2',
           genes = c("VEGF"),
           #col_option = "maroon",
           col_option =c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"),
           highlight = "red",
           keep_significant_only=T) +
   theme(axis.text  = element_text(size = 10, color = 'black'))
 dev.off()
 
 pdf(paste0("./","/cellphoneDB_Endo_EPi_myeloid.",".pdf"),width = 10,height = 20)
 plot_cpdb(cell_type1 = "Endo", cell_type2 = c("Myeloid","Epithelial"), scdata = sce_500,
           idents = 'subcelltype2', means = means, pvals = pvals, 
           #gene.family = 'costimulatory',
           #split.by = 'subcelltype2',
           genes = c("VEGFA"),
           #col_option = "maroon",
           col_option =c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"),
           highlight = "red",
           keep_significant_only=T) +
   theme(axis.text  = element_text(size = 10, color = 'black'))
 dev.off()
 
 pdf(paste0("./","/cellphoneDB_Endo_Epithelial.",".pdf"),width = 10,height = 20)
 plot_cpdb(cell_type1 = "Endo", cell_type2 = c("Epithelial"), scdata = sce_500,
           idents = 'subcelltype2', means = means, pvals = pvals, 
           #gene.family = 'costimulatory',
           #split.by = 'subcelltype2',
           genes = c("VEGFA"),
           #col_option = "maroon",
           col_option =c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"),
           highlight = "red",
           keep_significant_only=T) +
   theme(axis.text  = element_text(size = 10, color = 'black'))
 dev.off()
 
 pdf(paste0("./","/cellphoneDB_Endo_Epithelial.",".pdf"),width =6,height = 4)
 plot_cpdb(cell_type1 = "Endo", cell_type2 ="Epithelial", scdata = sce_500,
           idents = 'subcelltype2', means = means, pvals = pvals, 
           #gene.family = 'costimulatory',
           #split.by = 'subcelltype2',
           genes = c("VEGFA"),
           #col_option = "maroon",
           col_option =c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"),
           highlight = "red",
           keep_significant_only=T) +
   theme(axis.text  = element_text(size = 10, color = 'black'))
 dev.off()
 
 pdf(paste0("./","/cellphoneDB_Endo_Myeloid.",".pdf"),width =6,height = 4)
 plot_cpdb(cell_type1 = "Endo", cell_type2 ="Myeloid", scdata = sce_500,
           idents = 'subcelltype2', means = means, pvals = pvals, 
           #gene.family = 'costimulatory',
           #split.by = 'subcelltype2',
           genes = c("VEGFA"),
           #col_option = "maroon",
           col_option =c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"),
           highlight = "red",
           keep_significant_only=T) +
   theme(axis.text  = element_text(size = 10, color = 'black'))
 dev.off()
 pdf(paste0("./","/cellphoneDB_Endo_Fibroblast.",".pdf"),width =6,height = 4)
 plot_cpdb(cell_type1 = "Endo", cell_type2 ="Fibroblast", scdata = sce_500,
           idents = 'subcelltype2', means = means, pvals = pvals, 
           #gene.family = 'costimulatory',
           #split.by = 'subcelltype2',
           genes = c("VEGFA"),
           #col_option = "maroon",
           col_option =c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"),
           highlight = "red",
           keep_significant_only=T) +
   theme(axis.text  = element_text(size = 10, color = 'black'))
 dev.off()
 
 pdf(paste0("./","/cellphoneDB_Endo_Myeloid_KDR.",".pdf"),width =6,height = 4)
 plot_cpdb(cell_type1 = "Endo", cell_type2 ="Myeloid", scdata = sce_500,
           idents = 'subcelltype2', means = means, pvals = pvals, 
           #gene.family = 'costimulatory',
           #split.by = 'subcelltype2',
           genes = c("KDR"),
           #col_option = "maroon",
           col_option =c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"),
           highlight = "red",
           keep_significant_only=T) +
   theme(axis.text  = element_text(size = 10, color = 'black'))
 dev.off()
 