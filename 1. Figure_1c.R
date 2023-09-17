library(dplyr)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(ggrepel) #用于注释文本
library(magrittr)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
# 读取3个数据集的差异基因文件
pbmc.markers <- read.csv("output_pbmc.markers.csv")
table(pbmc.markers$cluster)
# 读取修改后的要标注的基因名文件
top_label <- read.csv("easy_input_label.csv")
### 准备绘制暗灰色背景所需数据 
background_position <- pbmc.markers %>%
  group_by(cluster) %>%
  summarise(Min = min(avg_log2FC) - 0.2, Max = max(avg_log2FC) + 0.2) %>%
  as.data.frame()

background_position$cluster %<>% as.vector(.) %>% as.numeric(.)
background_position$start <- background_position$cluster - 0.4
background_position$end <- background_position$cluster + 0.4

### 准备绘制中间区域cluster彩色bar所需数据
cluster_bar_position <- background_position
cluster_bar_position$start <- cluster_bar_position$cluster - 0.5
cluster_bar_position$end <- cluster_bar_position$cluster + 0.5
cluster_bar_position$cluster %<>% 
  factor(., levels = c(0:max(as.vector(.))))

## 设置填充颜色
cols_thr_signi <- c("Up_Highly" = "#d7301f",
                    "Down_Highly" = "#225ea8",
                    "Up_Lowly" = "black",
                    "Down_Lowly" = "black")
cols_cluster <- c("0" = "#35978f",
                  "1" = "#8dd3c7",
                  "2" = "#ffffb3",
                  "3" = "#bebada",
                  "4" = "#fb8072",
                  "5" = "#80b1d3",
                  "6" = "#fdb462",
                  "7" = "#b3de69",
                  "8" = "#fccde5")

p <- ggplot() +
  geom_rect(data = background_position, aes(xmin = start, xmax = end, ymin = Min,
                                            ymax = Max),
            fill = "#525252", alpha = 0.1) + ###添加灰色背景色
  geom_jitter(data = pbmc.markers, aes(x = cluster, y = avg_log2FC, colour = thr_signi),
              size = 1,position = position_jitter(seed = 1)) +
  scale_color_manual(values = cols_thr_signi) +
  scale_x_continuous(limits = c(-0.5, max(pbmc.markers$cluster) + 0.5),
                     breaks = seq(0, max(pbmc.markers$cluster), 1),
                     label = seq(0, max(pbmc.markers$cluster),1)) + #修改坐标轴显示刻度
  
##根据top_label标注基因名
geom_text_repel(data = top_label, aes(x = cluster, y = avg_log2FC, label = gene),
                  position = position_jitter(seed = 1), show.legend = F, size = 2.5,
                  box.padding = unit(0, "lines")) +
  
  geom_rect(data = cluster_bar_position, aes(xmin = start, xmax = end, ymin = -0.4,
                                             ymax = 0.4, fill = cluster), color = "black", alpha = 1, show.legend = F) +
  scale_fill_manual(values = cols_cluster) +
  labs(x = "Cluster", y = "average log2FC") +
  theme_bw()

plot1 <- p + theme(panel.grid.minor = element_blank(), ##去除网格线
                   panel.grid.major = element_blank(),
                   axis.text.y = element_text(colour = 'black', size = 14),
                   axis.text.x = element_text(colour = 'black', size = 14, vjust = 50), #调整x轴坐标,vjust的值按照最终结果稍加调整
                   panel.border = element_blank(), ## 去掉坐标轴
                   axis.ticks.x = element_blank(), ## 去掉的坐标刻度线
                   axis.line.y = element_line(colour = "black")) #添加y轴坐标轴
ggsave(filename = "Marker_gene_pointplot.pdf", plot = plot1, width = 9, height = 6)

