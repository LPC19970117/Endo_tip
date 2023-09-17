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
