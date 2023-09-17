#提取endo特定亚群
sce_BVEC<-subset(sce_Endo,subcelltype2 %in% c("Tip","Immature","Artery","Capillary","Vein"))

write.table(sce_BVEC@assays$RNA@counts,file=paste0("./",sam.name,"/","sce_BVEC@assays$RNA@counts.txt"),sep="\t",quote = F)
write.table(sce_BVEC@meta.data,file=paste0("./",sam.name,"/","sce_BVEC@meta.data.txt"),sep="\t",quote = F)
table(sce_BVEC@meta.data$subcelltype2)
# 每个细胞亚群抽500 
sce_BVEC= SetIdent(sce_BVEC,value="subcelltype2")
allCells=names(Idents(sce_BVEC))
allType = levels(Idents(sce_BVEC))
choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(sce_BVEC)== x ]
  cg=sample(cgCells,500)
  cg
}))
cg_sce = sce_BVEC[, allCells %in% choose_Cells]
cg_sce
as.data.frame(table(Idents(cg_sce)))

write.table(cg_sce@assays$RNA@counts,file=paste0("./",sam.name,"/","cg_sce@assays$RNA@counts.txt"),sep="\t",quote = F)
write.table(cg_sce@meta.data,file=paste0("./",sam.name,"/","cg_sce@meta.data.txt"),sep="\t",quote = F)
table(cg_sce@meta.data$subcelltype2)
table(sce_BVEC@meta.data$subcelltype2)
