library(Seurat)
library(SeuratData)
data("ifnb")
ifnb.list <- SplitObject(ifnb, split.by = "stim")
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- subset(x, features =c(VariableFeatures(x) ))
})

control_data<-GetAssayData(ifnb.list[["CTRL"]])
control_data <- as.matrix(control_data)
control_data <- t(control_data)
data_path <- '/Users/zhongyuanke/data/seurat_data/ifnb/control_data_norm.csv'
write.csv(control_data, data_path)

ctrl_celltype <- ifnb.list[["CTRL"]]@meta.data[["seurat_annotations"]]
celltype_path <- '/Users/zhongyuanke/data/seurat_data/ifnb/control_cell_type.csv'
write.csv(ctrl_celltype, celltype_path)

stim_data<-GetAssayData(ifnb.list[["STIM"]])
stim_data <- as.matrix(stim_data)
stim_data <- t(stim_data)
stim_data_path <- '/Users/zhongyuanke/data/seurat_data/ifnb/stim_data_norm.csv'
write.csv(stim_data, stim_data_path)

stim_celltype <- ifnb.list[["STIM"]]@meta.data[["seurat_annotations"]]
stim_celltype_path <- '/Users/zhongyuanke/data/seurat_data/ifnb/stim_cell_type.csv'
write.csv(stim_celltype, stim_celltype_path)