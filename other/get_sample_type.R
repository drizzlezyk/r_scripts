library(DropSeq.util)
dge.path <- "/Users/zhongyuanke/data/scRNA-seq_mouse/Cerebellum/F_GRCm38.81.P60Cerebellum_ALT.raw.dge.txt.gz"
dge <- loadSparseDge(dge.path)

label_path <- '/Users/zhongyuanke/data/scRNA-seq_mouse/annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS'
cluster_path <-  '/Users/zhongyuanke/data/scRNA-seq_mouse/Cerebellum/F_GRCm38.81.P60Cerebellum_ALT.cell_cluster_outcomes.RDS'
cluster_sign_path <- '/Users/zhongyuanke/data/scRNA-seq_mouse/Cerebellum/F_GRCm38.81.P60Cerebellum_ALT.cluster.assign.RDS'
celltype_path <- '/Users/zhongyuanke/data/atac/p0_BrainCortex/ref_barcodes.csv'
label <- readRDS(label_path)
clusters <- readRDS(cluster_path)
cluster_sign <- readRDS(cluster_sign_path)