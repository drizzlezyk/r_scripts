library(Seurat)
library(ggplot2)
library(patchwork)
base_path <- '/Users/zhongyuanke/data/seurat_data/'
atac_path <- paste(base_path, 'sc_atac/atac_pbmc_10k_v1_filtered_peak_bc_matrix.h5' ,sep = "")
peaks <- Read10X_h5(atac_path)

activity.matrix <- CreateGeneActivityMatrix(peak.matrix = peaks, 
                                            annotation.file = "/Users/zhongyuanke/data/seurat_data/sc_atac/Homo_sapiens.GRCh37.87.gtf", 
                                            seq.levels = c(1:22, "X", "Y"), upstream = 2000, verbose = TRUE)
pbmc.atac <- CreateSeuratObject(counts = peaks, assay = "ATAC", project = "10x_ATAC")
pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
meta <- read.table("/Users/zhongyuanke/data/sc_atac/atac_v1_pbmc_10k_singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                   stringsAsFactors = FALSE)
meta <- meta[colnames(pbmc.atac), ]
pbmc.atac <- AddMetaData(pbmc.atac, metadata = meta)
pbmc.atac <- subset(pbmc.atac, subset = nCount_ATAC > 5000)
pbmc.atac$tech <- "atac"

DefaultAssay(pbmc.atac) <- "ACTIVITY"
pbmc.atac <- FindVariableFeatures(pbmc.atac)
pbmc.atac <- NormalizeData(pbmc.atac)
pbmc.atac <- ScaleData(pbmc.atac)

atac_scale <- pbmc.atac[["ACTIVITY"]]@scale.data
atac_scale <- t(atac_scale)
atac_path <- '/Users/zhongyuanke/data/sc_atac/atac_gene_activity.csv'
write.csv(atac_scale, atac_path)


pbmc.rna <- readRDS("/Users/zhongyuanke/data/sc_atac/pbmc_10k_v3.rds")
pbmc.rna$tech <- "rna"
rna_scale <- pbmc.rna[["RNA"]]@scale.data
rna_scale <- t(rna_scale)
celltype <- pbmc.rna@meta.data[["celltype"]]
celltype_path <- '/Users/zhongyuanke/data/sc_atac/pbmc_10k_v3_celltype.csv'
rna_path <- '/Users/zhongyuanke/data/sc_atac/pbmc_10k_v3_scale.csv'
write.csv(celltype, celltype_path)
write.csv(rna_scale, rna_path)

