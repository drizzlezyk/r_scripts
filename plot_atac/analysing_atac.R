library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
library(zoo)
library(dplyr)
library(gggenes)
set.seed(1234)

base_path <- '/Users/zhongyuanke/data/seurat_data/'
counts_path <- paste(base_path, 'sc_atac/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5' ,sep = "")
meta_path <- paste(base_path, 'sc_atac/atac_v1_pbmc_10k_singlecell.csv' ,sep = "")
rna_path <- paste(base_path, 'sc_atac/pbmc_10k_v3.rds' ,sep = "")
fragment.path <- paste(base_path, 'seurat_data/sc_atac/atac_v1_pbmc_10k_fragments.tsv.gz' ,sep = "")

counts <- Read10X_h5(filename = counts_path)

metadata <- read.csv(
  file = meta_path,
  header = TRUE,
  row.names = 1
)

pbmc <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = metadata
)



pbmc <- SetFragments(
  object = pbmc,
  file = fragment.path
)

pbmc <- NucleosomeSignal(object = pbmc)
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')

gene.ranges <- genes(EnsDb.Hsapiens.v75)
seqlevelsStyle(gene.ranges) <- 'UCSC'
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')

tss.ranges <- resize(gene.ranges, width = 1, fix = "start")
seqlevelsStyle(tss.ranges) <- 'UCSC'
tss.ranges <- keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse')

# to save time use the first 2000 TSSs
pbmc <- TSSEnrichment(object = pbmc, tss.positions = tss.ranges[1:2000])
pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')

pbmc <- subset(
  x = pbmc,
  subset = peak_region_fragments > 3000 &
  peak_region_fragments < 20000 &
  pct_reads_in_peaks > 15 &
  blacklist_ratio < 0.05 &
  nucleosome_signal < 10 &
  TSS.enrichment > 2
)

# pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(
  object = pbmc,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)


pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)

genebodyandpromoter.coords <- Extend(x = gene.ranges, upstream = 2000, downstream = 0)

# create a gene by cell matrix
gene.activities <- FeatureMatrix(
  fragments = fragment.path,
  features = genebodyandpromoter.coords,
  cells = colnames(pbmc),
  chunk = 20
)

# convert rownames from chromsomal coordinates into gene names
gene.key <- genebodyandpromoter.coords$gene_name
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
rownames(gene.activities) <- gene.key[rownames(gene.activities)]

# add the gene activity matrix to the Seurat object as a new assay, and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
# pbmc <- NormalizeData(
#   object = pbmc,
#   assay = 'RNA',
#   normalization.method = 'LogNormalize',
#   scale.factor = median(pbmc$nCount_RNA)
# )

DefaultAssay(pbmc) <- 'RNA'

pbmc_anchor <- readRDS(anchor_path)

transfer.anchors <- FindTransferAnchors(
  reference = pbmc_anchor,
  query = pbmc,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_anchor$celltype,
  weight.reduction = pbmc[['lsi']],
  dims = 2:30
)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)

pbmc <- subset(pbmc,idents = 14, invert = TRUE)
pbmc <- RenameIdents(
  object = pbmc,
  '0' = 'CD14 Mono',
  '1' = 'CD4 Memory',
  '2' = 'CD8 Effector',
  '3' = 'CD4 Naive',
  '4' = 'CD14 Mono',
  '5' = 'CD8 Naive',
  '6' = 'DN T',
  '7' = 'NK CD56Dim',
  '8' = 'pre-B',
  '9' = 'CD16 Mono',
  '10' = 'pro-B',
  '11' = 'DC',
  '12' = 'NK CD56bright',
  '13' = 'pDC'
)

DefaultAssay(pbmc) <- 'peaks'
da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = "CD4 Naive",
  #ident.2 = "CD14 Mono",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)


open_cd4naive <- rownames(da_peaks[da_peaks$avg_logFC > 0.5, ])
open_cd14mono <- rownames(da_peaks[da_peaks$avg_logFC < -0.5, ])

closest_genes_cd4naive <- ClosestFeature(
  regions = open_cd4naive,
  annotation = gene.ranges,
  sep = c(':', '-')
)
closest_genes_cd14mono <- ClosestFeature(
  regions = open_cd14mono,
  annotation = gene.ranges,
  sep = c(':', '-')
)
rownames(da_peaks)[c(1,5)]

levels(pbmc) <- c("CD4 Naive","CD4 Memory","CD8 Naive","CD8 Effector","DN T","NK CD56bright","NK CD56Dim","pre-B",'pro-B',"pDC","DC","CD14 Mono",'CD16 Mono')

chr1:6896024-6929974

a<- c("chr12:69742121-69748015","chr2:85912298-85925978","chr11:60223225-60238234",
      "chr14:99635624-99737862","chr2:87011729-87035520","chr11:118175260-118186891",
      "chr12:6896024-6929975","chr5:140011313-140013260")

CoveragePlot(
  object = pbmc,
  region = a,
  sep = c(":", "-"),
  peaks = StringToGRanges(rownames(pbmc), sep = c(":", "-")),
  annotation = gene.ranges,
  extend.upstream = 2000,
  extend.downstream = 0,
  ncol = 8
)

# marker_index <- which(gene.ranges@ranges@NAMES == "ENSG00000156738")
# start<-gene.ranges@ranges@start[marker_index]
# end<-gene.ranges@ranges@width[marker_index]
# start
# start+end

# 
# lyz: ENSG00000090382:chr12:69742121-69748015
# gnly: ENSG00000115523: chr2:85912298-85925978
# ms4a1: ENSG00000156738:chr11:60223225-60238234
# bcl11b: ENSG00000127152: chr14:99635624-99737862
# cd8a: ENSG00000153563: chr2:87011729-87035520
# cd3e: ENSG00000198851: chr11:118175260-118186891
# cd4: ENSG00000010610: chr12:6896024-6929975
# cd14: ENSG00000170458: chr5:140011313-140013260


# save(pbmc,file='/Users/zhongyuanke/Desktop/r_para/pbmc.Rdata')
# save(gene.ranges,file='/Users/zhongyuanke/Desktop/r_para/geneRanges.Rdata')
# 
# load('/Users/zhongyuanke/Desktop/r_para/geneRanges.Rdata')
# load('/Users/zhongyuanke/Desktop/r_para/pbmc.Rdata')

