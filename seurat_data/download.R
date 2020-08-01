
#### Download 
urls<-c(
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126074&format=file&file=GSE126074%5FP0%5FBrainCortex%5FSNAREseq%5FcDNA%2Ebarcodes%2Etsv%2Egz",
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126074&format=file&file=GSE126074%5FP0%5FBrainCortex%5FSNAREseq%5FcDNA%2Ecounts%2Emtx%2Egz",
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126074&format=file&file=GSE126074%5FP0%5FBrainCortex%5FSNAREseq%5FcDNA%2Egenes%2Etsv%2Egz",
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126074&format=file&file=GSE126074%5FP0%5FBrainCortex%5FSNAREseq%5Fchromatin%2Ebarcodes%2Etsv%2Egz",
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126074&format=file&file=GSE126074%5FP0%5FBrainCortex%5FSNAREseq%5Fchromatin%2Ecounts%2Emtx%2Egz",
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126074&format=file&file=GSE126074%5FP0%5FBrainCortex%5FSNAREseq%5Fchromatin%2Epeaks%2Etsv%2Egz"
)
fnames<-c(
  "GSE126074_P0_BrainCortex_SNAREseq_cDNA/barcodes.tsv.gz",
  "GSE126074_P0_BrainCortex_SNAREseq_cDNA/matrix.mtx.gz",
  "GSE126074_P0_BrainCortex_SNAREseq_cDNA/features.tsv.gz",
  "GSE126074_P0_BrainCortex_SNAREseq_chromatin/barcodes.tsv.gz",
  "GSE126074_P0_BrainCortex_SNAREseq_chromatin/matrix.mtx.gz",
  "GSE126074_P0_BrainCortex_SNAREseq_chromatin/features.tsv.gz"
)
for(i in 1:6){
 	download.file(url = urls[i],destfile = paste0("./data/atacseq/snareseq/",fnames[i]))
}
mm10gtf="http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.refGene.gtf.gz"
download.file(mm10gtf,"./data/atacseq/snareseq/mm10_refgene.gtf.gz")
