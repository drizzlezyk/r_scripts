library(Seurat)
base_path <- '/Users/zhongyuanke/data/'
batch_path <- paste(base_path, 'merge_result/293t_jurkat_batch.csv', sep = "")
origi_csv <- paste(base_path, 'csv/origi.csv', sep = "") 
mse_csv <- paste(base_path, 'csv/mse.csv', sep = "") 
dca_csv <- paste(base_path, 'csv/dca.csv', sep = "") 
scxx_csv <- paste(base_path, 'csv/scxx.csv', sep = "") 
scan_csv <- paste(base_path, 'csv/scan.csv', sep = "") 

x_origi <-read.csv(origi_csv, header =FALSE)
x_mse <-read.csv(mse_csv, header =FALSE)
x_dca <-read.csv(dca_csv, header =FALSE)
x_scxx <-read.csv(scxx_csv, header =FALSE)
x_scan <-read.csv(scan_csv, header =FALSE)

x_origi<-t(x_origi)
x_mse<-t(x_mse)
x_dca<-t(x_dca)
x_scxx<-t(x_scxx)
x_scan<-t(x_scan)

x_origi<-object(x_origi)
count = 1
data <- ReadH5AD(file_dca, verbose = TRUE)
arr1 <- array(data = NA, dim=1)
for (i in 1:3){
  for(j in 1:test[1,i]+1)
  {
    arr1[count]<- i
    count = count+1
  }
}
print(arr1)
adata_dca = ReadH5AD(file_dca)
x_dca = adata_dca.obsm['mid']

MixingMetric <- function(
  object,
  grouping.var,
  reduction = "pca",
  dims = 1:2,
  k = 5,
  max.k = 300,
  eps = 0,
  verbose = TRUE
) {
  my.sapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pbsapply,
    no = future_sapply
  )
  nn <- nn2(
    data = object,
    k = max.k,
    eps = eps
  )
  group.info <- object[[grouping.var, drop = TRUE]]
  groups <- unique(x = group.info)
  mixing <- my.sapply(
    X = 1:ncol(x = object),
    FUN = function(x) {
      sapply(X = groups, FUN = function(y) {
        which(x = group.info[nn$nn.idx[x, ]] == y)[k]
      })
    }
  )
  mixing[is.na(x = mixing)] <- max.k
  mixing <- apply(
    X = mixing,
    MARGIN = 2,
    FUN = median
  )
  return(mixing)
}
m=MixingMetric(
  x_origi,
  arr1,
  reduction = "pca",
  dims = 2,
  k = 5,
  max.k = 300,
  eps = 0,
  verbose = TRUE
)
