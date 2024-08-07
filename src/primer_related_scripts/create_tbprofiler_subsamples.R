library(data.table)
library(RcppHungarian)
# code credit : https://stackoverflow.com/users/9463489/jblood94
# read data
tb_profiler_barcodes <- fread("mnt/tb_seqs/lineage_barcodes/tbdb.barcode.bed")
# get the ideal spacing
r <- range(tb_profiler_barcodes$V2)
n <- length(unique(tb_profiler_barcodes$V4))
ideal <- seq(0.5, n - 0.5)*diff(r)/n + r[1] # ideal "even" spacing
# For each value in V4 (lineage), get the best candidate for each ideal location 
d <- outer(tb_profiler_barcodes$V2, ideal, \(x, y) abs(x - y))
idx <- as.matrix(
  cbind(data.table(lineage = tb_profiler_barcodes$V4)[,ID := .I], d)[
    ,lapply(.SD, \(x) ID[which.min(x)]), lineage, .SDcols = 3:(n + 2)
  ][,lineage := NULL]
)
# Get the distance for each row index.
mindists <- idx
mindists[] <- d[cbind(c(idx), c(col(idx)))]
#Solve the assignment problem and take the samples
samples <- tb_profiler_barcodes[
  idx[RcppHungarian::HungarianSolver(mindists)$pairs],
]

