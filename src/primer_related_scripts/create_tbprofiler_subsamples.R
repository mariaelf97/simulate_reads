library(data.table)
library(RcppHungarian)
library(tidyverse)
### uniformity - uniform stratified sampling ###
# code credit for uniformity : https://stackoverflow.com/users/9463489/jblood94
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

###### optimize number of amplicons ####
# Get the distance matrix - pairwise distance between all SNPs
d <- outer(tb_profiler_barcodes$V2, tb_profiler_barcodes$V2, \(x, y) abs(x - y))
# get only uppper triangle values, remove redundant distances
upper_tri <- upper.tri(d)
upper_values <- d[upper_tri]
indices <- which(upper_tri, arr.ind = TRUE)
# convert it to a long format df
long_df <- data.frame(
  row = indices[, "row"],
  col = indices[, "col"],
  value = upper_values
)
# filter out distant values, remove diagonal values
d_neighbours <- long_df %>% filter(value > 0 & value <300) 
# use row ids to merge distance values 
tb_profiler_barcodes$rowids <- as.integer(rownames(tb_profiler_barcodes))
rep <- d_neighbours %>% inner_join(tb_profiler_barcodes, by=c("row"="rowids")) %>%
  inner_join(tb_profiler_barcodes, by=c("col"="rowids"))%>% select(V2.x,V2.y,V4.x,V4.y)%>%
  mutate(var_distance =abs(V2.x - V2.y))
# find unique lineages included
uniqs <- unique(c(rep$V4.x,rep$V4.y))
# find lineages not included
non_inc <- tb_profiler_barcodes %>% filter(!(V4 %in% uniqs))
# create a list of lineages to be included in addition
to_be_included <- non_inc %>% group_by(V4) %>% slice_sample(n=1) %>% select(V2,V4)
# process the new values to be able to merge them to rep df
to_be_included$V2.y <- as.integer("")
to_be_included$V4.y <- NA
rep$V2.y <- as.integer(rep$V2.y)
colnames(to_be_included) <- c("V2.x","V4.x","V2.y","V4.y")
# merge two dfs
all_rep <- rbind(to_be_included,rep,deparse.level = 1)
# change column names
colnames(all_rep) <- c("var1","lin1","var2","lin2","pairwise_distance")
# aggregate mutations that are in proximity of each variant  
mutation_per_position <- aggregate(var2 ~ var1, all_rep, list) 
variant_per_position <- aggregate(lin2~ var1, all_rep, list)
paiwisedist_per_position <- aggregate(pairwise_distance~ var1, all_rep, list)
# merge all values ( non-aggregated and aggregated)
final_df <- all_rep %>% 
  left_join(mutation_per_position, by= "var1") %>%
  left_join(variant_per_position, by= "var1")%>%
  left_join(paiwisedist_per_position,by="var1")%>% 
  select(var1,lin1,lin2.y,var2.y,pairwise_distance.y)
# change column names
colnames(final_df) <- c("main_variant","main_variant_lineage", "lineages_included",
                        "variants_included","pairwise_dist_main_variant")
# change column data types to string
final_df$lineages_included <- as.character(final_df$lineages_included)
final_df$variants_included <- as.character(final_df$variants_included)
final_df$pairwise_dist_main_variant <- as.character(final_df$pairwise_dist_main_variant)
# separate variants included into different columns for easier operation
final_df <- final_df %>% separate(lineages_included,into=c("lin_inc1","lin_inc2"),sep=",") 
# get unique lineage names 
lineages <- unique(tb_profiler_barcodes$V4)
# create a df to append unique variant and lineages
result_df <- data.frame(matrix(ncol = 6 , nrow = 0))
# name the columns same as the final_df
colnames(result_df) <- colnames(final_df)
result_df$lin_inc1 <- as.character(result_df$lin_inc1)
result_df$lin_inc2 <- as.character(result_df$lin_inc2)
result_df$variants_included <- as.character(result_df$variants_included)
result_df$main_variant_lineage <- as.character(result_df$lineages_included)
result_df$pairwise_dist_main_variant <- as.character(result_df$pairwise_dist_main_variant)
# create a list for lineages to be removed
removed_lineages <- c()
# process values to remove redundant variants and ensure 1 lineage defining SNP is included
while(length(lineages) > 0) {
  for (x in 1:length(lineages)) {
    lineage <- lineages[x]
    # If the lineage already exist, skip the loop
    if (nrow(result_df %>% filter(main_variant_lineage == lineage | lin_inc1 == lineage | lin_inc2 == lineage)) != 0) {
      removed_lineages <- c(removed_lineages, lineage)
      lineages <- setdiff(lineages, removed_lineages)
    } else {
      record_to_keep <- final_df %>% 
       filter(main_variant_lineage == lineage | lin_inc1 == lineage | lin_inc2 == lineage) %>% head(1)
      result_df <- bind_rows(result_df, record_to_keep)
      removed_lineages <- c(removed_lineages, record_to_keep$main_variant_lineage, record_to_keep$lineages_included)
      lineages <- setdiff(lineages, removed_lineages)
    }
  }
}
# process the final output
result_df <- result_df %>% separate(variants_included, into=c("var2","var3"), sep = ",")
result_df <- result_df %>% separate(pairwise_dist_main_variant, into=c("dist1","dist2"), sep = ",")
result_df$var2 <- gsub('^c\\(|\\)$', '', result_df$var2)
result_df$var3 <- gsub('\\)', '', result_df$var3)
result_df$dist1 <- gsub('^c\\(|\\)$', '', result_df$dist1)
result_df$dist2 <- gsub('\\)', '', result_df$dist2)
colnames(result_df) <- c("var1","var1_lin","var2_lin","var3_lin","var2","var3","dist1","dist2")
# change empty values to NA
result_df[] <- lapply(result_df, function(x) {
  x[x == "" | is.null(x) | x=="NULL"] <- NA
  return(x)
})
# plot the distribution of SNPs
result_df %>% ggplot(aes(var1)) + geom_histogram(bins=200)
# rearrange columns
result_df <- result_df %>% arrange(var1)%>%
  select(var1,var2,var3,var1_lin,var2_lin,var3_lin,dist1,dist2)

