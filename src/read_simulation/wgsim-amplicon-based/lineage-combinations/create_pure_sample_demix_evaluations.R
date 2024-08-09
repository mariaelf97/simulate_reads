library(tidyverse)
library(data.table)
library(rjson)

# read tbprofiler json file
read_tbprofiler_json <- function(isolate){
  tb_profiler_json <- fromJSON(file = paste("mnt/tb_seqs/seq_simulation/",
                                            "amplicon_read_simulations/string_match_wgsim/",
                                            "primer_v3/pure_samples/",
                                            isolate,"/reverse_complement/25000/",
                                            "tb-profiler/results/tbprofiler.results.json",
                                            sep = ""))
  if(length(tb_profiler_json$lineage) >0)
  {tb_profiler_df <- as.data.frame(tb_profiler_json$lineage)
  tb_profiler_df$isolate <- isolate 
  tb_profiler_df <- tb_profiler_df %>% select(lineage,support.target_allele_percent)
  return(tb_profiler_df)}
  
}
read_tbprofiler_json_wgs <- function(isolate){
  tb_profiler_json <- fromJSON(file = paste("mnt/tb_seqs/seq_simulation/amplicon_read_simulations/string_match_wgsim/primer_v3/",isolate,
                                            "/whole_genome/25000/tb-profiler/results/tbprofiler.results.json",
                                            sep = ""))
  if(length(tb_profiler_json$lineage) >0)
  {tb_profiler_df <- as.data.frame(tb_profiler_json$lineage)
  tb_profiler_df$isolate <- isolate 
  tb_profiler_df <- tb_profiler_df %>% select(lineage,support.target_allele_percent)
  return(tb_profiler_df)}
  
}
isolate <- fread("mnt/tb_seqs/assemblies/long_read_assemblies/pacbio-RSII/isolates.txt", header = FALSE)
json_files <- mapply(read_tbprofiler_json, isolate$V1, SIMPLIFY = FALSE)
json_files_unzipped <- do.call(rbind, json_files)
json_files_unzipped$isolate <- rownames(json_files_unzipped)

json_files_wgs <- mapply(read_tbprofiler_json_wgs, isolate$V1, SIMPLIFY = FALSE)
json_files_unzipped_wgs <- do.call(rbind, json_files_wgs)
json_files_unzipped_wgs$isolate <- rownames(json_files_unzipped_wgs)


results<-read.table("mnt/measles_seqs/seq_simulation/combinations/aggregated_result.tsv", fill = TRUE, sep = "\t", h=T)
results<-as.data.frame(sapply(results, function(x) str_replace_all(x, "[',()\\]\\[]", ""))) # Removed the unwanted character: [], () and commas
results<-as.data.frame(sapply(results, function(x) trimws(gsub("\\s+", " ", x)))) # Removed double spaces

lineages <- fread("mnt/measles_seqs/assemblies/lineages.tsv")

results_comb_names <- results %>% separate(lineages, into = c("obs_lin1","obs_lin2"), sep = " ")
results_comb_names <- results_comb_names %>% separate(abundances, into = c("obs_abun1","obs_abun2"), sep = " ")
results_comb_names <-results_comb_names[,c(1,3,4,5,6)]

results_comb_names <- results_comb_names %>% 
  separate(X, into = as.character(1:6),sep="_")

results_comb_names$isolate1 <- paste(results_comb_names$`1`,results_comb_names$`2`, sep ="_")
results_comb_names <- results_comb_names[,-c(1,2)]

results_comb_names$isolate2 <- paste(results_comb_names$`4`,results_comb_names$`5`, sep ="_")
results_comb_names <- results_comb_names[,-c(2,3)]

results_comb_names$exp_abun1 <- as.numeric(str_remove(results_comb_names$'6', ".vcf"))

results_comb_names<-results_comb_names[,-2]
colnames(results_comb_names) <- c("exp_abund2","obs_lin1","obs_lin2","obs_abun1","obs_abun2","isolate1","isolate2","exp_abund1")

results_comb_names <- results_comb_names %>% mutate(real_obs_abun1 = case_when(obs_abun1 > obs_abun2~ obs_abun1,
                        obs_abun2 > obs_abun1 ~ obs_abun2, obs_abun1 == obs_abun2~ obs_abun1))
results_comb_names <- results_comb_names %>% mutate(real_obs_abun2 = case_when(obs_abun1 < obs_abun2~ obs_abun1,
                                                                               obs_abun2 < obs_abun1 ~ obs_abun2, obs_abun1 == obs_abun2~ obs_abun1))
results_comb_names <- results_comb_names %>% mutate(exp_abund1 = as.numeric(exp_abund1))
results_comb_names <- results_comb_names %>% mutate(exp_abund2 = as.numeric(exp_abund2))
results_comb_names <- results_comb_names %>% mutate(real_exp_abun1 = case_when(exp_abund1 > exp_abund2~ exp_abund1,
                                                                               exp_abund2 > exp_abund1 ~ exp_abund2, exp_abund1 == exp_abund2~ exp_abund1))
results_comb_names <- results_comb_names %>% mutate(real_exp_abun2 = case_when(exp_abund1 < exp_abund2~ exp_abund1,
                                                                               exp_abund2 < exp_abund1 ~ exp_abund2, exp_abund1 == exp_abund2~ exp_abund1))
results_comb_names <- results_comb_names %>% mutate(res1 = abs(as.numeric(real_obs_abun1) - as.numeric(real_exp_abun1)))
results_comb_names <- results_comb_names %>% mutate(res2 = abs(as.numeric(real_obs_abun2) - as.numeric(real_exp_abun2)))

results_comb_names %>% ggplot(aes(as.numeric(real_obs_abun2),as.numeric(real_exp_abun2))) + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")+
  theme_minimal() + xlab("Lineage observed abundance") + ylab("Lineage expected abundance")
rss <- sum((as.numeric(results_comb_names$real_obs_abun1) - as.numeric(results_comb_names$real_exp_abun1)) ^ 2,na.rm = TRUE)  ## residual sum of squares
tss <- sum((results_comb_names$real_exp_abun1 - mean(results_comb_names$real_exp_abun1)) ^ 2,na.rm = TRUE)  ## total sum of squares
rsq <- 1 - rss/tss

results_comb_names <- lineages %>% select(seqName,clade) %>%
  inner_join(results_comb_names, by =c("seqName" ="isolate1")) 
colnames(results_comb_names)[1:2] <- c("isolate1","exp_lin1")

results_comb_names <- lineages %>% select(seqName,clade) %>%
  inner_join(results_comb_names, by =c("seqName" ="isolate2")) 
colnames(results_comb_names)[1:2] <- c("isolate2","exp_lin2")
results_comb_names <- results_comb_names %>% select(isolate1,isolate2,exp_lin1,exp_lin2,obs_lin1,obs_lin2,real_exp_abun1,real_exp_abun2,real_obs_abun1,real_obs_abun2,res1,res2)
colnames(results_comb_names)[7:10] <-c("exp_abun1","exp_abun2","obs_abun1","obs_abun2")
