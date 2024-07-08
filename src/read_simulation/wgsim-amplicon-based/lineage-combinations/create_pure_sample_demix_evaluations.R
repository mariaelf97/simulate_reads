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


results<-read.table("mnt/tb_seqs/seq_simulation/amplicon_read_simulations/string_match_wgsim/primer_v3/pure_samples/aggregated_result.tsv", fill = TRUE, sep = "\t", h=T)
results<-as.data.frame(sapply(results, function(x) str_replace_all(x, "[',()\\]\\[]", ""))) # Removed the unwanted character: [], () and commas
results<-as.data.frame(sapply(results, function(x) trimws(gsub("\\s+", " ", x)))) # Removed double spaces

lineages <- fread("mnt/tb_seqs/assemblies/long_read_assemblies/pacbio-RSII/lineages.tsv")

results_comb_names <- results %>% separate(lineages, into = c("obs_lin"), sep = " ")
results_comb_names <- results_comb_names %>% separate(abundances, into = c("obs_abun"), sep = " ")
results_comb_names <-results_comb_names[,c(1,3,4)]

results_comb_names <- results_comb_names %>% 
  separate(X, into = as.character(1:4))

results_comb_names$isolate <- paste(results_comb_names$`1`,results_comb_names$`2`, sep ="_")
results_comb_names$isolate <- paste(results_comb_names$isolate,".","1",sep="")
results_comb_names <- results_comb_names[,-c(1,2,3)]

results_comb_names <- lineages %>% select(Isolate,freschi2020) %>%
  inner_join(results_comb_names, by =c("Isolate" ="isolate")) 
colnames(results_comb_names) <- c("isolate","exp_lin","read_cnt","obs_lin","obs_abun")

results_comb_names <- results_comb_names %>% relocate(isolate,exp_lin,obs_lin,read_cnt,obs_abun)

results_comb_names <- results_comb_names %>% mutate(exp_abun = 1) %>%
  mutate(resiual = as.numeric(exp_abun) - as.numeric(obs_abun)) %>% 
  mutate(exp_abun = round(exp_abun,2)) %>%
  mutate(obs_abun = round(as.numeric(obs_abun),2))

results_comb_names <- results_comb_names %>% mutate(exp_lin = paste("lineage",exp_lin,sep = ""))
results_merged <- results_comb_names %>% full_join(json_files_unzipped, by = "isolate")%>%
  full_join(json_files_unzipped_wgs, by = "isolate")
colnames(results_merged) <- c("isolate","tb_profiler_exp_lin_whole_genome","freyja_obs_lin","read_cnt",
                              "freyja_obs_abun","exp_abun","residual","tb_profiler_amplicon_lineage",
                              "tb_profiler_amplicon_abun","tb_profiler_wgs_lineage",
                              "tb_profiler_wgs_abun")
results_merged %>% write_csv("Downloads/demix_eval_rev.csv")
