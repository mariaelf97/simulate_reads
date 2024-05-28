library(tidyverse)
library(data.table)


results<-read.table("mnt/tb_seqs/seq_simulation/amplicon_read_simulations/string_match_wgsim/primer_v2/pure_samples/aggregated_result.tsv", fill = TRUE, sep = "\t", h=T)
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
summary(results_comb_names)
results_comb_names %>% write_csv("mnt/tb_seqs/seq_simulation/amplicon_read_simulations/amplicon_alignment_wgsim/primer_v2/pure_samples/demix_output_eval.csv")
