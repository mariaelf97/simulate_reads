library(tidyverse)
library(data.table)


results<-read.table("mnt/tb_seqs/seq_simulation/amplicon_read_simulations/amplicon_alignment_wgsim/primer_v2/mixed_samples_diff_cov/aggregated_result.tsv", fill = TRUE, sep = "\t", h=T)
results<-as.data.frame(sapply(results, function(x) str_replace_all(x, "[',()\\]\\[]", ""))) # Removed the unwanted character: [], () and commas
results<-as.data.frame(sapply(results, function(x) trimws(gsub("\\s+", " ", x)))) # Removed double spaces

lineages <- fread("mnt/tb_seqs/assemblies/long_read_assemblies/pacbio-RSII/lineages.tsv")


results_comb_names <- results %>% separate(lineages, into = c("obs_lin1","obs_lin2"), sep = " ")
results_comb_names <- results_comb_names %>% separate(abundances, into = c("obs_abun1","obs_abun2"), sep = " ")
results_comb_names <-results_comb_names[,c(1,3,4,5,6)]

results_comb_names <- results_comb_names %>% 
  separate(X, into = as.character(1:8))

results_comb_names$isolate1 <- paste(results_comb_names$`1`,results_comb_names$`2`, sep ="_")
results_comb_names$isolate1 <- paste(results_comb_names$isolate1,".","1",sep="")
results_comb_names$isolate2 <- paste(results_comb_names$`5`,results_comb_names$`6`, sep ="_")
results_comb_names$isolate2 <- paste(results_comb_names$isolate2,".","1",sep="")

results_comb_names <- results_comb_names[,-c(1,2,3,5,6,7)]

colnames(results_comb_names) <- c("exp_abun2","exp_abun1","obs_lin1","obs_lin2","obs_abun1","obs_abun2","isolate1","isolate2")

results_comb_names$obs_abun1 <- as.numeric(results_comb_names$obs_abun1)
results_comb_names$exp_abun1 <- as.numeric(results_comb_names$exp_abun1)
results_comb_names$obs_abun2 <- as.numeric(results_comb_names$obs_abun2)
results_comb_names$exp_abun2 <- as.numeric(results_comb_names$exp_abun2)


results_comb_names <- results_comb_names%>%
  mutate(exp_abun1_ratio = exp_abun1 / (exp_abun1 + exp_abun2))
results_comb_names <- results_comb_names%>%
  mutate(exp_abun2_ratio = exp_abun2 / (exp_abun1 + exp_abun2))

results_comb_names$exp_abun1 <- NULL
results_comb_names$exp_abun2 <- NULL


colnames(results_comb_names) <- c("obs_lin1","obs_lin2","obs_abun1",
                                  "obs_abun2","isolate1","isolate2",
                                  "exp_abun1","exp_abun2")



results_comb_names %>% ggplot() +
  geom_point(mapping=aes(exp_abun1,obs_abun1,color="green"))+
  geom_point(mapping=aes(exp_abun2,obs_abun2,color="red")) 


results_comb_names <- lineages %>% select(Isolate,freschi2020) %>%
  inner_join(results_comb_names, by =c("Isolate" ="isolate1")) 
colnames(results_comb_names) <- c("isolate1","exp_lin1","obs_lin1","obs_lin2","obs_abun1",
                                 "obs_abun2","isolate2",
                                 "exp_abun1","exp_abun2")
results_comb_names <- lineages %>% select(Isolate,freschi2020) %>%
  inner_join(results_comb_names, by =c("Isolate" ="isolate2")) 
colnames(results_comb_names) <- c("isolate2","exp_lin1","isolate1",
                                  "exp_lin2","obs_lin1","obs_lin2","obs_abun1",
                                  "obs_abun2", "exp_abun1","exp_abun2")
results_comb_names$exp_lin2 <- paste("lineage",results_comb_names$exp_lin2,sep = "")
results_comb_names$exp_lin1 <- paste("lineage",results_comb_names$exp_lin1,sep = "")
results_comb_names <- results_comb_names %>% 
  mutate(dev1 = abs(obs_abun1 - exp_abun1)) %>%
  mutate(dev2 = abs(obs_abun2 - exp_abun2))
