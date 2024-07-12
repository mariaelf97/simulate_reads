library(data.table)
library(tidyverse)
library(fuzzyjoin)
library(IRanges)

coll_primers <- fread("mnt/tb_seqs/primers/primer.bed")
colnames(coll_primers) <- c("chr","start","end","name","pool","strand","seq")
coll_primers <- coll_primers %>%
  separate(name, into = c("name","number","side"))

# remove repeated name variable not needed
coll_primers$name <- NULL
coll_primers_joined <- coll_primers %>% filter(side == "LEFT") %>%
  inner_join((coll_primers %>% filter(side == "RIGHT")), by = "number")

coll_primers_joined <- coll_primers_joined %>%
  mutate(lower_bound = start.x, upper_bound = end.y)

read_depth_files <- function(isolate) {
  path_to_file <- paste("mnt/tb_seqs/NICD_TB_seq_data/",isolate,
                        "/depths.tsv",sep = "")
  depth_file <- fread(path_to_file)
  depth_file$isolate <- isolate
  return(depth_file)
}

isolates <- fread("mnt/tb_seqs/NICD_TB_seq_data/samples.txt",header = FALSE)
depth_files <- mapply(read_depth_files,isolates$V1, SIMPLIFY = FALSE)
depth_files_unzipped <- do.call(rbind,depth_files)
depth_files_unzipped<- depth_files_unzipped %>% select(V2,V4,isolate)%>%
  mutate(lower_bound = V2, upper_bound = V2 )
overlaps <- coll_primers_joined %>% 
  interval_inner_join(depth_files_unzipped,
                      by = c("lower_bound","upper_bound")) %>% 
  filter(V4 > 10) %>% 
  group_by(number,lower_bound.x,upper_bound.x,isolate) %>% 
  summarise(mean_depth = mean(V4))
colnames(overlaps) <- c("amplicon_number","amplicon_start","amplicon_end"
                        ,"sample", "mean_coverage_depth")
overlaps$amplicon_number <- as.numeric(overlaps$amplicon_number)
ggplot() +
  geom_line(depth_files_unzipped,
            mapping =aes(x=V2,y=V4, color = "red",linewidth =0.00001,
                         , alpha = 0.15))+
  geom_vline(coll_primers_joined,mapping = aes(xintercept = start.x,
                                                          color="green"))+
  theme_minimal()+ facet_wrap(.~isolate,scales="free_y") +ylab("Coverage_depth")+
  xlab("Genomic position")
