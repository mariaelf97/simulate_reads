library(data.table)
library(tidyverse)
# function to read amplicon length files
  read_variant_files <- function(isolate) {
    path_to_file <- paste("mnt/tb_seqs/seq_simulation/amplicons/string_match_approach/primer_V3/",isolate,
                          "/reverse_complement/amplicon_stats.csv", sep = "")
    variant_file <- fread(path_to_file)
    variant_file$isolate <- isolate
    return(variant_file)
  }
# read isolate names and lineages
isolate <- fread("mnt/tb_seqs/assemblies/long_read_assemblies/pacbio-RSII//isolates.txt", header = FALSE)
lineages <- fread("mnt/tb_seqs/assemblies/long_read_assemblies/pacbio-RSII/lineages.tsv")
# read amplicon length files
variant_files <- mapply(read_variant_files, isolate$V1, SIMPLIFY = FALSE)
variant_files_unzipped <- do.call(rbind, variant_files)
primer_24 <- variant_files_unzipped %>% filter(amplicon_number == 24)
# import primer locations in H37Rv
primer_file_H37Rv <- fread("mnt/tb_seqs/primers/primer_v3.bed")
colnames(primer_file_H37Rv) <- c("chr","start","end","name","pool","strand","seq")
# divide amplicon name to name, number and side (e.g. left or right)
primer_file_H37Rv <- primer_file_H37Rv %>%
  separate(name, into = c("name","number","side"))
# remove repeated name variable not needed
primer_file_H37Rv$name <- NULL

primer_file_H37Rv$name <- "H37Rv_primer_loc"
H37Rv_24 <- primer_file_H37Rv%>% filter(number =="24")
ggplot() +
  geom_dumbbell(primer_24,
                mapping = aes(
                  x = primer_start,
                  xend = primer_end,
                  y = isolate
                ),
                colour = "pink",
                colour_xend = "#049000",
                size = 0.9
  )+   geom_dumbbell(H37Rv_24,
                     mapping = aes(
                       x = start,
                       xend = end,
                       y= name
                     ),
                     colour = "yellow",
                     colour_xend = "red",
                     size = 0.8
  )+
  
  theme_minimal()
# find primers with more than one match
  variant_files_unzipped %>% group_by(isolate,amplicon_number)%>%
    mutate(counts=n()) %>% view()
# find counts of amplicon dropout for each isolate
  variant_files_unzipped%>%
    inner_join(lineages, by=c("isolate"="Isolate")) %>%
    mutate(isolate_lineage = paste(isolate,coll2014,sep = "_"))%>%
    group_by(isolate_lineage) %>%
    summarise(ampl_dropout_cnt= sum(amplicon_length == "0"))%>% view()
# create amplicon length matrix
variant_files_unzipped%>%
  inner_join(lineages, by=c("isolate"="Isolate")) %>%
  mutate(isolate_lineage = paste(isolate,freschi2020,sep = "_"))%>%
           select(isolate_lineage,amplicon_number,amplicon_length)%>%
           reshape(idvar = "isolate_lineage", timevar = "amplicon_number", direction = "wide")%>%
  write_csv("amplicon_length_v3_inv.csv")
         
                             