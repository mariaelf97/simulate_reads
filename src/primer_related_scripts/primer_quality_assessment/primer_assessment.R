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
         
                             