library(data.table)
library(tidyverse)

primer_bed <- fread("mnt/tb_seqs/primers/primer_v2.bed")
primer_bed <- primer_bed %>% filter(V6 == "+") %>% mutate(value = 20)
read_variant_files <- function(isolate) {
  path_to_file <- paste("mnt/tb_seqs/structural-variants/",isolate,
                        "/OUT.Assemblytics_structural_variants.bed", sep = "")
  variant_file <- fread(path_to_file)
  variant_file$isolate <- isolate
  return(variant_file)
}
isolate <- fread("mnt/tb_seqs/assemblies/long_read_assemblies/pacbio-RSII/isolates.txt", header = FALSE)
variant_files <- mapply(read_variant_files, isolate$V1, SIMPLIFY = FALSE)
variant_files_unzipped <- do.call(rbind, variant_files)

failed_isolates <- c("GCA_014884625.1","GCA_014899425.1","GCA_014899745.1",
                     "GCA_014899885.1","GCA_014899905.1", "GCA_014899945.1",
                     "GCA_014900095.1","GCA_014900195.1")

variant_files_unzipped <- variant_files_unzipped %>%
  mutate(isolate_type = ifelse(isolate %in% failed_isolates, "failed_isolates","normal_isolate"))

variant_files_unzipped %>% group_by(isolate_type, type) %>%
  summarise(mean_sv_length = mean(size), median_sv_length = median(size),
            sd_sv_length = sd(size), sv_counts = n(), normalized_mean = mean(size)/n())%>% view()


variant_files_unzipped %>%
  filter(type %in% c("Insertion","Deletion"))%>%
  filter(size > 2000) %>%
           ggplot() + geom_histogram(aes(ref_start, fill = type),bins = 100) + facet_wrap(.~isolate_type)
         