library(data.table)
library(tidyverse)
library(fuzzyjoin)
library(IRanges)

# import primer locations in H37Rv
primer_file_H37Rv <- fread("mnt/tb_seqs/primers/primer_v3.bed")
colnames(primer_file_H37Rv) <- c("chr","start","end","name","pool","strand","seq")
# divide amplicon name to name, number and side (e.g. left or right)
primer_file_H37Rv <- primer_file_H37Rv %>%
  separate(name, into = c("name","number","side"))
# remove repeated name variable not needed
primer_file_H37Rv$name <- NULL

# import primer length files 
read_primer_files <- function(isolate) {
  path_to_file <- paste("mnt/tb_seqs/seq_simulation/amplicons/",
  "string_match_approach/primer_V3/",
                        isolate,"/amplicon_stats.csv", sep = "")
  primer_file <- fread(path_to_file)
  primer_file$isolate <- isolate
  return(primer_file)
}
# read SV files
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

primer_assess_files <-mapply(read_primer_files, isolate$V1, SIMPLIFY = FALSE)
primer_assess_unzipped <- do.call(rbind, primer_assess_files)

# defining isolates with high number of amplicon dropouts 
failed_isolates <- c("GCA_014884625.1","GCA_014899425.1","GCA_014899745.1",
                     "GCA_014899885.1","GCA_014899905.1", "GCA_014899945.1",
                     "GCA_014900095.1","GCA_014900195.1")
# Inner join the primer dataset with itself to get letf and right primer on one row
amplicon_coords <- primer_file_H37Rv%>% filter(side == "LEFT") %>% 
  inner_join(primer_file_H37Rv[primer_file_H37Rv$side == "RIGHT",], by = "number") %>%
  select(start.x,end.y,number)
colnames(amplicon_coords)<- c("start","end","number")
# select and rename columns in SV dataset
sv_selected <- variant_files_unzipped %>% select(ref_start,ref_stop,size,isolate,type)
colnames(sv_selected) <- c("start","end","size","isolate","type")
# join by SV and primer intervals to see if amplicons fall within any SVs
# inner join with amplicon length dataset to check amplicon dropout
sv_selected %>% 
  interval_inner_join(amplicon_coords,
                      by = c("start","end")) %>% 
  mutate(isolate_type = ifelse(isolate %in% failed_isolates,
                               "failed_isolates","normal_isolate")) %>% 
  mutate(number = as.numeric(number))%>%
  inner_join(primer_assess_unzipped, by= c("isolate"="isolate","number"="amplicon_number"))%>%
  select(-c(primer_start,primer_end,primer_seq_y))%>%
  rename("start.x" = "sv.start", "end.x" = "sv.end","size"="sv_size",
         "type" ="sv_type","start.y" = "start_amplicon","end.y"="end_amplicon",
         "number" ="amplicon_number")%>% view()

## look at the distribution of Svs across the genome and wether or not 
# they overlap an amplicon
primer_file_H37Rv$name <- "H37Rv_primer_loc"
library(ggalt)
png("image.png", width = 17000, height = 12000, res= 1100)
p <-  ggplot() +  geom_vline(xintercept = 10)+
  geom_dumbbell(variant_files_unzipped,
    mapping = aes(
      x = ref_start,
      xend = ref_stop,
      y = isolate
    ),
    colour = "pink",
    colour_xend = "#049000",
    size = 0.9
  )+   geom_dumbbell(primer_file_H37Rv,
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
print(p)

dev.off()
