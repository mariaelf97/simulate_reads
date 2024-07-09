library(tidyverse)
library(data.table)
library(fuzzyjoin)
library(IRanges)

# who 2023 mutation catalog
who_2023 <- fread("mnt/tb_seqs/who_drug_resistance/amr-tb-who-2023.csv")
# tbprofiler barcodes
tb_profiler_barcodes <- fread("mnt/tb_seqs/lineage_barcodes/tbdb.barcode.bed")
# read who drug resistance mutation positions
mutation_positions <- fread("mnt/tb_seqs/who_drug_resistance/genomic_location.csv")
# merge specific mutation locations in each drug resistance genes
who_2023_variants_joined <- who_2023 %>% inner_join(mutation_positions,by = "variant")

# create a df containing all the genomic positions from who drug resistance 
# keep one copy per position
who_pos <- who_2023_variants_joined %>% select(position,variant) %>% distinct(position,.keep_all = TRUE)
# create intervals for each position, 600 bp amplicon length 
who_pos$lower_bound <- who_pos$position 
who_pos$upper_bound <- who_pos$position + 600
# create a df containing all the genomic positions from tb-profiler lineage typing
tb_prof_pos <- tb_profiler_barcodes %>% select(V2,V4) 
colnames(tb_prof_pos) <- c("position","variant")
tb_prof_pos$lower_bound <- tb_prof_pos$position 
tb_prof_pos$upper_bound <- tb_prof_pos$position + 600

#merge AMR and lineage SNPs dataframes
all_data <- rbind(tb_prof_pos,who_pos)

# merge the dataframe with itself to get the overlapping intervals 
interval_joins<-all_data%>%
  interval_inner_join(all_data,
                      by =c("lower_bound","upper_bound"))
# get counts of how many genomic position will be covered in each interval
interval_cnts <- interval_joins %>% group_by(position.x) %>% summarise(cnt = n()) %>%
  inner_join(interval_joins, by= "position.x")
# get a list of variants and their positions included in each interval
mutation_per_position <- aggregate(position.y ~ position.x, interval_joins, list)
variant_per_position <- aggregate(variant.y ~ position.x, interval_joins, list)

# Merge variants and their positions included in intervals with tb-profiler barcodes
# Find how many variants overlap in each interval
final_df <- tb_profiler_barcodes %>% inner_join(variant_per_position, by= c("V2" = "position.x"))%>%
  inner_join(mutation_per_position, by= c("V2" = "position.x"))%>%rowwise() %>% 
  mutate(overlapping_variant_cnt = length(unlist(position.y))) %>% ungroup() 
# For each lineage, keep only variants with the maximum number of variants covered in one amplicon
max_overlaps <- final_df %>% group_by(V4) %>% mutate(max_cnt_per_lineage = max(overlapping_variant_cnt)) %>%
  filter(max_cnt_per_lineage == overlapping_variant_cnt) 
# When subsampling, make sure we are not deleting any overlapping variants
sub_sample_lineages <- function(i){
  overlapping_elements <- unlist(max_overlaps[i,10][[1]])
  df <- final_df %>% filter(V2 %in% overlapping_elements)
  return(df)
}
# apply the function
all_variants <- mapply(sub_sample_lineages,1:1111,SIMPLIFY = FALSE)
# remove duplicate rows
all_variants_unzipped <- do.call(rbind,all_variants) %>% distinct(V2, .keep_all = TRUE)
# remove unnecessary variable
all_variants_unzipped$V1 <- NULL
# rename column names
colnames(all_variants_unzipped) <- c("interval_start","interval_end","lineage",
                                     "expected_base","type1","type2","RD",
                                     "overlapping_variants","overlapping_variant_positions",
                                     "overlapping_variant_count")

### Drug resistance ###
# which drug resistance barcodes will we miss if we include counts higher than the median count value?
final_AMR <- who_2023_variants_joined %>% distinct(position,.keep_all = TRUE) %>%
  inner_join(variant_per_position, by= c("position" = "position.x"))%>%
  inner_join(mutation_per_position, by= c("position" = "position.x"))%>%
  rowwise() %>% 
  mutate(overlapping_variant_cnt = length(unlist(position.y))) %>% 
  ungroup() %>%
  mutate(variant_size = abs(nchar(alternative_nucleotide) - nchar(reference_nucleotide)))%>%
  mutate(interval_lower_bound = position, interval_upper_bound = position + 600)
  

final_AMR_grouped <-final_AMR %>% mutate(variant_type = case_when(
  nchar(reference_nucleotide) == nchar(alternative_nucleotide) ~ "SNP",
  nchar(reference_nucleotide) > nchar(alternative_nucleotide) ~ "deletion",
  nchar(reference_nucleotide) < nchar(alternative_nucleotide) ~ "insertion"))

final_AMR_grouped %>% filter(variant_size < 100) %>% 
  group_by(drug,gene) %>% mutate(max_cnt_per_gene = max(overlapping_variant_cnt)) %>%
  filter(max_cnt_per_gene == overlapping_variant_cnt) %>% view()

binned_data <- bin_data(data.table(final_AMR), binCol="position", bins=unlist(bins[[1]])) 

# gene names included in 600bps amplicon scheme by quick's lab
deeplex_genes <- fread("mnt/tb_seqs/primers/tb-amr-panel/600/v1.0.0/work/freschi.gene.bed")


# get gene names for drug resistance in quick's lab scheme
deeplex_amr <- deeplex_genes %>% filter(!(grepl('freschi', V4))) 
# which genes are in deeplex not covered by who?
deeplex_amr %>% filter(!(V4 %in% who_2023_drug_res_mut$gene))
# which genes are in who not covered by deeplex?
who_2023_drug_res_mut %>% filter(!(gene %in% deeplex_amr$V4))
