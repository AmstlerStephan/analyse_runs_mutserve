library(plyr)
library(tidyr)
library(readr)
library(jsonlite)
library(tidyverse)
library(argparser)

parser <- arg_parser("Commandline parser")
parser <- add_argument(
  parser,
  "--run_path",
  help = "Run directory"
)
parser <- add_argument(
  parser,
  "--nanostat_summary",
  help = "Parsed nanostat summary of the run"
)
parser <- add_argument(
  parser,
  "--mutserve_summary",
  help = "Raw mutserve summary file"
)
parser <- add_argument(
  parser,
  "--ngs_data",
  help = "File containing the NGS reference data"
)

argv <- parse_args(parser)
run_path <- argv$run_path
nanostat_summary <- argv$nanostat_summary
mutserve_summary <- argv$mutserve_summary
ngs_data <- argv$ngs_data

### define parameters 
STR_start <- 2472
STR_end <- 2505
UMI_cutoff <- 0.003

sample_sets <- c("AK|SAPHIR")
SAPHIR_samples <- c("4612|4901|4451|4624|4864|5248|5538|5400")
AK_samples <- c("AK03|AK07|AK14|AK17|AK33")
# 
# run_path <- "~/UMI_LPA_KIV2/run7/"
# nanostat_summary <- "~/post_pipeline_analysis/QC/Nanostat_parsed_merged/run7/run7_0_0.tsv"
# mutserve_summary <- "~/UMI_LPA_KIV2/run7/mutserve/run7_summary_mutserve.txt"
# bed_file <- "~/UMI_LPA_KIV2/run7/ont_pl/barcode01/targets.bed"
# ngs_data <- "~/UMI_LPA_KIV2/data_ngs/data_ngs/20221122_NGS_reference_data_SAPHIR_mutation_classification.csv"

run <- str_extract(run_path, "run\\d*_*[[:alnum:]]*")

mutserve_summary <-
  read_tsv(mutserve_summary, na = c('', 'NA', '-'))

barcodes <-
  read_tsv(nanostat_summary) %>% 
  select(run:Sample, number_of_reads, mean_qual) %>% 
  mutate(sample = str_sub(Sample, end = -6),
         fragment = str_sub(Sample, start = -4))

NGS <- read_csv(ngs_data) %>% 
  filter(pos < STR_start | pos > STR_end) %>% 
  mutate(fragment = as.character(fragment))

### filter mutserve data 
### filter for full conversions (called variant is not the reference AND has no minor variant level OR minor variant level is below a certain threshold)
mutserve_raw_full_conversions <- mutserve_summary %>% 
  filter(REF != `TOP-REV` & (is.na(`MINOR-REV`) | `MINOR-REV-PERCENT` < UMI_cutoff)) %>% 
  mutate(Variant_level_UMI = 1, 
         Variant_UMI = `TOP-REV`)

### drop all NA values, where no minor variant was found (Either full conversion or no variant at that position)
### Variant level can be over 50% -> take Percentage and variant accordingly
mutserve_raw_variants <- mutserve_summary %>% 
  drop_na(`MINOR-REV`) %>% 
  mutate(Variant_level_UMI = ifelse((`REF` == `MINOR-REV`), `TOP-REV-PERCENT`, `MINOR-REV-PERCENT`),
         Variant_UMI = ifelse((`REF` == `MINOR-REV`), `TOP-REV`, `MINOR-REV`))

mutserve_combined <- bind_rows(mutserve_raw_full_conversions, mutserve_raw_variants) %>% 
  mutate( barcode = str_extract(SAMPLE, "barcode\\d\\d")) %>% 
  dplyr::rename(pos = POS)

### Join Barcodes and mutserve data

UMI <- mutserve_combined %>%
  inner_join(barcodes, by = c('barcode'))

UMI_readable <- UMI %>%
  filter(grepl('A_B', sample)) %>%
  separate(sample,
           c(NA, NA, 'Percent_A', 'Percent_B'),
           sep = '_',
           remove = FALSE) %>%
  mutate(
    Percent_A = as.numeric(Percent_A) / 10 ,
    Percent_B = as.numeric(Percent_B) / 10 ,
    Sample_readable_PL = paste(Percent_A, Percent_B, sep = ':')
  ) %>% 
  select(barcode, Sample_readable_PL) %>% 
  right_join(UMI) %>% 
  mutate( sample_readable = ifelse(is.na(Sample_readable_PL), sample, Sample_readable_PL))

UMI_readable_filtered <- UMI_readable %>% 
  filter(Variant_level_UMI > UMI_cutoff) %>%
  filter(pos < STR_start | pos > STR_end)

UMI_Samples <- UMI_readable %>% 
  filter(grepl(sample_sets, sample))

### getting a dataframge with AKs according fragments included in the sequencing run 
UMI_Samples_groups <- UMI_Samples %>% group_by(sample_readable, fragment) %>% summarize()

### the groups are used to filter the NGS data before joining both dataframes
### exclude all samples that are not covered with the UMI run

NGS_Samples <- NGS %>%
  filter(sample %in% UMI_Samples_groups$sample_readable & fragment %in% UMI_Samples_groups$fragment)

NGS_UMI_Samples <- NGS_Samples %>% 
  full_join(UMI_Samples,
            by = c('fragment', 'sample', 'pos')) %>%
  dplyr::rename(
    position = pos,
    ref_UMI = REF,
    variant_NGS = variant,
    variant_level_NGS = variant_level,
    ref_NGS = ref, 
    num_of_consensus_sequences = `COV-TOTAL`,
    variant_UMI = Variant_UMI,
    variant_level_UMI = Variant_level_UMI 
  ) %>%
  select(
    sample,
    fragment,
    run,
    position,
    ref_UMI,
    variant_UMI,
    variant_level_UMI,
    ref_NGS,
    variant_NGS,
    variant_level_NGS, 
    run, 
    number_of_reads,
    num_of_consensus_sequences, 
    mean_qual
  ) %>%
  mutate(
    variant_level_UMI = coalesce(variant_level_UMI, 0),
    variant_level_NGS = coalesce(variant_level_NGS, 0),
    variance_level_absolute_difference = variant_level_NGS - variant_level_UMI,
    # Variance_level_relative_difference = (Variant_level_NGS / Variant_level_UMI - 1)
  )

NGS_UMI_Samples_filtered <- NGS_UMI_Samples %>% 
  filter(!is.na(variant_NGS) | (variant_level_UMI > UMI_cutoff | variant_level_UMI == 0)) %>%
  filter(position < STR_start | position > STR_end) %>%
  filter(variant_NGS != 'D')

write.csv(UMI_readable, 'UMI_sequencing_mutserve.csv')
write.csv(UMI_readable_filtered, 'UMI_sequencing_filtered_mutserve.csv')
#write.csv(mutserve_combined, paste(run_path, 'mutserve_UMI_combined.csv', sep = ''))
write.csv(NGS_UMI_Samples, paste(run_path, 'NGS_UMI_samples.csv', sep = ''))
write.csv(NGS_UMI_Samples_filtered, paste(run_path, 'NGS_UMI_samples_filtered.csv', sep = ''))




