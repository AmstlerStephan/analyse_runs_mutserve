library(plyr)
library(tidyr)
library(readr)
library(jsonlite)
library(tidyverse)
library(argparser)

parser <- arg_parser("Commandline parser")
parser <- add_argument(
  parser,
  "--run",
  help = "Run"
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
parser <- add_argument(
  parser,
  "--umi_cutoff_R9",
  help = "cutoff value for filtering reads"
)
parser <- add_argument(
  parser,
  "--umi_cutoff_V14",
  help = "cutoff value for filtering reads"
)
parser <- add_argument(
  parser,
  "--corresponding_mutations",
  help = "corresponding mutations for the fragments"
)


argv <- parse_args(parser)
run <- argv$run
nanostat_summary <- argv$nanostat_summary
mutserve_summary <- argv$mutserve_summary
ngs_data <- argv$ngs_data
corresponding_mutations <- argv$corresponding_mutations
umi_cutoff <- ifelse(
  str_detect(run, "V14"),
  argv$umi_cutoff_V14,
  argv$umi_cutoff_R9
)

### define parameters
STR_start <- 2472
STR_end <- 2505
overlap_2645_end <- 982
overlap_2645_start <- 1431


sample_sets <- c("AK|SAPHIR")
SAPHIR_samples <- c("4612|4901|4451|4624|4864|5248|5538|5400")
AK_samples <- c("AK03|AK07|AK14|AK17|AK33")

mutserve_summary <-
  read_tsv(mutserve_summary, na = c("", "NA", "-"))

barcodes <-
  read_tsv(nanostat_summary) %>%
  select(run:Sample, number_of_reads, mean_qual) %>%
  mutate(
    sample = str_sub(Sample, end = -6),
    fragment = str_sub(Sample, start = -4)
  ) %>%
  dplyr::rename(
    sample_fragment = Sample,
    Q_score = mean_qual
  )

NGS <- read_csv(ngs_data) %>%
  filter(pos < STR_start | pos > STR_end) %>%
  mutate(fragment = as.character(fragment))

mutserve_summary_parsed <- mutserve_summary %>%
  mutate(minor_variant_umi = ifelse(is.na(`MINOR-REV`), `MINOR-FWD`, `MINOR-REV`),
         minor_variant_level_umi = ifelse(is.na(`MINOR-REV`), `MINOR-FWD-PERCENT`, `MINOR-REV-PERCENT`),
         top_variant_umi = ifelse(is.na(`TOP-REV`), `TOP-FWD`, `TOP-REV`),
         top_variant_level_umi = ifelse(is.na(`TOP-REV`), `TOP-FWD-PERCENT`, `TOP-REV-PERCENT`),
         coverage = `COV-TOTAL`,
         ref_umi = `REF`,
         pos = POS)
  
corresponding_position <- 
  read_csv(corresponding_position)

### filter mutserve data
### filter for full conversions (called variant is not the reference AND has no minor variant level OR minor variant level is below a certain threshold)
mutserve_raw_full_conversions <- mutserve_summary_parsed %>%
  filter(ref_umi != top_variant_umi & (is.na(minor_variant_umi) | minor_variant_level_umi < umi_cutoff)) %>%
  mutate(
    variant_level_umi = 1,
    variant_umi = top_variant_umi
  )

### drop all NA values, where no minor variant was found (Either full conversion or no variant at that position)
### Variant level can be over 50% -> take Percentage and variant accordingly
mutserve_raw_variants <- mutserve_summary_parsed %>%
  drop_na(minor_variant_umi) %>%
  mutate(
    variant_level_umi = ifelse((ref_umi == minor_variant_umi), top_variant_level_umi, minor_variant_level_umi),
    variant_umi = as.character(ifelse((ref_umi == minor_variant_umi), top_variant_umi, minor_variant_umi))
  )

mutserve_combined <- bind_rows(mutserve_raw_full_conversions, mutserve_raw_variants) %>%
  mutate(barcode = str_extract(SAMPLE, "barcode\\d\\d"),
        umi_cutoff = umi_cutoff)

### Join Barcodes and mutserve data

UMI <- mutserve_combined %>%
  inner_join(barcodes, by = c("barcode"))

UMI_plasmids <- UMI %>%
  filter(grepl("A_B", sample)) %>%
  separate(sample,
    c(NA, NA, "Percent_A", "Percent_B"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    Percent_A = as.numeric(Percent_A) / 10,
    Percent_B = as.numeric(Percent_B) / 10,
    Sample_readable = paste(Percent_A, Percent_B, sep = ":")
  )

UMI_plasmids_filtered <- UMI_plasmids %>%
  filter(variant_level_umi > umi_cutoff) %>%
  filter(pos < STR_start | pos > STR_end)


### Creating the NGS data
UMI_Samples_temp <- UMI %>%
  filter(grepl(sample_sets, sample))

UMI_Samples_2645 <- UMI_Samples_temp %>% 
  filter(fragment == "2645") %>% 
  left_join(corresponding_position, by = c("pos" = "pos_2645")) %>% 
  mutate(corresponding_position = pos_5104) %>% 
  select(!pos_5104)

UMI_Samples_5104 <- UMI_Samples_temp %>% 
  filter(fragment == "5104") %>% 
  left_join(corresponding_position, by = c("pos" = "pos_5104")) %>% 
  mutate(corresponding_position = pos_2645) %>% 
  select(!pos_2645)

UMI_Samples <- bind_rows(UMI_Samples_2645, UMI_Samples_5104)

NGS_Samples_group <- NGS %>%
  transmute(sample_fragment = paste(sample, fragment, sep = "_")) %>% 
  unique()

UMI_Samples_missing <- UMI_Samples %>%
  anti_join(NGS_Samples_temp, by = "sample_fragment") %>% 
  mutate(original_position = pos,
         pos = corresponding_position,
         original_fragment = fragment,
         fragment = ifelse(fragment == "5104", "2645", "5104"),
         sample_fragment = paste(sample, fragment, sep = "_"))

UMI_Samples_available <- UMI_Samples %>%
  inner_join(NGS_Samples_temp, by = "sample_fragment") %>% 
  mutate(original_fragment = fragment,
         original_position = pos,
         sample_fragment = paste(sample, fragment, sep = "_"))

UMI_Samples_parsed <- bind_rows(UMI_Samples_available, UMI_Samples_missing)

### the groups are used to filter the NGS data before joining both dataframes
### exclude all samples that are not covered with the UMI run

NGS_Samples <- NGS %>%
  unite("sample_fragment", sample, fragment, sep = "_") %>% 
  filter(sample_fragment %in% UMI_Samples_parsed$sample_fragment) 

NGS_UMI_Samples <- NGS_Samples %>%
  full_join(UMI_Samples_parsed,
    by = c("sample_fragment", "pos")
  ) %>%
  dplyr::rename(
    position = pos,
    variant_ngs = variant,
    variant_level_ngs = variant_level,
    ref_ngs = ref
  ) %>%
  select(
    sample,
    fragment,
    run,
    position,
    ref_umi,
    variant_umi,
    variant_level_umi,
    ref_ngs,
    variant_ngs,
    variant_level_ngs,
    run,
    number_of_reads,
    coverage,
    Q_score,
    umi_cutoff,
    original_position,
    original_fragment,
    corresponding_position
  ) %>%
  mutate(
    variant_level_umi = coalesce(variant_level_umi, 0),
    variant_level_ngs = coalesce(variant_level_ngs, 0),
    variance_level_absolute_difference = variant_level_ngs - variant_level_umi,
    # Variance_level_relative_difference = (variant_level_ngs / variant_level_umi - 1)
  ) %>% 
  filter(fragment == original_fragment | original_position > overlap_2645_start | original_position < overlap_2645_end)

NGS_UMI_Samples_filtered <- NGS_UMI_Samples %>%
  filter(!is.na(variant_ngs) | (variant_level_umi > umi_cutoff | variant_level_umi == 0)) %>%
  filter(position < STR_start | position > STR_end) %>%
  filter(variant_ngs != "D")


write_tsv(UMI, paste0("UMI_sequencing_mutserve_all_", run, ".tsv"))
write_tsv(UMI_plasmids, "UMI_sequencing_mutserve_plasmids.tsv")
write_tsv(UMI_plasmids_filtered, "UMI_sequencing_mutserve_plasmids_filtered.tsv")
write_tsv(NGS_UMI_Samples, "NGS_UMI_samples.tsv")
write_tsv(NGS_UMI_Samples_filtered, "NGS_UMI_samples_filtered.tsv")
