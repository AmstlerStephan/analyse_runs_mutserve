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


argv <- parse_args(parser)
run <- argv$run
nanostat_summary <- argv$nanostat_summary
mutserve_summary <- argv$mutserve_summary
ngs_data <- argv$ngs_data
umi_cutoff <- ifelse(
  str_detect(run, "V14"),
  argv$umi_cutoff_V14,
  argv$umi_cutoff_R9
)

### define parameters
STR_start <- 2472
STR_end <- 2506 ### adapted to 2506 instead of 2505

sample_sets <- c("AK|SAPHIR")

mutserve_summary <-
  read_tsv(mutserve_summary, na = c("", "NA", "-"))

pipeline <- ifelse(str_detect(mutserve_summary$SAMPLE[1], "lpa"), "snakemake", "nextflow")

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
  mutate(fragment = as.character(fragment), 
         sample_fragment = paste(sample, fragment, sep = "_"))

fwd_muts <- mutserve_summary %>% 
  filter(!is.na(`TOP-FWD`)) %>%
  nrow()
rev_muts <- mutserve_summary %>% 
  filter(!is.na(`TOP-REV`)) %>%
  nrow()

if(fwd_muts < rev_muts){
  mutserve_summary_parsed <- mutserve_summary %>%
    mutate(minor_variant_umi = `MINOR-REV`,
           minor_variant_level_umi = `MINOR-REV-PERCENT`,
           top_variant_umi = `TOP-REV`,
           top_variant_level_umi = `TOP-REV-PERCENT`,
           coverage = `COV-TOTAL`,
           ref_umi = `REF`,
           pos = POS) %>% 
    select(minor_variant_umi,
           minor_variant_level_umi,
           top_variant_umi,
           top_variant_level_umi,
           coverage,
           ref_umi,
           pos,
           SAMPLE)
} else {
  mutserve_summary_parsed <- mutserve_summary %>%
    mutate(minor_variant_umi = `MINOR-FWD`,
           minor_variant_level_umi = `MINOR-FWD-PERCENT`,
           top_variant_umi = `TOP-FWD`,
           top_variant_level_umi = `TOP-FWD-PERCENT`,
           coverage = `COV-TOTAL`,
           ref_umi = `REF`,
           pos = POS) %>% 
    select(minor_variant_umi,
           minor_variant_level_umi,
           top_variant_umi,
           top_variant_level_umi,
           coverage,
           ref_umi,
           pos,
           SAMPLE)

}

### filter mutserve data
### filter for full conversions (called variant is not the reference AND has no minor variant level OR minor variant level is below a certain threshold)
## excluded () | minor_variant_level_umi < umi_cutoff)
mutserve_raw_full_conversions <- mutserve_summary_parsed %>%
  filter(ref_umi != top_variant_umi & is.na(minor_variant_umi)) %>%
  mutate(
    variant_level_umi = 1,
    variant_umi = top_variant_umi,
  )

### drop all NA values, where no minor variant was found (Either full conversion or no variant at that position)
### Variant level can be over 50% -> take Percentage and variant accordingly
mutserve_raw_variants <- mutserve_summary_parsed %>%
  drop_na(minor_variant_umi) %>%
  mutate(
    variant_level_umi = ifelse((ref_umi == minor_variant_umi), top_variant_level_umi, minor_variant_level_umi),
    variant_umi = as.character(ifelse((ref_umi == minor_variant_umi), top_variant_umi, minor_variant_umi))
  )

mutserve_combined <- rbind(mutserve_raw_full_conversions, mutserve_raw_variants) %>%
  mutate(barcode = str_extract(SAMPLE, "barcode\\d\\d"),
         umi_cutoff = umi_cutoff)

### Join Barcodes and mutserve data

UMI <- mutserve_combined %>%
  inner_join(barcodes, by = c("barcode"))

UMI_plasmids <- UMI %>%
  filter(grepl("A_B", sample))

if(nrow(UMI_plasmids) != 0){
  UMI_plasmids_parsed <- UMI_plasmids %>%
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
  
  UMI_plasmids_filtered <- UMI_plasmids_parsed %>%
    filter(variant_level_umi >= umi_cutoff) %>%
    filter(pos < STR_start | pos > STR_end)
  
  write_tsv(UMI_plasmids_parsed, "UMI_sequencing_mutserve_plasmids.tsv")
  write_tsv(UMI_plasmids_filtered, "UMI_sequencing_mutserve_plasmids_filtered.tsv")
  
}



### Creating the NGS data
UMI_samples <- UMI %>%
  filter(grepl(sample_sets, sample))

if(nrow(UMI_samples) != 0){
  
  NGS_sample_fragment_groups <- unique(NGS$sample_fragment)
  UMI_sample_fragment_groups <- unique(UMI_samples$sample_fragment) 
  
  available_sample_fragments <- intersect(NGS_sample_fragment_groups, UMI_sample_fragment_groups)
  available_sample_fragments_parsed <- paste(available_sample_fragments, collapse = "|")

  available_UMI_samples <- UMI_samples %>%
    filter(str_detect(sample_fragment, available_sample_fragments_parsed))
  
  available_NGS_samples <- NGS %>% 
    filter(str_detect(sample_fragment, available_sample_fragments_parsed))
  
  NGS_UMI_samples <- available_UMI_samples %>% 
    full_join(available_NGS_samples)
  
  NGS_UMI_samples_parsed <- NGS_UMI_samples %>% 
    rename(
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
      number_of_reads,
      coverage,
      Q_score,
      umi_cutoff,
    ) %>%
    mutate(
      variant_level_umi = coalesce(variant_level_umi, 0),
      variant_level_ngs = coalesce(variant_level_ngs, 0),
      variance_level_absolute_difference = variant_level_ngs - variant_level_umi,
    )
  
  NGS_UMI_samples_parsed_filtered <- 
    NGS_UMI_samples_parsed %>% 
    mutate(variant_level_umi = ifelse(variant_level_umi < umi_cutoff & !is.na(variant_ngs), 0, variant_level_umi)) %>%
    filter(variant_level_umi >= umi_cutoff | variant_level_umi == 0) %>% 
    filter(position < STR_start | position > STR_end)
  
    # filter(variant_level_umi != 1 | !is.na(variant_ngs)) %>% 
    # filter(position != 1659)
    # 
  write_tsv(UMI_samples, paste0("UMI_sequencing_samples_corresponding_position_", run, ".tsv"))
  write_tsv(NGS_UMI_samples_parsed, "NGS_UMI_samples.tsv")
  write_tsv(NGS_UMI_samples_parsed_filtered, "NGS_UMI_samples_filtered.tsv")
  
}

write_tsv(UMI, paste0("UMI_sequencing_mutserve_all_", run, ".tsv"))
