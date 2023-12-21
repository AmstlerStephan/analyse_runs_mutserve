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
  "--sample_sheet",
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

mutserve_summary <-
  read_tsv(mutserve_summary, na = c("", "NA", "-")) %>% 
  mutate(
    fragment = 5104
  )

barcodes <-
  read_tsv(sample_sheet) %>%
  mutate(
    barcode = str_replace(Barcode, "NB", "barcode"), 
    sample = Sample, 
    fragment = 5104
  ) %>% 
  select(!c(Barcode, Sample)) 

NGS <- read_tsv(ngs_data) %>%
  filter(
    Filter == "PASS" & VariantLevel >= 0.0085
  ) %>% 
  rowwise() %>% 
  mutate(
    sample = str_split(ID, "\\.")[[1]][1],
    fragment = 5104,
    major_minor = paste(MajorBase, MinorBase, sep = "/"), 
  ) %>% 
  dplyr::rename(
    pos = Pos, 
    ref = Ref, 
    variant = Variant, 
    variant_level = VariantLevel
  )

mutserve_summary_parsed <- mutserve_summary %>%
  dplyr::rename(minor_variant_umi = `MINOR-FWD`,
                minor_variant_level_umi = `MINOR-FWD-PERCENT`,
                top_variant_umi = `TOP-FWD`,
                top_variant_level_umi = `TOP-FWD-PERCENT`,
                coverage = `COV-TOTAL`,
                ref_umi = `REF`,
                pos = POS) %>% 
  select(
    SAMPLE,
    pos,
    coverage,
    ref_umi,
    top_variant_umi,
    minor_variant_umi,
    top_variant_level_umi,
    minor_variant_level_umi,
  )

### filter mutserve data
mutserve_combined <- 
  mutserve_summary_parsed %>% 
  filter(minor_variant_level_umi > 0 | ref_umi != top_variant_umi) %>% 
  mutate(
    variant_level_umi = ifelse(ref_umi != top_variant_umi, top_variant_level_umi, minor_variant_level_umi),
    variant_umi = as.character(ifelse(ref_umi != top_variant_umi, top_variant_umi, minor_variant_umi))
  ) %>%
  mutate(barcode = str_extract(SAMPLE, "barcode\\d\\d"),
         umi_cutoff = umi_cutoff)

### Join Barcodes and mutserve data

UMI_samples <- mutserve_combined %>%
  inner_join(barcodes, by = c("barcode"))

if(nrow(UMI_samples) != 0){
  
  NGS_sample_fragment_groups <- unique(NGS$sample)
  UMI_sample_fragment_groups <- unique(UMI_samples$sample) 
  
  available_samples <- intersect(NGS_sample_fragment_groups, UMI_sample_fragment_groups)
  available_samples_parsed <- paste(available_samples, collapse = "|")
  
  available_UMI_samples <- UMI_samples %>%
    filter(str_detect(sample, available_samples_parsed))
  
  available_NGS_samples <- NGS %>% 
    filter(str_detect(sample, available_samples_parsed))
  
  NGS_UMI_samples_0 <- available_UMI_samples %>% 
    merge(available_NGS_samples, by = c("pos", "sample"), all = TRUE)
  
  
  ### CONTINUE HERE -> Too many matches
  NGS_UMI_samples <- available_UMI_samples %>% 
    full_join(available_NGS_samples, relationship = "one-to-one")
  
  NGS_UMI_samples_parsed <- NGS_UMI_samples %>% 
    dplyr::rename(
      position = pos,
      variant_ngs = variant,
      variant_level_ngs = variant_level,
      ref_ngs = ref
    ) %>%
    select(
      sample,
      fragment,
      position,
      ref_umi,
      variant_umi,
      variant_level_umi,
      ref_ngs,
      variant_ngs,
      variant_level_ngs,
      coverage,
      umi_cutoff,
    ) %>%
    mutate(
      variant_level_umi = coalesce(variant_level_umi, 0),
      variant_level_ngs = coalesce(variant_level_ngs, 0),
      variance_level_absolute_difference = variant_level_ngs - variant_level_umi,
    )
  
  # filter for:
  # umi variant below the threshold, but keep all ngs variants
  # STR positions
  # full conversions that are not considered as variants in the ngs dataset
  NGS_UMI_samples_parsed_filtered <- 
    NGS_UMI_samples_parsed %>% 
    filter(variant_level_umi >= umi_cutoff | !is.na(variant_ngs)) %>% 
    filter(position < STR_start | position > STR_end) %>% 
    filter(variance_level_absolute_difference > -1 )
  
  write_tsv(UMI_samples, paste0("UMI_sequencing_samples_corresponding_position_", run, ".tsv"))
  write_tsv(NGS_UMI_samples_parsed, "NGS_UMI_samples.tsv")
  write_tsv(NGS_UMI_samples_parsed_filtered, "NGS_UMI_samples_filtered.tsv")
  
}

write_tsv(UMI, paste0("UMI_sequencing_mutserve_all_", run, ".tsv"))

