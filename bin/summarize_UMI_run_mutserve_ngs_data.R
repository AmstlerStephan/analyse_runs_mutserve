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
STR_end <- 2505

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
  select(!Sample)

NGS <- read_csv(ngs_data) %>%
  filter(pos < STR_start | pos > STR_end) %>%
  mutate(fragment = as.character(fragment))

### filter mutserve data
### filter for full conversions (called variant is not the reference AND has no minor variant level OR minor variant level is below a certain threshold)
mutserve_raw_full_conversions <- mutserve_summary %>%
  filter(REF != `TOP-REV` & (is.na(`MINOR-REV`) | `MINOR-REV-PERCENT` < umi_cutoff)) %>%
  mutate(
    variant_level_umi = 1,
    variant_umi = `TOP-REV`
  )

### drop all NA values, where no minor variant was found (Either full conversion or no variant at that position)
### Variant level can be over 50% -> take Percentage and variant accordingly
mutserve_raw_variants <- mutserve_summary %>%
  drop_na(`MINOR-REV`) %>%
  mutate(
    variant_level_umi = ifelse((`REF` == `MINOR-REV`), `TOP-REV-PERCENT`, `MINOR-REV-PERCENT`),
    variant_umi = ifelse((`REF` == `MINOR-REV`), `TOP-REV`, `MINOR-REV`)
  )

mutserve_combined <- bind_rows(mutserve_raw_full_conversions, mutserve_raw_variants) %>%
  mutate(barcode = str_extract(SAMPLE, "barcode\\d\\d"),
        umi_cutoff = umi_cutoff) %>%
  dplyr::rename(pos = POS)

### Join Barcodes and mutserve data

UMI <- mutserve_combined %>%
  inner_join(barcodes, by = c("barcode")) %>%
  dplyr::rename(
    coverage = `COV-TOTAL`,
    Q_score = mean_qual,
    ref_UMI = REF
  )

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
UMI_Samples <- UMI %>%
  filter(grepl(sample_sets, sample))

UMI_Samples_groups <- UMI_Samples %>%
  group_by(sample, fragment) %>%
  summarize()

### the groups are used to filter the NGS data before joining both dataframes
### exclude all samples that are not covered with the UMI run

NGS_Samples <- NGS %>%
  filter(sample %in% UMI_Samples_groups$sample & fragment %in% UMI_Samples_groups$fragment)

NGS_UMI_Samples <- NGS_Samples %>%
  full_join(UMI_Samples,
    by = c("fragment", "sample", "pos")
  ) %>%
  dplyr::rename(
    position = pos,
    variant_ngs = variant,
    variant_level_ngs = variant_level,
    ref_NGS = ref
  ) %>%
  select(
    sample,
    fragment,
    run,
    position,
    ref_UMI,
    variant_umi,
    variant_level_umi,
    ref_NGS,
    variant_ngs,
    variant_level_ngs,
    run,
    number_of_reads,
    coverage,
    Q_score,
    umi_cutoff
  ) %>%
  mutate(
    variant_level_umi = coalesce(variant_level_umi, 0),
    variant_level_ngs = coalesce(variant_level_ngs, 0),
    variance_level_absolute_difference = variant_level_ngs - variant_level_umi,
    # Variance_level_relative_difference = (variant_level_ngs / variant_level_umi - 1)
  )

NGS_UMI_Samples_filtered <- NGS_UMI_Samples %>%
  filter(!is.na(variant_ngs) | (variant_level_umi > umi_cutoff | variant_level_umi == 0)) %>%
  filter(position < STR_start | position > STR_end) %>%
  filter(variant_ngs != "D")


write_tsv(UMI, paste0("UMI_sequencing_mutserve_all_", run, ".tsv"))
write_tsv(UMI_plasmids, "UMI_sequencing_mutserve_plasmids.tsv")
write_tsv(UMI_plasmids_filtered, "UMI_sequencing_mutserve_plasmids_filtered.tsv")
write_tsv(NGS_UMI_Samples, "NGS_UMI_samples.tsv")
write_tsv(NGS_UMI_Samples_filtered, "NGS_UMI_samples_filtered.tsv")
