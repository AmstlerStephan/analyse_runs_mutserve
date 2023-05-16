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


# mutserve_summary <- "M:/Masterarbeit_GENEPI/20230324_mutserve/indel/HAC/run10_V14/barcode12/barcode12_mutserve.txt"
# umi_cutoff <- 0.0085
# ngs_data <- "M:/Masterarbeit_GENEPI/data_ngs/20221122_NGS_reference_data_SAPHIR.csv"
# nanostat_summary <- "M:/Masterarbeit_GENEPI/20230324_QC_ALL/HAC/Nanostat_parsed_merged/all_runs_1000_9.tsv"

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
  mutate(fragment = as.character(fragment), 
         sample_fragment = paste(sample, fragment, sep = "_"))


mutserve_combined <- mutserve_summary %>%
  mutate(
    variant_umi = Variant, 
    variant_level_umi = VariantLevel,
    coverage = Coverage,
    ref_umi = Ref,
    pos = Pos, 
    barcode = str_extract(ID, "barcode\\d\\d"), 
    umi_cutoff = umi_cutoff) %>% 
  select(
    barcode, 
    pos, 
    ref_umi,
    variant_umi, 
    variant_level_umi, 
    coverage,
    umi_cutoff,
    )

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
     merge(available_NGS_samples, by = c("pos", "sample", "fragment"), all = TRUE)
  
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
  
  # filter for:
  # umi variant below the threshold, but keep all ngs variants
  # STR positions
  # full conversions that are not considered as variants in the ngs dataset
  NGS_UMI_samples_parsed_filtered <- 
    NGS_UMI_samples_parsed %>% 
    filter(variant_level_umi >= umi_cutoff | variant_level_ngs > 0) %>% 
    filter(position < STR_start | position > STR_end) %>% 
    filter(!(is.na(variant_ngs) & variant_level_umi == 1 ))

  write_tsv(UMI_samples, paste0("UMI_sequencing_samples_corresponding_position_", run, ".tsv"))
  write_tsv(NGS_UMI_samples_parsed, "NGS_UMI_samples.tsv")
  write_tsv(NGS_UMI_samples_parsed_filtered, "NGS_UMI_samples_filtered.tsv")
  
}

write_tsv(UMI, paste0("UMI_sequencing_mutserve_all_", run, ".tsv"))
