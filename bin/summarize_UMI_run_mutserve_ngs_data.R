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
  "--corresponding_positions",
  help = "corresponding mutations for the fragments"
)


argv <- parse_args(parser)
run <- argv$run
nanostat_summary <- argv$nanostat_summary
mutserve_summary <- argv$mutserve_summary
ngs_data <- argv$ngs_data
corresponding_positions <- argv$corresponding_positions
umi_cutoff <- ifelse(
  str_detect(run, "V14"),
  argv$umi_cutoff_V14,
  argv$umi_cutoff_R9
)

join_NGS_reference_data_UMI_missing_samples <- function(UMI_data) {
    NGS %>% 
      filter(grepl(missing_samples_parsed_corresponding_samples_available_parsed, sample_fragment)) %>%
      filter(pos <= overlap_5104_ending_position_front | 
               pos >= overlap_5104_starting_position_end) %>% 
      full_join(UMI_data,
                by = c("sample", "fragment", "pos"))
}
join_NGS_reference_data_UMI_available_samples <- function(UMI_data) {
    NGS %>% 
      filter(grepl(available_sample_fragments_parsed, sample_fragment)) %>%
      full_join(UMI_data,
                by = c("sample", "fragment", "pos"))
}
parse_NGS_UMI_samples <- function(UMI_data){
  UMI_data %>%
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
      corresponding_position,
      reference_available
    ) %>%
    mutate(
      variant_level_umi = coalesce(variant_level_umi, 0),
      variant_level_ngs = coalesce(variant_level_ngs, 0),
      variance_level_absolute_difference = variant_level_ngs - variant_level_umi,
    )
}


run <- "run12_V14"
mutserve_summary <- "run12_V14/ont_pl/mutserve/run12_V14_summary_mutserve.txt"
nanostat_summary <- "~/post_pipeline_analysis/QC/Nanostat_parsed_merged/run12_V14/run12_V14_1000_9.tsv"
ngs_data <- "data_ngs/data_ngs/20221122_NGS_reference_data_SAPHIR.csv"
corresponding_positions <- "data_ngs/data_ngs/20221129_corresponding_positions.csv"
umi_cutoff <- 0.005

### define parameters
STR_start <- 2472
STR_end <- 2506 ### adapted to 2506 instead of 2505
overlap_2645_ending_position_front <- 982 # inclusive 982 position are overlapping
overlap_2645_starting_position_end <- 1431 # from inclusive 1431 positions are overlapping 
overlap_5104_ending_position_front <- 2645 - overlap_2645_starting_position_end + 1 # inclusive 1215 positions are overlapping
overlap_5104_starting_position_end <- 5104 - overlap_2645_ending_position_front + 1 # form inclusive 4123 positions are overlapping


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
  mutate(fragment = as.character(fragment), 
         sample_fragment = paste(sample, fragment, sep = "_"))

mutserve_summary_parsed <- mutserve_summary %>%
  mutate(minor_variant_umi = ifelse(is.na(`MINOR-REV`), `MINOR-FWD`, `MINOR-REV`),
         minor_variant_level_umi = ifelse(is.na(`MINOR-REV`), `MINOR-FWD-PERCENT`, `MINOR-REV-PERCENT`),
         top_variant_umi = ifelse(is.na(`TOP-REV`), `TOP-FWD`, `TOP-REV`),
         top_variant_level_umi = ifelse(is.na(`TOP-REV`), `TOP-FWD-PERCENT`, `TOP-REV-PERCENT`),
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
  
corresponding_positions <- 
  read_csv(corresponding_positions)

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
UMI_samples_temp <- UMI %>%
  filter(grepl(sample_sets, sample))

if(nrow(UMI_samples_temp) != 0){
  
  UMI_samples_2645 <- UMI_samples_temp %>% 
    filter(fragment == "2645") %>% 
    left_join(corresponding_positions, by = c("pos" = "pos_2645")) %>% 
    mutate(corresponding_position = pos_5104) %>% 
    select(!pos_5104)
  
  UMI_samples_5104 <- UMI_samples_temp %>% 
    filter(fragment == "5104") %>% 
    left_join(corresponding_positions, by = c("pos" = "pos_5104")) %>% 
    mutate(corresponding_position = pos_2645) %>% 
    select(!pos_2645)
  
  UMI_samples <- rbind(UMI_samples_2645, UMI_samples_5104)
  
  NGS_sample_fragment_groups <- unique(NGS$sample_fragment)
  UMI_sample_fragment_groups <- unique(UMI_samples$sample_fragment) 
  
  available_sample_fragments <- intersect(NGS_sample_fragment_groups, UMI_sample_fragment_groups)
  missing_sample_fragments <- setdiff(UMI_sample_fragment_groups, available_sample_fragments)
  
  ### get missing samples and check if the corresponding fragment exists in the NGS data
  if(!is_empty(missing_sample_fragments)){
    
    missing_samples <- str_sub(missing_sample_fragments, end = -6)
    missing_samples_parsed <- paste(missing_samples, collapse = "|")
    missing_samples_parsed_corresponding_samples_available <- 
      NGS_sample_fragment_groups[
        grepl(missing_samples_parsed, 
              NGS_sample_fragment_groups)
      ]
    ### filtet missing samples for available corresponding samples of the NGS data
    missing_samples_parsed_corresponding_samples_available_parsed <- 
      paste(missing_samples_parsed_corresponding_samples_available, collapse =  "|")
    
    missing_samples_parsed_corresponding_samples_available_samplename <- 
      str_sub(missing_samples_parsed_corresponding_samples_available, end = -6) 
    
    missing_samples_parsed_corresponding_samples_available_samplename_parsed <- 
      paste(missing_samples_parsed_corresponding_samples_available_samplename, collapse = "|")
    
    missing_samples_corresponding_samples_available_filtered <-
      missing_sample_fragments[grepl(
        missing_samples_parsed_corresponding_samples_available_samplename_parsed,
        missing_sample_fragments)
      ]
    missing_samples_corresponding_samples_available_filtered_parsed <- 
      paste(missing_samples_corresponding_samples_available_filtered, collapse = "|")
    
    UMI_samples_missing <- UMI_samples %>% 
      filter( grepl(missing_samples_corresponding_samples_available_filtered_parsed, sample_fragment) ) %>% 
      mutate(original_position = pos,
             pos = corresponding_position,
             original_fragment = fragment,
             fragment = ifelse(fragment == "5104", "2645", "5104"),
             sample_fragment = paste(sample, fragment, sep = "_"),
             reference_available = FALSE,
      ) %>% 
      drop_na(corresponding_position)
    
    UMI_samples_missing_filtered <- UMI_samples_missing %>% 
      filter(variant_level_umi >= umi_cutoff)
    
    
    NGS_UMI_missing_samples <- join_NGS_reference_data_UMI_missing_samples(UMI_samples_missing)
    NGS_UMI_missing_samples_filtered <- join_NGS_reference_data_UMI_missing_samples(UMI_samples_missing_filtered)   
  } else {
    NGS_UMI_missing_samples <- NA
    NGS_UMI_missing_samples_filtered <- NA
  }
  if (!is_empty(available_sample_fragments)){
    available_sample_fragments_parsed <- paste(available_sample_fragments, collapse = "|")
    
    UMI_samples_available <- UMI_samples %>% 
      filter(grepl(available_sample_fragments_parsed, sample_fragment)) %>% 
      mutate(original_fragment = fragment,
             original_position = pos,
             sample_fragment = paste(sample, fragment, sep = "_"),
             reference_available = TRUE)
    
    UMI_samples_available_filtered <- UMI_samples_available %>% 
      filter(variant_level_umi >= umi_cutoff)
    
    NGS_UMI_available_samples <- join_NGS_reference_data_UMI_available_samples(UMI_samples_available)
    NGS_UMI_available_samples_filtered <- join_NGS_reference_data_UMI_available_samples(UMI_samples_available_filtered)

  } else {
    NGS_UMI_available_samples <- NA
    NGS_UMI_available_samples_filtered <- NA
  }
 
  NGS_UMI_samples <- rbind(NGS_UMI_available_samples, NGS_UMI_missing_samples)
  NGS_UMI_samples_filtered <- rbind(NGS_UMI_available_samples_filtered, NGS_UMI_missing_samples_filtered)

  NGS_UMI_samples_parsed <- parse_NGS_UMI_samples(NGS_UMI_samples)
  NGS_UMI_samples_parsed_filtered <- parse_NGS_UMI_samples(NGS_UMI_samples_filtered) %>%
    filter(position < STR_start | position > STR_end) %>%
    filter(variant_ngs != "D" | is.na(variant_ngs)) %>% 
    filter(variant_level_umi != 1 | !is.na(variant_ngs)) %>% 
    filter(variant_umi != "D" | is.na(variant_umi)) %>% 
    filter(position != 1659)
  
  write_tsv(UMI_samples, paste0("UMI_sequencing_samples_corresponding_position_", run, ".tsv"))
  write_tsv(NGS_UMI_samples_parsed, "NGS_UMI_samples.tsv")
  write_tsv(NGS_UMI_samples_parsed_filtered, "NGS_UMI_samples_filtered.tsv")

}

write_tsv(UMI, paste0("UMI_sequencing_mutserve_all_", run, ".tsv"))
