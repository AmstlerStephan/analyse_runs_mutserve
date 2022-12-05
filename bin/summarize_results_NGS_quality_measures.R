library(plyr)
library(tidyr)
library(readr)
library(jsonlite)
library(tidyverse)
library(BlandAltmanLeh)
library(argparser)

parser <- arg_parser("Commandline parser")
parser <- add_argument(
  parser,
  "--ngs_umi_samples",
  help = "Merged UMI summary files"
)

argv <- parse_args(parser)
ngs_umi_samples <- argv$ngs_umi_samples

# ngs_umi_samples <- "~/UMI_LPA_KIV2/results_adapted_nf_pipeline_20221205/umi_summary_files_per_run/run12_V14/NGS_UMI_samples.tsv"

### load data
ngs_data <-
  read_tsv(ngs_umi_samples)

groups <- ngs_data %>%
  group_by(sample, original_fragment, run) %>%
  summarize() %>%
  drop_na()

number_of_groups <- nrow(groups)

detected_mutations <- tibble(
  Sample = character(),
  Fragment = integer(),
  Run = character(),
  num_of_muts_UMI = numeric(),
  num_of_muts_NGS = numeric(),
  num_of_deletions_NGS = numeric(),
  num_of_observations = numeric(),
  num_of_reads = numeric(),
  coverage = numeric(),
  Q_score = numeric(),
  positve = numeric(),
  negative = numeric(),
  true_positive = numeric(),
  true_negative = numeric(),
  false_positive = numeric(),
  false_negative = numeric(),
  sensitivity_true_positive_rate = numeric(),
  specificity_true_negative_rate = numeric(),
  precision_positive_predictive_value = numeric(),
  f1_score = numeric(),
  f1_score_control = numeric()
)

# For loop to create data per Sample and Run
for (i in 1:number_of_groups) {
  if (number_of_groups < 1) {
    break()
  }

  Sample <- groups[[1]][i]
  Fragment <- groups[[2]][i]
  Run <- groups[[3]][i]

  data_filtered <- ngs_data %>%
    filter(sample == Sample, original_fragment == Fragment)

  num_of_observations <- nrow(data_filtered)

  num_of_reads <- mean(data_filtered$number_of_reads, na.rm = TRUE)

  Q_score <- mean(data_filtered$Q_score, na.rm = TRUE)

  coverage <- ceiling(mean(data_filtered$coverage, na.rm = TRUE))

  # View(data_filtered %>% filter(as.character(variant_umi) != as.character(variant_ngs)))

  positive <- data_filtered %>%
    filter(!is.na(variant_ngs)) %>%
    nrow()

  # OR all positions - positions with SNP
  negative <- as.numeric(Fragment) - positive

  # Number of positions that are recognized of having the same SNP
  true_positive <- data_filtered %>%
    filter(!is.na(variant_ngs)) %>%
    filter(as.character(variant_umi) == as.character(variant_ngs)) %>%
    nrow()

  # Number of positions where a SNP was found in the UMI data, but not or a different in the NGS data
  false_positive <- data_filtered %>%
    filter(is.na(variant_ngs) |
      as.character(variant_umi) != as.character(variant_ngs)) %>%
    nrow()

  # Number of positions where a SNP was found in the NGS data, but not in the UMI data
  false_negative <- data_filtered %>%
    filter(is.na(variant_umi)) %>%
    nrow()

  # All positions - Number of positions where a variant was found in the NGS data (!is.na(NGS))
  true_negative <-
    as.numeric(Fragment) - false_positive - false_negative - true_positive

  # Specificity, Precision, Recall and Sensitivity ( + F1-score)
  sensitivity_true_positive_rate <- true_positive / (true_positive + false_negative)
  specificity_true_negative_rate <- true_negative / (true_negative + false_positive)
  precision_positive_predictive_value <-
    true_positive / (true_positive + false_positive)
  f1_score <-
    2 * precision_positive_predictive_value * sensitivity_true_positive_rate / (precision_positive_predictive_value + sensitivity_true_positive_rate)
  f1_score_control <-
    true_positive / (true_positive + 0.5 * (false_positive + false_negative))

  detected_mutations <- detected_mutations %>% add_row(
    Sample = Sample,
    Fragment = Fragment,
    Run = Run,
    num_of_muts_UMI = length(which(data_filtered$variant_level_umi > 0)),
    num_of_muts_NGS = length(which(data_filtered$variant_level_ngs > 0)),
    num_of_deletions_NGS = length(which(data_filtered$variant_ngs == "D")),
    num_of_observations = num_of_observations,
    positve = positive,
    negative = negative,
    true_positive = true_positive,
    true_negative = true_negative,
    false_positive = false_positive,
    false_negative = false_negative,
    sensitivity_true_positive_rate = sensitivity_true_positive_rate,
    specificity_true_negative_rate = specificity_true_negative_rate,
    precision_positive_predictive_value = precision_positive_predictive_value,
    f1_score = f1_score,
    f1_score_control = f1_score_control,
    num_of_reads = num_of_reads,
    coverage = coverage,
    Q_score = Q_score
  )
}

write_tsv(detected_mutations, "umi_ngs_quality_parameters.tsv")
