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

### load data
ngs_data <-
  read_tsv(ngs_umi_samples)

groups <- ngs_data %>%
  group_by(sample, fragment, run) %>%
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
  num_of_consensus_sequences = numeric(),
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
  Sample <- groups[[1]][i]
  Fragment <- groups[[2]][i]
  Run <- groups[[3]][i]

  data_filtered <- ngs_data %>%
    filter(sample == Sample, fragment == Fragment)

  num_of_observations <- nrow(data_filtered)

  num_of_reads <- mean(data_filtered$number_of_reads, na.rm = TRUE)

  Q_score <- mean(data_filtered$Q_score, na.rm = TRUE)

  num_of_consensus_sequences <- ceiling(mean(data_filtered$num_of_consensus_sequences, na.rm = TRUE))

  # View(data_filtered %>% filter(as.character(variant_UMI) != as.character(variant_NGS)))

  positive <- data_filtered %>%
    filter(!is.na(variant_NGS)) %>%
    nrow()

  # OR all Positions - positions with SNP
  negative <- as.numeric(Fragment) - positive

  # Number of positions that are recognized of having the same SNP
  true_positive <- data_filtered %>%
    filter(!is.na(variant_NGS)) %>%
    filter(as.character(variant_UMI) == as.character(variant_NGS)) %>%
    nrow()

  # Number of positions where a SNP was found in the UMI data, but not or a different in the NGS data
  false_positive <- data_filtered %>%
    filter(is.na(variant_NGS) |
      as.character(variant_UMI) != as.character(variant_NGS)) %>%
    nrow()

  # Number of positions where a SNP was found in the NGS data, but not in the UMI data
  false_negative <- data_filtered %>%
    filter(is.na(variant_UMI)) %>%
    nrow()

  # All positions - Number of positions where a variant was found in the NGS data (!is.na(NGS))
  true_negative <-
    as.numeric(Fragment) - false_positive - false_negative - true_positive

  # Specificity, Precision, Recall and Sensitivity ( + F1-score)
  sensitivity_true_positive_rate <- true_positive / positive
  specificity_true_negative_rate <- true_negative / negative
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
    num_of_muts_UMI = length(which(data_filtered$variant_level_UMI > 0)),
    num_of_muts_NGS = length(which(data_filtered$variant_level_NGS > 0)),
    num_of_deletions_NGS = length(which(data_filtered$variant_NGS == "D")),
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
    num_of_consensus_sequences = num_of_consensus_sequences,
    Q_score = Q_score
  )
}

write_tsv(detected_mutations, "umi_ngs_quality_parameters.tsv")
