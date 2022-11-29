library(tidyverse)
library(argparser)

parser <- arg_parser("Commandline parser")
parser <- add_argument(
  parser,
  "--umi_plasmid_samples",
  help = "Merged UMI summary files"
)
parser <- add_argument(
  parser,
  "--mutation_classification",
  help = "Merged UMI summary files"
)

argv <- parse_args(parser)
umi_plasmid_samples <- argv$umi_plasmid_samples
mutation_classification <- argv$mutation_classification

### load data

umi_data <-
  read_tsv(umi_plasmid_samples)
plasmid_expected_mutations <-
  read.csv(mutation_classification) %>%
  mutate(Position = as.numeric(Position))

### define parameters
STR_start <- 2472
STR_end <- 2505

detected_mutations <- tibble(
  Sample = character(),
  Fragment = integer(),
  Run = character(),
  num_of_muts_UMI = numeric(),
  num_of_muts_Ref = numeric(),
  num_of_reads = numeric(),
  coverage = numeric(),
  Q_score = numeric(),
  num_of_observations = numeric(),
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

groups <- umi_data %>%
  group_by(sample, fragment, run) %>%
  summarize() %>%
  drop_na()
number_of_groups <- nrow(groups)

for (i in 1:number_of_groups) {
  Sample <- groups[[1]][i]
  Fragment <- groups[[2]][i]
  Run <- groups[[3]][i]

  plasmid_filtered <- plasmid_expected_mutations %>%
    filter(fragment == Fragment)

  data_filtered <- umi_data %>%
    filter(sample == Sample, fragment == Fragment, run == Run) %>%
    select(
      sample,
      Percent_A,
      Percent_B,
      run,
      variant_level_umi,
      variant_umi,
      pos,
      number_of_reads,
      coverage,
      Q_score,
      Sample_readable
    ) %>%
    full_join(plasmid_filtered, by = c("pos" = "Position")) %>%
    filter(pos < STR_start | pos > STR_end)

  if (100 %in% data_filtered$Percent_A) {
    data_filtered <- data_filtered %>%
      filter(type_annot != "type_b")
  }
  if (100 %in% data_filtered$Percent_B) {
    data_filtered <- data_filtered %>%
      filter(type_annot != "type_a_mut")
  }

  data_filtered <- data_filtered %>%
    mutate(mutation_type = ifelse(
      as.character(variant_umi) == TypeA,
      as.character(type_annot),
      ifelse(as.character(variant_umi) == TypeB,
        "type_b", "undefined"
      )
    )) %>%
    mutate(variant_level_umi = coalesce(variant_level_umi, 0.4))

  num_of_muts_UMI <- data_filtered %>%
    filter(!is.na(variant_umi)) %>%
    nrow()

  num_of_muts_Ref <- data_filtered %>%
    filter(!is.na(Ref)) %>%
    nrow()

  num_of_observations <- data_filtered %>%
    nrow()

  num_of_reads <- mean(data_filtered$number_of_reads, na.rm = TRUE)

  Q_score <- mean(data_filtered$Q_score, na.rm = TRUE)

  coverage <- ceiling(mean(data_filtered$coverage, na.rm = TRUE))

  positive <- num_of_muts_Ref

  # OR all positions - positions with SNP
  negative <- as.numeric(Fragment) - positive

  # Number of positions that are recognized of having the same SNP
  true_positive <- data_filtered %>%
    filter(!is.na(Ref)) %>%
    filter(as.character(mutation_type) == as.character(type_annot)) %>%
    nrow()

  # Number of positions where a SNP was found in the UMI data, but not or a different in the NGS data
  false_positive <- data_filtered %>%
    filter(is.na(Ref) |
      as.character(mutation_type) != as.character(type_annot)) %>%
    nrow()

  # Number of positions where a SNP was found in the NGS data, but not in the UMI data
  false_negative <- data_filtered %>%
    filter(is.na(variant_umi)) %>%
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

  num_of_detected_mutations <- data_filtered %>%
    select(variant_level_umi)
  num_of_detected_mutations <-
    colSums(num_of_detected_mutations > 0)

  detected_mutations <- detected_mutations %>% add_row(
    Sample = Sample,
    Fragment = Fragment,
    Run = Run,
    num_of_muts_UMI = num_of_muts_UMI,
    num_of_muts_Ref = num_of_muts_Ref,
    # num_of_deletions_Ref = length(which(data_filtered$variant_ngs == "D")), DISCUSS WITH STEFAN!
    num_of_observations = num_of_observations,
    num_of_reads = num_of_reads,
    coverage = coverage,
    Q_score = Q_score,
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
    f1_score_control = f1_score_control
  )

  write_tsv(data_filtered, paste(Sample, Fragment, Run, "detected_mutations.tsv", sep = "_"))
}

write_tsv(detected_mutations, "umi_plasmid_quality_parameters.tsv")
