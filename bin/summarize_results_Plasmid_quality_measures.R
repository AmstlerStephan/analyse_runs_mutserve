library(tidyverse)
library(argparser)

parser <- arg_parser("Commandline parser")
parser <- add_argument(
  parser,
  "--UMI_Plasmid_Samples",
  help = "Merged UMI summary files"
)
parser <- add_argument(
  parser,
  "--mutation_classification",
  help = "Merged UMI summary files"
)

argv <- parse_args(parser)
UMI_Plasmid_Samples <- argv$UMI_Plasmid_Samples
mutation_classification <- argv$mutation_classification

### load data

umi_data <-
  read_tsv(UMI_Plasmid_Samples)
plasmid_expected_mutations <-
  read.csv(mutation_classification) %>%
  mutate(
    Position = as.numeric(as.character(Position)),
    Corresponding_Position = as.numeric(as.character(Corresponding_Position))
  )

get_groups <- function(data) {
  groups <- data %>%
    group_by(sample, fragment, run) %>%
    summarize() %>%
    drop_na()
  return(groups)
}

detected_mutations <- tibble(
  Sample = character(),
  Fragment = integer(),
  Run = character(),
  num_of_muts_UMI = numeric(),
  num_of_muts_Ref = numeric(),
  num_of_reads = numeric(),
  num_of_consensus_sequences = numeric(),
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

groups <- get_groups(umi_data)
number_of_groups <- nrow(groups)
path <- paste(result_folder,
  dir_figures,
  data_type,
  dir_Plasmid,
  sep = "/"
)

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
      Variant_level_UMI,
      Variant_UMI,
      pos,
      number_of_reads,
      `COV.TOTAL`,
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
      as.character(Variant_UMI) == TypeA,
      as.character(type_annot),
      ifelse(as.character(Variant_UMI) == TypeB,
        "type_b", "undefined"
      )
    )) %>%
    mutate(Variant_level_UMI = coalesce(Variant_level_UMI, 0.4))

  num_of_muts_UMI <- data_filtered %>%
    filter(!is.na(Variant_UMI)) %>%
    nrow()

  num_of_muts_Ref <- data_filtered %>%
    filter(!is.na(Ref)) %>%
    nrow()

  num_of_observations <- data_filtered %>%
    nrow()

  num_of_reads <- mean(data_filtered$number_of_reads, na.rm = TRUE)

  Q_score <- mean(data_filtered$Q_score, na.rm = TRUE)

  num_of_consensus_sequences <- ceiling(mean(data_filtered$`COV.TOTAL`, na.rm = TRUE))

  positive <- num_of_muts_Ref

  # OR all Positions - positions with SNP
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
    filter(is.na(Variant_UMI)) %>%
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
    select(Variant_level_UMI)
  num_of_detected_mutations <-
    colSums(num_of_detected_mutations > 0)

  detected_mutations <- detected_mutations %>% add_row(
    Sample = Sample,
    Fragment = Fragment,
    Run = Run,
    num_of_muts_UMI = num_of_muts_UMI,
    num_of_muts_Ref = num_of_muts_Ref,
    # num_of_deletions_Ref = length(which(data_filtered$Variant_NGS == "D")), DISCUSS WITH STEFAN!
    num_of_observations = num_of_observations,
    num_of_reads = num_of_reads,
    num_of_consensus_sequences = num_of_consensus_sequences,
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

  write_tsv(data_filtered, paste(Sample, Fragment, Run, "detected_mutations.csv", sep = "_"))
}

write_tsv(detected_mutations, "umi_plasmid_quality_parameters.tsv")
