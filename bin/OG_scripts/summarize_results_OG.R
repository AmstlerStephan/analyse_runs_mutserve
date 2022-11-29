library(plyr)
library(tidyr)
library(readr)
library(jsonlite)
library(tidyverse)
library(BlandAltmanLeh)

### set result folder (output folder)
args <- commandArgs(TRUE)
if (!is.null(args[1])) {
  result_folder <- args[1]
  print(result_folder)
}

result_folder <- "results_all/"

### Create dirs
dir_figures <- "figures"
dir_NGS <- "NGS"
dir_Plasmid <- "Plasmid"
dir_filtered <- "filtered"
dir_raw <- "raw"
dir_reads_vs_cons_consensus_sequences <- "reads_vs_cons_consensus_sequences"
dir_Bland_Altman <- "bland_Altman"
dir_compare_NGS_vs_UMI_per_position <- "compare_variant_level_NGS_vs_UMI_per_position"
dir_compare_NGS_vs_UMI <- "compare_variant_level_NGS_vs_UMI"
dir_density_NGS_vs_UMI <- "density_mutation_levels_NGS_vs_UMI"
dir_density <- "density_mutation_levels"
dir_variance_per_pos <- "variance_per_pos"
dir_detected_mutations <- "detected_mutations"

if(!dir.exists(paste(result_folder, dir_figures, sep = "/"))) {
  
  #reads_vs_cons
  dir.create(paste(result_folder, dir_figures, dir_reads_vs_cons_consensus_sequences, sep = "/"), recursive = TRUE)
  
  # NGS
  dir.create(paste(result_folder, dir_figures, dir_raw, dir_NGS, dir_Bland_Altman, sep = "/"), recursive = TRUE)
  dir.create(paste(result_folder, dir_figures, dir_raw, dir_NGS, dir_compare_NGS_vs_UMI_per_position, sep = "/"), recursive = TRUE)
  dir.create(paste(result_folder, dir_figures, dir_raw, dir_NGS, dir_compare_NGS_vs_UMI, sep = "/"), recursive = TRUE)
  dir.create(paste(result_folder, dir_figures, dir_raw, dir_NGS, dir_density_NGS_vs_UMI, sep = "/"), recursive = TRUE)
  dir.create(paste(result_folder, dir_figures, dir_filtered, dir_NGS, dir_Bland_Altman, sep = "/"), recursive = TRUE)
  dir.create(paste(result_folder, dir_figures, dir_filtered, dir_NGS, dir_compare_NGS_vs_UMI_per_position, sep = "/"), recursive = TRUE)
  dir.create(paste(result_folder, dir_figures, dir_filtered, dir_NGS, dir_compare_NGS_vs_UMI, sep = "/"), recursive = TRUE)
  dir.create(paste(result_folder, dir_figures, dir_filtered, dir_NGS, dir_density_NGS_vs_UMI, sep = "/"), recursive = TRUE)
  
  #Plasmid
  dir.create(paste(result_folder, dir_figures, dir_raw, dir_Plasmid, dir_variance_per_pos, sep = "/"), recursive = TRUE)
  dir.create(paste(result_folder, dir_figures, dir_raw, dir_Plasmid, dir_density, sep = "/"), recursive = TRUE)
  dir.create(paste(result_folder, dir_figures, dir_raw, dir_Plasmid, dir_detected_mutations, sep = "/"), recursive = TRUE)
  dir.create(paste(result_folder, dir_figures, dir_filtered, dir_Plasmid, dir_variance_per_pos, sep = "/"), recursive = TRUE)
  dir.create(paste(result_folder, dir_figures, dir_filtered, dir_Plasmid, dir_density, sep = "/"), recursive = TRUE)
  dir.create(paste(result_folder, dir_figures, dir_filtered, dir_Plasmid, dir_detected_mutations, sep = "/"), recursive = TRUE)
}

### Set parameters
STR_start <- 2472
STR_end <- 2505

### load data
NGS_data <-
  read.csv(paste(result_folder, 'temp_all_NGS.csv', sep = "/"))
NGS_data_filtered <-
  read.csv(paste(result_folder, 'temp_all_NGS_filtered.csv', sep = "/"))
UMI_data <-
  read.csv(paste(result_folder, 'temp_all_UMI.csv', sep = "/"))
UMI_data_filtered <-
  read.csv(paste(result_folder, 'temp_all_UMI_filtered.csv', sep = "/"))
run_info <- read.csv('info/Status_overview.csv', header = TRUE) %>%
  select('Sequenced_with', 'Run_directory') %>%
  dplyr::rename(device = Sequenced_with)
plasmid_expected_mutations <-
  read.csv("data_ngs/plasmid_expected_muts.csv") %>%
  mutate(
    Position = as.numeric(as.character(Position)),
    Corresponding_Position = as.numeric(as.character(Corresponding_Position))
  )

### functions



reads_vs_cons_consensus_sequences_plot <- function(split_by) {
  reads_vs_final_bam_files %>%
    ggplot(aes(x = number_of_reads, y = mean_cov, color = split_by)) +
    geom_point() +
    geom_smooth(method = lm, se = FALSE)
}

get_groups <- function(data) {
  groups <- data %>%
    group_by(sample, fragment, run) %>%
    summarize() %>%
    drop_na()
  return(groups)
}

save_plot <- function(plot, path, Fragment, Sample, Run) {
  ggsave(filename =
           paste(
             paste(Fragment,
                   Sample,
                   Run,
                   substitute(plot),
                   sep = '_'),
             "jpg",
             sep = "."
           ),
         plot,
         path = path)
}

create_bland_altman <- function(data, path, Fragment, Sample, Run) {
  jpeg(
    file =  paste(path,
                  paste(
                    paste(Fragment,
                          Sample,
                          Run,
                          "bland_altman",
                          sep = '_'),
                    "jpg",
                    sep = "."
                  ),
                  sep = "/"),
    width = 10,
    height = 10,
    units = "in",
    res = 300
  )
  
  
  bland_stats <-
    bland.altman.stats(data$Variant_level_NGS,
                       data$Variant_level_UMI)
  
  print(
    bland.altman.plot(
      data$Variant_level_NGS,
      data$Variant_level_UMI,
      main = paste(Sample, Fragment, "Variant levels", sep = "_"),
      xlab = "Means",
      ylab = "Differences",
      graph.sys = "ggplot2",
    ) +
      geom_hline(aes(yintercept = 0)) +
      coord_cartesian(xlim = c(0, 1), ylim = c(-0.2, 0.2)) +
      annotate(
        "text",
        x = 0.7,
        y = 0.15,
        label = paste(
          "cv =",
          bland_stats$mean.diffs / mean(bland_stats$means) * 100,
          sep = " "
        )
      )
  )
  dev.off()
}

analyze_NGS_data <- function(data, data_type) {
  groups <- get_groups(data)
  number_of_groups <- nrow(groups)
  path <- paste(result_folder,
                dir_figures,
                data_type,
                dir_NGS,
                sep = "/")
  
  print(groups)
  print(number_of_groups)
  
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
    
    print(Sample)
    print(Fragment)
    print(Run)
    
    data_filtered <- data %>%
      filter(sample == Sample , fragment == Fragment , run == Run)
    
    num_of_observations <- nrow(data_filtered)
    
    num_of_reads <- mean(data_filtered$number_of_reads, na.rm = TRUE)
    
    Q_score <- mean(data_filtered$Q_score, na.rm = TRUE)
    
    num_of_consensus_sequences <- ceiling(mean(data_filtered$num_of_consensus_sequences, na.rm = TRUE))
    
    #View(data_filtered %>% filter(as.character(Variant_UMI) != as.character(Variant_NGS)))
    
    positive <- data_filtered %>%
      filter(!is.na(Variant_NGS)) %>%
      nrow()
    
    # OR all Positions - positions with SNP
    negative <- as.numeric(Fragment) - positive
    
    # Number of positions that are recognized of having the same SNP
    true_positive <- data_filtered %>%
      filter(!is.na(Variant_NGS)) %>%
      filter(as.character(Variant_UMI) == as.character(Variant_NGS)) %>%
      nrow()
    
    # Number of positions where a SNP was found in the UMI data, but not or a different in the NGS data
    false_positive <- data_filtered %>%
      filter(is.na(Variant_NGS) |
               as.character(Variant_UMI) != as.character(Variant_NGS)) %>%
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
    
    detected_mutations <- detected_mutations %>% add_row(
      Sample = Sample,
      Fragment = Fragment,
      Run = Run,
      num_of_muts_UMI = length(which(data_filtered$Variant_level_UMI > 0)),
      num_of_muts_NGS = length(which(data_filtered$Variant_level_NGS > 0)),
      num_of_deletions_NGS = length(which(data_filtered$Variant_NGS == "D")),
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
    
    # NGS variant level vs UMI variant level
    
    r_squared <-
      data_filtered %>% lm(Variant_level_UMI ~ Variant_level_NGS, data = .)
    r_squared <- summary(r_squared)$r.squared
    
    comparison_variant_levels_ngs_umi <- data_filtered %>%
      ggplot(aes(x = Variant_level_UMI, y = Variant_level_NGS)) +
      geom_abline() +
      geom_smooth(method = 'lm', show.legend = TRUE) +
      geom_point(position = 'jitter') +
      labs(
        x = 'relative variance level UMI',
        y = 'relative variance level NGS',
        title = paste(
          Sample,
          Fragment,
          Run,
          'Variant levels of both Sequencing technology',
          sep = "_"
        )
      ) +
      annotate(
        'text',
        x = 0.8,
        y = 0.05,
        label = paste0('R Squared = ', round(r_squared, digits = 3))
      ) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      scale_x_continuous(breaks = seq(0, 1, by = 0.1))
    
    # Compare variance levels of different Sequencing methods
    
    comparison_variant_levels_ngs_umi_per_position <-
      data_filtered %>%
      ggplot(aes(x = Position)) +
      geom_point(aes(y = Variant_level_NGS, color = 'NGS')) +
      geom_point(aes(y = Variant_level_UMI, color = 'UMI')) +
      labs(
        x = 'position relative to reference sequence',
        y = 'relative variant Level',
        title = paste(
          Sample,
          Fragment,
          'Variant levels per position and Sequencing technology',
          sep = "_"
        )
      ) +
    scale_x_continuous(breaks = seq(0, 5200, by = 200))
    
    density_plot_variant_levels <-
      data_filtered %>% 
      ggplot() +
      geom_density(
        aes(Variant_level_UMI),
        fill = "green",
        color = "grey",
        alpha = 0.4
      ) +
      geom_density(
        aes(Variant_level_NGS),
        fill = "blue",
        color = "grey",
        alpha = 0.4
      ) +
      labs(
        x = "relative variant level",
        y = "number of variants per level",
        title = paste(Sample,
                      Fragment,
                      'number of variants per level of UMI vs NGS',
                      sep = " ")
      ) +
      coord_cartesian(xlim = c(0, 1)) +
      scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
      annotate(
        geom = "text",
        x = 0.2,
        y = -1,
        label =
          "UMI = green NGS = blue"
      )
    
    save_plot(
      density_plot_variant_levels,
      paste(path, dir_density_NGS_vs_UMI,
            sep = "/"),
      Fragment,
      Sample,
      Run
    )
    
    
    # Bland Altman Plot NGS vs UMI variant levels
    create_bland_altman(
      data_filtered,
      paste(path,
        dir_Bland_Altman,
        sep = "/"
      ),
      Fragment,
      Sample,
      Run
    )
    
    # save plots
    save_plot(
      comparison_variant_levels_ngs_umi,
      paste(
        path,
        dir_compare_NGS_vs_UMI,
        sep = "/"
      ),
      Fragment,
      Sample,
      Run
    )
  
    
    save_plot(
      comparison_variant_levels_ngs_umi_per_position,
      paste(
        path, dir_compare_NGS_vs_UMI_per_position,
        sep = "/"
      ),
      Fragment,
      Sample,
      Run
    )
  }
  
  write.csv(detected_mutations, paste(
    result_folder,
    paste("NGS_quality_parameters_", data_type, ".csv", sep = ""),
    sep = "/"
  ))
  
}

analyze_Plasmid_data <- function(data, data_type) {
  
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
  
  groups <- get_groups(data)
  number_of_groups <- nrow(groups)
  path <- paste(result_folder,
                dir_figures,
                data_type,
                dir_Plasmid,
                sep = "/")
  
  for (i in 1:number_of_groups) {
    Sample <- groups[[1]][i]
    Fragment <- groups[[2]][i]
    Run <- groups[[3]][i]

    print(Sample)
    print(Fragment)
    print(Run)
    
    plasmid_filtered <- plasmid_expected_mutations %>%
      filter(fragment == Fragment) 
    
    data_filtered <- data %>%
      filter(sample == Sample , fragment == Fragment , run == Run) %>%
      select(sample,
             Percent_A,
             Percent_B,
             run,
             Variant_level_UMI,
             Variant_UMI,
             pos, 
             number_of_reads,
             `COV.TOTAL`,
             Q_score, 
             Sample_readable) %>%
      full_join(plasmid_filtered, by = c('pos' = 'Position')) %>%
      filter(pos < STR_start | pos > STR_end)
    
    if(100 %in% data_filtered$Percent_A){
      data_filtered <- data_filtered %>% 
        filter(type_annot != "type_b")
    }
    if(100 %in% data_filtered$Percent_B){
      data_filtered <- data_filtered %>% 
        filter(type_annot != "type_a_mut" )
    }
    
    data_filtered <- data_filtered %>%
      mutate(mutation_type = ifelse(
        as.character(Variant_UMI) == TypeA,
        as.character(type_annot) ,
        ifelse(as.character(Variant_UMI) == TypeB,
               "type_b", "undefined")
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
    
    plot_variance_level_per_sample <-
      ggplot(data_filtered,
             aes(x = pos, y = Variant_level_UMI, color = mutation_type)) +
      geom_point() +
      labs(
        x = 'Position compared to reference genome',
        y = 'Variant level',
        title = paste(
          Sample,
          Fragment,
          'Variant level across the whole Amplicon',
          sep = "_"
        )
      ) +
      coord_cartesian(ylim = c(0, 1))+
      scale_y_continuous(breaks = seq(0, 1, by = 0.1))
    
    plot_variance_level_per_sample_zoomed <-
      ggplot(data_filtered,
             aes(x = pos, y = Variant_level_UMI, color = mutation_type)) +
      geom_point() +
      labs(
        x = 'Position compared to reference genome',
        y = 'Variant level',
        title = paste(
          Sample,
          Fragment,
          'Variant level across the whole amplicon',
          sep = "_"
        )
      ) +
      coord_cartesian(ylim = c(0.00, 0.1)) +
      scale_y_continuous(breaks = seq(0, 0.1, by = 0.02))
    
    Percent_A <- median(data_filtered$Percent_A, na.rm = TRUE)  / 100
    Percent_B <- median(data_filtered$Percent_B, na.rm = TRUE) / 100
    Sample_readable <- unique(data_filtered$Sample_readable, nmax = 1)[1]
    # Sample_readable <- str_extract(data_filtered$Sample_readable, "\\d*:\\d*")
    
    print(Percent_A)
    print(Percent_B)
    print(Sample_readable)
    
    density_plot_variant_levels <-
      data_filtered %>% 
      ggplot() +
      geom_density(
        aes(Variant_level_UMI),
        fill = "green",
        color = "grey",
        alpha = 0.4
      ) +
      geom_vline(xintercept = Percent_A, color = "red", alpha = 0.5) +
      geom_vline(xintercept = Percent_B, color = "red", alpha = 0.5) +
      labs(
        x = "relative variant level",
        y = "number of variants per level",
        title = paste('number of variants per level',
                      '(',
                      'TypeA : TypeB', 
                      Sample_readable,
                      Fragment,
                      ')',
                      sep = " ")
      ) +
      annotate(
        geom = "text",
        x = Percent_A,
        y = -0.5,
        label = "Exp."
      ) +
      annotate(
        geom = "text",
        x = Percent_B,
        y = -0.5,
        label = "Exp."
      ) +
      coord_cartesian(xlim = c(0, 1)) +
      scale_x_continuous(breaks = seq(0, 1, by = 0.1))
    
    save_plot(
      density_plot_variant_levels,
      paste(path, 
            dir_density,
            sep = "/"),
      Fragment,
      Sample,
      Run
    )
    
    
    save_plot(
      plot_variance_level_per_sample,
      paste(path,
            dir_variance_per_pos, 
            sep = "/"),
      Fragment,
      Sample,
      Run
    )
    save_plot(
      plot_variance_level_per_sample_zoomed,
      paste(path, dir_variance_per_pos, sep = "/"),
      Fragment,
      Sample,
      Run
    )
    
    write.csv(data_filtered,
              paste(
                path,
                dir_detected_mutations,
                paste(
                  paste(Sample,
                        Fragment,
                        Run,
                        dir_detected_mutations,
                        data_type,
                        sep = "_"),
                  ".csv",
                  sep = ""
                ),
                sep = "/"
              ))
    
  }
  
  write.csv(detected_mutations, paste(
    result_folder,
    paste("Plasmid_quality_parameters_", data_type, ".csv", sep = ""),
    sep = "/"
  ))
  
}

### Analyze data

reads_vs_final_bam_files <- UMI_data %>%
  group_by(sample, fragment, run) %>%
  inner_join(run_info, by = c("run" = "Run_directory")) %>%
  mutate(mean_cov = mean(COV.TOTAL)) %>%
  group_by(mean_cov, number_of_reads, sample, fragment, run, device) %>%
  summarize()

### Number of reads versus number of final consensus sequences
### split by Run or sequencing device

ggsave(
  "num_of_reads_vs_final_consensus_sequences_plot_split_by_run.jpg",
  reads_vs_cons_consensus_sequences_plot(reads_vs_final_bam_files$run),
  path = paste(
    result_folder,
    dir_figures,
    dir_reads_vs_cons_consensus_sequences,
    sep = "/"
  )
)
ggsave(
  "num_of_reads_vs_final_consensus_sequences_plot_split_by_device.jpg",
  reads_vs_cons_consensus_sequences_plot(reads_vs_final_bam_files$device),
  path = paste(
    result_folder,
    dir_figures,
    dir_reads_vs_cons_consensus_sequences,
    sep = "/"
  )
)

write.csv(
  reads_vs_final_bam_files,
  paste(result_folder, "reads_vs_final_bam_files.csv", sep = "/")
)


analyze_NGS_data(NGS_data, "raw")
analyze_NGS_data(NGS_data_filtered, "filtered")
analyze_Plasmid_data(UMI_data %>% filter(grepl('A_B', sample)), "raw")
analyze_Plasmid_data(UMI_data_filtered %>% filter(grepl('A_B', sample)), "filtered")
