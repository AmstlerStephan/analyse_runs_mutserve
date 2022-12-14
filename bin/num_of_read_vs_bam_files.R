library(tidyverse)
library(argparser)

parser <- arg_parser("Commandline parser")
parser <- add_argument(
  parser,
  "--umi_merged",
  help = "Merged UMI summary files"
)

argv <- parse_args(parser)
umi_merged <- argv$umi_merged

# umi_merged <- "~/UMI_LPA_KIV2/results_adapted_nf_pipeline_20221130_no_del/umi_summary_file_all_runs/UMI_sequencing_mutserve_merged.tsv"

umi_data <-
  read_tsv(umi_merged)

### Analyze data

reads_vs_final_bam_files <- umi_data %>%
  group_by(sample, fragment, run, is_V14) %>%
  summarize(mean_cov = mean(coverage),
            number_of_reads = number_of_reads,
            min_read_length = min_read_length,
            min_qscore = min_qscore) %>%
  unique()

### Number of reads versus number of final consensus sequences
### split by Run or sequencing device

for ( Run in unique(reads_vs_final_bam_files$run)){
  print(Run)
  
  if( Run == "run3" | Run == "run4"){
    next
  }
  
  reads_vs_final_bam_files_filtered <- reads_vs_final_bam_files %>%
    filter(run == Run)

  correlation = cor.test(x = reads_vs_final_bam_files_filtered$number_of_reads, 
                     y = reads_vs_final_bam_files_filtered$mean_cov)


  reads_vs_final_bam_files_filtered %>% 
    ggplot(aes(x = number_of_reads, y = mean_cov, color = run)) +
    geom_point(show.legend = FALSE) +
    geom_smooth(method = lm, se = FALSE, show.legend = FALSE) +
    theme_light() +
    labs(
      x = "",
      y = ""
    ) +
    annotate("text", x = 15000, y = 100, 
             label = paste0("p: ", round(correlation$p.value, digits = 2), "\n", 
                            "R: ", round(correlation$estimate, digits = 2)),
             size = 7) +
    annotate("text", x = 30000, y = 990, 
             label = str_extract(Run, "run\\d*"), 
             size = 15) +
    coord_cartesian(xlim = c(10, 190000), ylim = c(0, 1000)) +
    scale_x_continuous(breaks = seq(0, 150000, by = 50000)) +
    scale_y_continuous(breaks = seq(0, 1000, by = 200))
  
  ggsave(
    paste("umi_num_of_reads_vs_final_consensus_sequences_plot_split_by_run.jpg", Run, sep = "_"),
    device = "jpg"
  )
  
}


for ( Chemistry in unique(reads_vs_final_bam_files$is_V14)){
  print(Chemistry)
  reads_vs_final_bam_files_filtered <- reads_vs_final_bam_files %>%
    filter(is_V14 == Chemistry) %>% 
    filter(number_of_reads < 200000)

  if(Chemistry){
    Chemistry <- "V14"
  } else {
    Chemistry <- "R9"
  }

  correlation = cor.test(x = reads_vs_final_bam_files_filtered$number_of_reads, 
                     y = reads_vs_final_bam_files_filtered$mean_cov)


  reads_vs_final_bam_files_filtered %>% 
    ggplot(aes(x = number_of_reads, y = mean_cov)) +
    geom_point(show.legend = FALSE) +
    geom_smooth(method = lm, se = FALSE, show.legend = FALSE) +
    theme_light() +
    labs(
      x = "",
      y = ""
    ) +
    annotate("text", x = 15000, y = 100, 
             label = paste0("p: ", round(correlation$p.value, digits = 2), "\n", 
                            "R: ", round(correlation$estimate, digits = 2)),
             size = 7) +
    annotate("text", x = 30000, y = 990, 
             label = Chemistry, 
             size = 15) +
    coord_cartesian(xlim = c(10, 190000), ylim = c(0, 1000)) +
    scale_x_continuous(breaks = seq(0, 150000, by = 50000)) +
    scale_y_continuous(breaks = seq(0, 1000, by = 200))
  
  ggsave(
    paste("umi_num_of_reads_vs_final_consensus_sequences_plot_split_by_run.jpg", Chemistry, sep = "_"),
    device = "jpg"
  )
  
}


write_tsv(
  reads_vs_final_bam_files,
  "umi_reads_vs_final_bam_files.tsv"
)
