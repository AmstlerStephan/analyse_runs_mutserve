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

umi_data <-
  read_tsv(umi_merged)

### functions

reads_vs_cons_consensus_sequences_plot <- function(split_by) {
  reads_vs_final_bam_files %>%
    ggplot(aes(x = number_of_reads, y = mean_cov, color = split_by)) +
    geom_point() +
    geom_smooth(method = lm, se = FALSE)
}

### Analyze data

reads_vs_final_bam_files <- umi_data %>%
  mutate(mean_cov = mean(COV.TOTAL)) %>%
  group_by(mean_cov, number_of_reads, sample, fragment, run, is_V14) %>%
  summarize()

### Number of reads versus number of final consensus sequences
### split by Run or sequencing device

reads_vs_cons_consensus_sequences_plot(reads_vs_final_bam_files$run)

ggsave(
  "umi_num_of_reads_vs_final_consensus_sequences_plot_split_by_run.jpg",
  device = "jpg"
)

reads_vs_cons_consensus_sequences_plot(reads_vs_final_bam_files$is_V14)

ggsave(
  "umi_num_of_reads_vs_final_consensus_sequences_plot_split_by_kit.jpg",
  device = "jpg"
)

write_tsv(
  reads_vs_final_bam_files,
  "umi_reads_vs_final_bam_files.tsv"
)
