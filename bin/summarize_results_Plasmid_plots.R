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


umi_data <-
  read_tsv(umi_plasmid_samples)
plasmid_expected_mutations <-
  read.csv(mutation_classification) %>%
  mutate(
    position = as.numeric(as.character(Position)),
    corresponding_position = as.numeric(as.character(Corresponding_Position))
  )

### define parameters
STR_start <- 2472
STR_end <- 2505

### create output dirs
umi_density_plot_variant_levels_dir <-
  "umi_density_plot_variant_levels"
umi_plot_variance_level_per_sample_dir <-
  "umi_plot_variance_level_per_sample"
umi_plot_variance_level_per_sample_zoomed_dir <-
  "umi_plot_variance_level_per_sample_zoomed"

dir.create(umi_density_plot_variant_levels_dir)
dir.create(umi_plot_variance_level_per_sample_dir)
dir.create(umi_plot_variance_level_per_sample_zoomed_dir)

groups <- umi_data %>%
  group_by(sample, fragment, run) %>%
  summarize() %>%
  drop_na()

number_of_groups <- nrow(groups)


for (i in 1:number_of_groups) {
  if (number_of_groups < 1) {
    break()
  }
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

  write_tsv(data_filtered, "debugging_joined_table.tsv")

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


  plot_variance_level_per_sample <-
    ggplot(
      data_filtered,
      aes(x = pos, y = variant_level_umi, color = mutation_type)
    ) +
    geom_point() +
    labs(
      x = "Position compared to reference genome",
      y = "Variant level",
      title = paste(
        Sample,
        Fragment,
        "Variant level across the whole Amplicon",
        sep = "_"
      )
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1))

  plot_variance_level_per_sample_zoomed <-
    ggplot(
      data_filtered,
      aes(x = pos, y = variant_level_umi, color = mutation_type)
    ) +
    geom_point() +
    labs(
      x = "Position compared to reference genome",
      y = "Variant level",
      title = paste(
        Sample,
        Fragment,
        "Variant level across the whole amplicon",
        sep = "_"
      )
    ) +
    coord_cartesian(ylim = c(0.00, 0.1)) +
    scale_y_continuous(breaks = seq(0, 0.1, by = 0.02))

  Percent_A <- median(data_filtered$Percent_A, na.rm = TRUE) / 100
  Percent_B <- median(data_filtered$Percent_B, na.rm = TRUE) / 100
  Sample_readable <- unique(data_filtered$Sample_readable, nmax = 1)[1]


  # density_plot_variant_levels <-
  #   data_filtered %>%
  #   ggplot() +
  #   geom_density(
  #     aes(variant_level_umi),
  #     fill = "green",
  #     color = "grey",
  #     alpha = 0.4
  #   ) +
  #   geom_vline(xintercept = Percent_A, color = "red", alpha = 0.5) +
  #   geom_vline(xintercept = Percent_B, color = "red", alpha = 0.5) +
  #   labs(
  #     x = "relative variant level",
  #     y = "number of variants per level",
  #     title = paste("number of variants per level",
  #       "(",
  #       "TypeA : TypeB",
  #       Sample_readable,
  #       Fragment,
  #       ")",
  #       sep = " "
  #     )
  #   ) +
  #   annotate(
  #     geom = "text",
  #     x = Percent_A,
  #     y = -0.5,
  #     label = "Exp."
  #   ) +
  #   annotate(
  #     geom = "text",
  #     x = Percent_B,
  #     y = -0.5,
  #     label = "Exp."
  #   ) +
  #   coord_cartesian(xlim = c(0, 1)) +
  #   scale_x_continuous(breaks = seq(0, 1, by = 0.1))

  # ggsave(
  #   filename =
  #     paste0(
  #       paste(Fragment, Sample, Run, sep = "_"),
  #       ".jpeg"),
  #   path = umi_density_plot_variant_levels_dir,
  #   device = "jpg",
  #   density_plot_variant_levels
  # )

  ggsave(
    filename =
      paste0(
        paste(Fragment, Sample, Run, sep = "_"),
        ".jpeg"),
    path = umi_plot_variance_level_per_sample_dir,
    device = "jpg",
    plot_variance_level_per_sample
  )

  ggsave(
    filename =
      paste0(
        paste(Fragment, Sample, Run, sep = "_"),
        ".jpeg"),
    path = umi_plot_variance_level_per_sample_zoomed_dir,
    device = "jpg",
    plot_variance_level_per_sample_zoomed
  )
}
