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

umi_comparison_variant_levels_ngs_dir <-
  "umi_comparison_variant_levels_ngs"
umi_comparison_variant_levels_ngs_per_position_dir <-
  "umi_comparison_variant_levels_ngs_per_position"
umi_density_plot_variant_levels_dir <-
  "umi_density_plot_variant_levels"
umi_bland_altman_dir <-
  "umi_bland_altman"

dir.create(umi_comparison_variant_levels_ngs_dir)
dir.create(umi_comparison_variant_levels_ngs_per_position_dir)
dir.create(umi_density_plot_variant_levels_dir)
dir.create(umi_bland_altman_dir)

### functions

create_bland_altman <- function(data, path, Fragment, Sample) {
  jpeg(
    file = paste(
      path,
      paste(Fragment, Sample, "bland_altman.jpg", sep = "_"),
      sep = "/"
    ),
    width = 10,
    height = 10,
    units = "in",
    res = 300
  )

  bland_stats <-
    bland.altman.stats(
      data$variant_level_ngs,
      data$variant_level_umi
    )

  print(
    bland.altman.plot(
      data$variant_level_ngs,
      data$variant_level_umi,
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

groups <- ngs_data %>%
  group_by(sample, fragment) %>%
  summarize() %>%
  drop_na()

number_of_groups <- nrow(groups)

# For loop to create data per Sample and Run
for (i in 1:number_of_groups) {
  if (number_of_groups < 1) {
    break()
  }

  Sample <- groups[[1]][i]
  Fragment <- groups[[2]][i]

  print(Sample)
  print(Fragment)

  data_filtered <- ngs_data %>%
    filter(sample == Sample, fragment == Fragment)

  r_squared <-
    data_filtered %>% lm(variant_level_umi ~ variant_level_ngs, data = .)
  r_squared <- summary(r_squared)$r.squared

  comparison_variant_levels_ngs_umi <- data_filtered %>%
    ggplot(aes(x = variant_level_umi, y = variant_level_ngs)) +
    geom_abline() +
    geom_smooth(method = "lm", show.legend = TRUE) +
    geom_point(position = "jitter") +
    labs(
      x = "relative variance level UMI",
      y = "relative variance level NGS",
      title = paste(
        Sample,
        Fragment,
        "Variant levels of both Sequencing technology",
        sep = "_"
      )
    ) +
    annotate(
      "text",
      x = 0.8,
      y = 0.05,
      label = paste0("R Squared = ", round(r_squared, digits = 3))
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1))

  # Compare variance levels of different Sequencing methods

  comparison_variant_levels_ngs_umi_per_position <-
    data_filtered %>%
    ggplot(aes(x = position)) +
    geom_point(aes(y = variant_level_ngs, color = "NGS")) +
    geom_point(aes(y = variant_level_umi, color = "UMI")) +
    labs(
      x = "position relative to reference sequence",
      y = "relative variant Level",
      title = paste(
        Sample,
        Fragment,
        "variant levels per position and Sequencing technology",
        sep = "_"
      )
    ) +
    scale_x_continuous(breaks = seq(0, 5200, by = 200))

  density_plot_variant_levels <-
    data_filtered %>%
    ggplot() +
    geom_density(
      aes(variant_level_umi),
      fill = "green",
      color = "grey",
      alpha = 0.4
    ) +
    geom_density(
      aes(variant_level_ngs),
      fill = "blue",
      color = "grey",
      alpha = 0.4
    ) +
    labs(
      x = "relative variant level",
      y = "number of variants per level",
      title = paste(Sample,
        Fragment,
        "number of variants per level of UMI vs NGS",
        sep = " "
      )
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

  # Bland Altman Plot NGS vs UMI variant levels
  create_bland_altman(
    data_filtered,
    umi_bland_altman_dir,
    Fragment,
    Sample
  )

  ggsave(
    filename =
      paste0(
        paste(Fragment, Sample, sep = "_"),
        ".jpeg"),
    path = umi_comparison_variant_levels_ngs_dir,
    device = "jpg",
    comparison_variant_levels_ngs_umi
  )

  ggsave(
    filename =
      paste0(
        paste(Fragment, Sample, sep = "_"),
        ".jpeg"),
    path = umi_comparison_variant_levels_ngs_per_position_dir,
    device = "jpg",
    comparison_variant_levels_ngs_umi_per_position
  )

  ggsave(
    filename =
      paste0(
        paste(Fragment, Sample, sep = "_"),
        ".jpeg"),
    path = umi_density_plot_variant_levels_dir,
    device = "jpg",
    density_plot_variant_levels
  )
}
