library(tidyverse)
library(BlandAltmanLeh)
library(argparser)

parser <- arg_parser("Commandline parser")
parser <- add_argument(
  parser,
  "--NGS_UMI_Samples",
  help = "Merged UMI summary files"
)

argv <- parse_args(parser)
NGS_UMI_Samples <- argv$NGS_UMI_Samples

### load data
NGS_data <-
  read_tsv(NGS_UMI_Samples)

### functions

get_groups <- function(data) {
  groups <- data %>%
    group_by(sample, fragment, run) %>%
    summarize() %>%
    drop_na()
  return(groups)
}

create_bland_altman <- function(data, path, Fragment, Sample, Run) {
  jpeg(
    file = paste(Fragment, Sample, Run, "bland_altman.jpg", sep = '_'),
    width = 10,
    height = 10,
    units = "in",
    res = 300
  )
  
  bland_stats <-
    bland.altman.stats(data$variant_level_NGS,
                       data$variant_level_UMI)
  
  print(
    bland.altman.plot(
      data$variant_level_NGS,
      data$variant_level_UMI,
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

groups <- get_groups(data)
number_of_groups <- nrow(groups)

  # For loop to create data per Sample and Run
for (i in 1:number_of_groups) {
  Sample <- groups[[1]][i]
  Fragment <- groups[[2]][i]
  Run <- groups[[3]][i]

  print(Sample)
  print(Fragment)
  
  data_filtered <- data %>%
    filter(sample == Sample , fragment == Fragment)

  r_squared <-
    data_filtered %>% lm(variant_level_UMI ~ variant_level_NGS, data = .)
  r_squared <- summary(r_squared)$r.squared
  
  comparison_variant_levels_ngs_umi <- data_filtered %>%
    ggplot(aes(x = variant_level_UMI, y = variant_level_NGS)) +
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
    geom_point(aes(y = variant_level_NGS, color = 'NGS')) +
    geom_point(aes(y = variant_level_UMI, color = 'UMI')) +
    labs(
      x = 'position relative to reference sequence',
      y = 'relative variant Level',
      title = paste(
        Sample,
        Fragment,
        'variant levels per position and Sequencing technology',
        sep = "_"
      )
    ) +
  scale_x_continuous(breaks = seq(0, 5200, by = 200))
  
  density_plot_variant_levels <-
    data_filtered %>% 
    ggplot() +
    geom_density(
      aes(variant_level_UMI),
      fill = "green",
      color = "grey",
      alpha = 0.4
    ) +
    geom_density(
      aes(variant_level_NGS),
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
  
  ggsave(filename =paste(Fragment, Sample, Run, "comparison_variant_levels_ngs_umi", sep = '_'),
         device = "jpg",
         comparison_variant_levels_ngs_umi)
  
  ggsave(filename =paste(Fragment, Sample, Run, "comparison_variant_levels_ngs_umi_per_position", sep = '_'),
         device = "jpg",
         comparison_variant_levels_ngs_umi_per_position)
  
  ggsave(filename =paste(Fragment, Sample, Run, "density_plot_variant_levels", sep = '_'),
         device = "jpg",
         density_plot_variant_levels)
 
  
}