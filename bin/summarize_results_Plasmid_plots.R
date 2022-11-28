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


UMI_data <-
  read_tsv(UMI_Plasmid_Samples)
plasmid_expected_mutations <-
  read.csv(mutation_classification) %>%
  mutate(
    Position = as.numeric(as.character(Position)),
    Corresponding_Position = as.numeric(as.character(Corresponding_Position))
  )

### functions


get_groups <- function(data) {
  groups <- data %>%
    group_by(sample, fragment, run) %>%
    summarize() %>%
    drop_na()
  return(groups)
}

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

  plasmid_filtered <- plasmid_expected_mutations %>%
    filter(fragment == Fragment) 
  
  data_filtered <- data %>%
    filter(sample == Sample , fragment == Fragment , run == Run) %>%
    select(sample,
           Percent_A,
           Percent_B,
           run,
           variant_level_UMI,
           variant_UMI,
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
      as.character(variant_UMI) == TypeA,
      as.character(type_annot) ,
      ifelse(as.character(variant_UMI) == TypeB,
             "type_b", "undefined")
    )) %>% 
    mutate(variant_level_UMI = coalesce(variant_level_UMI, 0.4))
  
  plot_variance_level_per_sample <-
    ggplot(data_filtered,
           aes(x = pos, y = variant_level_UMI, color = mutation_type)) +
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
           aes(x = pos, y = variant_level_UMI, color = mutation_type)) +
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
  
  print(Percent_A)
  print(Percent_B)
  print(Sample_readable)
  
  density_plot_variant_levels <-
    data_filtered %>% 
    ggplot() +
    geom_density(
      aes(variant_level_UMI),
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

  ggsave(filename =paste(Fragment, Sample, Run, "density_plot_variant_levels", sep = '_'),
         device = "jpg",
         density_plot_variant_levels)
  ggsave(filename =paste(Fragment, Sample, Run, "density_plot_variant_levels", sep = '_'),
         device = "jpg",
         density_plot_variant_levels)
  ggsave(filename =paste(Fragment, Sample, Run, "plot_variance_level_per_sample", sep = '_'),
         device = "jpg",
         plot_variance_level_per_sample)
  ggsave(filename =paste(Fragment, Sample, Run, "plot_variance_level_per_sample_zoomed", sep = '_'),
         device = "jpg",
         plot_variance_level_per_sample_zoomed)
  ggsave(filename =paste(Fragment, Sample, Run, "plot_variance_level_per_sample", sep = '_'),
         device = "jpg",
         plot_variance_level_per_sample)

}
