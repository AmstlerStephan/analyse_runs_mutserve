library(tidyverse)
library(argparser)

parser <- arg_parser("Commandline parser")
parser <- add_argument(
  parser,
  "--parsed_files",
  nargs = Inf,
  help = "parsed stat files"
)
parser <- add_argument(
  parser,
  "--run",
  help = "Sample sheet to join barcodes and samples"
)
parser <- add_argument(
  parser,
  "--min_qscore",
  help = "Sample sheet to join barcodes and samples"
)
parser <- add_argument(
  parser,
  "--min_read_length",
  help = "Sample sheet to join barcodes and samples"
)

argv <- parse_args(parser)
run <- argv$run
min_read_length <- argv$min_read_length
min_qscore <- argv$min_qscore
parsed_files <- argv$parsed_files

merged_files <- lapply(parsed_files, read_tsv) %>%
  bind_rows()

filename_merged_stats <- paste(run, min_read_length, min_qscore, sep = "_")
write_tsv(merged_files, paste0(filename_merged_stats, ".tsv"))