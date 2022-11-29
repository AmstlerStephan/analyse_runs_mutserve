library(tidyverse)
library(argparser)

parser <- arg_parser("Commandline parser")
parser <- add_argument(
  parser,
  "--parsed_files",
  nargs = Inf,
  help = "parsed stat files"
)

argv <- parse_args(parser)
parsed_files <- argv$parsed_files

merged_files <- lapply(parsed_files, read_tsv) %>%
  bind_rows()

write_tsv(merged_files, "UMI_sequencing_mutserve_merged.tsv")