#!/usr/bin/env Rscript
library(optparse)

# Read command line arguments
option_list <- list(
  make_option(c("-i", "--input"),
              type = "character",
              help = "Input file(s) to merge (comma-separated to specify multiple files)",
              metavar = "character"),
  make_option(c("-f", "--functions"),
              type = "character",
              help = "Path to functions file (/path/to/SCM.smk.R)",
              metavar = "character"),
  make_option(c("-o", "--output"),
              type = "character",
              help = "Output file name (<output_prefix>.h5)",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (length(commandArgs(trailingOnly = TRUE)) == 0) {
  print_help(OptionParser(option_list = option_list))
  quit(status = 1)
}

if (is.null(opt$input)) {
  stop("At least one file path must be provided using -i or --input.")
}

if (is.null(opt$output)) {
  stop("An output file name must be provided using -o or --output.")
}

source(opt$functions)

# Create h5 file vector and name each element by the second-to-last element
# when delimited by "."
input_files <- unlist(strsplit(opt$input, split = ","))
input_files <- setNames(input_files,
                        sapply(strsplit(input_files, "\\."),
                        function(x) x[length(x) - 1]))

# Copy each h5 file to the output file
merge_h5_files(input_files, opt$output)