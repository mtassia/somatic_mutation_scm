#!/usr/bin/env Rscript
library(optparse)

# Read command line arguments
option_list <- list(
  make_option(c("-i", "--input"),
              type = "character",
              help = "H5F path to write tree",
              metavar = "character"),
  make_option(c("-t", "--tree"),
              type = "character",
              help = "Phylogeny used for SCM",
              metavar = "character"),
  make_option(c("-o", "--output"),
              type = "character",
              help = "Output newick file",
              metavar = "character"),
  make_option(c("-f", "--functions"),
              type = "character",
              help = "Path to functions file (/path/to/SCM.smk.R)",
              metavar = "character"),
  make_option(c("-g", "--outgroup"),
              type = "character",
              help = "Outgroup for tree rooting",
              metavar = "character"),
  make_option(c("-m", "--merged"),
              action = "store_true",
              default = FALSE,
              help = paste(
                "Flag to indicate if the input is merged",
                "[default: %default]"
              )),
  make_option(c("-a", "--add_tree"),
              action = "store_true",
              default = FALSE,
              help = paste(
                "Flag to indicate if the tree should be added",
                "to the H5F file [default: %default]"
              ))
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Input tests
if (length(commandArgs(trailingOnly = TRUE)) == 0) {
  print_help(OptionParser(option_list = option_list))
  quit(status = 1)
}
if (is.null(opt$input)) {
  stop("An input file path must be provided using -i or --input.")
}
if (is.null(opt$tree)) {
  stop("A tree file path must be provided using -t or --tree.")
}
if (is.null(opt$functions)) {
  stop("A functions file path must be provided using -f or --functions.")
}
if (is.null(opt$outgroup)) {
  stop("An outgroup must be provided using -g or --outgroup.")
}
if (is.null(opt$output)) {
  stop("An output file name must be provided using -o or --output.")
}

# Load functions
suppressMessages(source(opt$functions))

# Load input
tr <- read.tree(file = opt$tree) %>%
  root(outgroup = opt$outgroup, resolve.root = TRUE) %>%
  drop.tip(opt$outgroup) %>%
  ladderize()

# Generate scm-scaled tree
tr_out <- add_scaled_tree_to_h5f(h5f_path = opt$input,
                                 phylo = tr,
                                 write = opt$add_tree,
                                 merged = opt$merged)

# Write tree to file
write.tree(phy = tr_out,
           file = opt$output)