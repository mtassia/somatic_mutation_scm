#!/usr/bin/env Rscript

cat(paste0(paste0(format(Sys.time(), "[%D %H:%M:%S]"),
                  " Loading libraries...\n")))
library(optparse)

# Read command line arguments
option_list <- list(
  make_option(c("-v", "--vcf"),
              type = "character",
              help = "VCF file (<vcf>.gz)",
              metavar = "character"),
  make_option(c("-f", "--functions"),
              type = "character",
              help = "Path to functions file (/path/to/SCM.smk.R)",
              metavar = "character"),
  make_option(c("-t", "--tree"),
              type = "character",
              help = "Tree file (<tree>.nwk)",
              metavar = "character"),
  make_option(c("-m", "--model"),
              type = "character",
              help = "Model file (<model>.model)",
              metavar = "character"),
  make_option(c("-p", "--prefix"),
              type = "character",
              help = "Output prefix for <output_prefix>.h5",
              metavar = "character"),
  make_option(c("-o", "--outgroup"),
              type = "character",
              help = "Outgroup for tree rooting",
              metavar = "character"),
  make_option(c("-n", "--threads"),
              type = "integer",
              default = 1,
              help = "Number of threads to use [default: %default]",
              metavar = "integer"),
  make_option(c("-i", "--iterations"),
              type = "integer",
              default = 100,
              help = "Number of scm iterations to perform [default: %default]",
              metavar = "integer"),
  make_option(c("-c", "--chromosome"),
              type = "character",
              default = "all",
              help = "Target chromosome(s) [default: %default]",
              metavar = "character"),
  make_option(c("-x", "--overwrite"),
              action = "store_true",
              default = FALSE,
              help = "Overwrite existing h5 file [default: %default]"),
  make_option(c("-d", "--dryrun"),
              action = "store_true",
              default = FALSE,
              help = "Dry run scm without writing output [default: %default]"),
  make_option(c("-b", "--brute"),
              action = "store_true",
              default = FALSE,
              help = "Skip singleton LRT and run SCM on all somatic SNPs through SCM [default: %default]")
)

if (length(commandArgs(trailingOnly = TRUE)) == 0) {
  print_help(OptionParser(option_list = option_list))
  quit(status = 1)
}

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

required_args <- c(opt$vcf, opt$tree, opt$model, opt$prefix, opt$outgroup)
if (any(sapply(required_args, is.null))) {
    print_help(opt_parser)
    stop("All required arguments (vcf, tree, model, prefix, and outgroup) must be provided.", call. = FALSE)
}

#### Load functions ####
cat(paste0(paste0(format(Sys.time(), "[%D %H:%M:%S]"),
                  " Loading functions...\n")))
source(opt$functions)

##### Data import #####
#Read tree
cat(paste0(paste0(format(Sys.time(), "[%D %H:%M:%S]"),
                  " Reading tree...\n")))
tr <- read.tree(file = opt$tree) %>%
  root(outgroup = opt$outgroup, resolve.root = T) %>%
  drop.tip(opt$outgroup) %>%
  ladderize()

#Read vcf
cat(paste0(paste0(format(Sys.time(), "[%D %H:%M:%S]"),
                  " Reading vcf...\n")))
dat <- read.vcfR(file = opt$vcf, verbose = FALSE)

#Compose genotype state list
cat(paste0(paste0(format(Sys.time(), "[%D %H:%M:%S]"),
                  " Generating genotype state list...\n")))
gt_list.snp <- compile_gt_states.snp(dat)

#Prep Q matrices
cat(paste0(paste0(format(Sys.time(), "[%D %H:%M:%S]"),
                  " Reading model...\n")))
snp.q <- read_cellphy_model(bestModel_path = opt$model)

##### Run SCM #####
# Run multi_scm
cat(paste0(paste0(format(Sys.time(), "[%D %H:%M:%S]"),
                  " Beginning SCM...\n")))
multi_scm(gt_state_list = gt_list.snp,
          tree = tr,
          Q = snp.q,
          scm_its = opt$iterations,
          cores = opt$threads,
          h5f_path = paste0(opt$prefix, ".h5"),
          chr = opt$chromosome,
          overwrite = opt$overwrite,
          dryrun = opt$dryrun,
          brute = opt$brute)
