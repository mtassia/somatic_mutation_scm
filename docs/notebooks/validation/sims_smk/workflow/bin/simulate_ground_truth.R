#!/usr/bin/env Rscript

## Run one hematopoiesis simulation (rsimpop) and write the ground-truth
## phylogeny and VCF that CellPhy can then be run on to fit genotype
## substitution-rate parameters. 
##
## USAGE
##   Rscript simulate_ground_truth.R --sample_size 100 --age 60 \
##       [--ism] [--ism_multiplier 10] [--seed 1212121] [--out_dir sims]
##
## ARGUMENTS
##   --sample_size      required. Number of cells sampled from the full
##                      population simulation (an outgroup tip, "s1", is added
##                      by rsimpop's sampler, so the tree has sample_size + 1
##                      tips).
##   --age              required. Years of simulated hematopoiesis (nyears).
##   --ism              flag. Introduce ISM violations. When set, the
##                      population is simulated under the discretized gamma
##                      DFE (--dfe_k deviates drawn once, then sampled uniformly
##                      per driver event); otherwise under the continuous
##                      gamma DFE.
##   --ism_multiplier   numeric > 0, default 1. With --ism, introduce
##                      multiplier x E[ISM violations] violations (rounded up),
##                      each with 2 origins, where E is
##                      estimate_ism_violations(tree). Ignored without --ism.
##   --seed             integer, default 1212121.
##   --pop_size         full-population target size, default 10000.
##   --dfe_k            number of discrete fitness deviates for the ISM DFE,
##                      default 10.
##   --driver_prop      optional, in [0, 1]. With --ism, force this proportion
##                      of violations to include one driver mutation (see
##                      introduce_ism_violations()). Default: not used.
##   --mask             REQUIRED. Path to a LOCAL accessibility-mask BED file
##                      (e.g. the 1000 Genomes GRCh38 strict mask, fetched by
##                      the Snakefile's download_mask rule) used to rescale
##                      branch lengths to substitutions/site. URLs are not
##                      accepted, so compute nodes need no network access.
##   --out_dir          output directory, default "ground_truth".
##   --prefix           output file prefix, default built from the arguments.
##
## OUTPUTS  (in --out_dir)
##   <prefix>.vcf       ground-truth VCF (build_vcf_df()/write_vcf_df()); the
##                      EDGE= INFO tag holds each mutation's true origin
##                      edge(s), ISM_VIOL/N_ORIGINS mark ISM violations, and
##                      the FORMAT field TG holds true presence/absence.
##   <prefix>.nwk       ground-truth tree, branch lengths in substitutions per
##                      site -- the --tree input for CellPhy.
##   <prefix>.truth.rds list(tr, tr_mol, params, summary): `tr` is the
##                      un-rescaled phylo whose edge numbering the VCF EDGE=
##                      values index (pass it as `sim_tree` to
##                      compare_scm_edges()), `tr_mol` is the rescaled tree.
##
## DESIGN PARAMETERS FIXED HERE:
##   read simulation: D_bar = 30, min_cov = 10, max_cov = 250,
##   site_sdlog = 0.1, sample_sdlog = 0.1, theta = 8, rho_max = 0.1,
##   site_conc_shape = 4, site_conc_rate = 0.3, error_rate = 0,
##   dropout_rate = 0; driver rate 2e-3 / 365 per cell per day; background
##   mutation rate 16.8 / 365 per day; DFE gamma(shape = 0.47, rate = 34).

## ---- locate and load the shared functions --------------------------------
script_dir <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1) dirname(normalizePath(f)) else getwd()
})
## simphy_functions.R lives with the notebook (docs/notebooks/validation/),
## several directories above this script (sims_smk/workflow/bin/): search the
## script's directory and each parent for it.
search_dirs <- Reduce(function(d, i) dirname(d), seq_len(5), script_dir, accumulate = TRUE)
functions_file <- Filter(file.exists, file.path(search_dirs, "simphy_functions.R"))
if (length(functions_file) == 0) stop("Could not find simphy_functions.R at or above ", script_dir)
suppressMessages(suppressWarnings(source(functions_file[1])))

## ---- argument parsing ----------------------------------------------------
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--sample_size", type = "double", default = NULL,
              help = "REQUIRED. Cells sampled from the full population simulation (an outgroup tip is added, so the tree has sample_size + 1 tips)."),
  make_option("--age", type = "double", default = NULL,
              help = "REQUIRED. Years of simulated hematopoiesis."),
  make_option("--ism", action = "store_true", default = FALSE,
              help = "Introduce ISM violations (discretized gamma DFE; see 'Alternative ISM-violation strategy')."),
  make_option("--ism_multiplier", type = "double", default = 1,
              help = "With --ism: introduce this multiple of the expected number of ISM violations, each with 2 origins [default %default]."),
  make_option("--seed", type = "double", default = 1212121,
              help = "RNG seed [default %default]."),
  make_option("--pop_size", type = "double", default = 1e4,
              help = "Full-population target size [default %default]."),
  make_option("--dfe_k", type = "double", default = 10,
              help = "With --ism: number of discrete fitness deviates drawn from the gamma DFE [default %default]."),
  make_option("--driver_prop", type = "double", default = NULL,
              help = "With --ism: proportion in [0, 1] of violations forced to include one driver mutation [default: not used]."),
  make_option("--mask", type = "character", default = NULL,
              help = "REQUIRED. Path to a local accessibility-mask BED file for rescaling branch lengths to substitutions/site."),
  make_option("--out_dir", type = "character", default = "ground_truth",
              help = "Output directory [default %default]."),
  make_option("--prefix", type = "character", default = NULL,
              help = "Output file prefix [default: built from the arguments].")
)
opt_parser <- OptionParser(
  usage = "usage: %prog --sample_size N --age YEARS [options]",
  option_list = option_list,
  description = "Simulate one hematopoiesis replicate and write the ground-truth VCF and tree for CellPhy."
)
opt <- parse_args(opt_parser)

for (req in c("sample_size", "age", "mask")) {
  if (is.null(opt[[req]])) {
    print_help(opt_parser)
    stop("Missing required argument --", req, call. = FALSE)
  }
}
if (grepl("^[a-zA-Z][a-zA-Z0-9+.-]*://", opt$mask)) {
  stop("--mask must be a local file path, not a URL (", opt$mask,
       "); download it first (e.g. with the Snakefile's download_mask rule).", call. = FALSE)
}
if (!file.exists(opt$mask)) stop("--mask file not found: ", opt$mask, call. = FALSE)
check_whole <- function(x, name) {
  if (is.na(x) || x != round(x)) stop("--", name, " must be a whole number", call. = FALSE)
}
check_whole(opt$sample_size, "sample_size"); check_whole(opt$seed, "seed"); check_whole(opt$dfe_k, "dfe_k")
if (opt$sample_size < 1) stop("--sample_size must be >= 1", call. = FALSE)
if (opt$age <= 0) stop("--age must be > 0", call. = FALSE)
if (opt$ism_multiplier <= 0) stop("--ism_multiplier must be > 0", call. = FALSE)
if (opt$dfe_k < 1) stop("--dfe_k must be >= 1", call. = FALSE)
if (opt$sample_size >= opt$pop_size) stop("--sample_size must be smaller than --pop_size", call. = FALSE)
if (!is.null(opt$driver_prop) && (opt$driver_prop < 0 || opt$driver_prop > 1)) {
  stop("--driver_prop must be in [0, 1]", call. = FALSE)
}
if (!is.null(opt$driver_prop) && !opt$ism) warning("--driver_prop is ignored without --ism")

sample_size    <- opt$sample_size
age            <- opt$age
ism            <- opt$ism
ism_multiplier <- opt$ism_multiplier
seed           <- opt$seed
pop_size       <- opt$pop_size
dfe_k          <- opt$dfe_k
driver_prop    <- opt$driver_prop
mask_path      <- opt$mask
out_dir        <- opt$out_dir
prefix <- if (is.null(opt$prefix)) {
  sprintf("sim_N%d_age%g_ism%d%s_seed%d", as.integer(sample_size), age, as.integer(ism),
          if (ism) sprintf("_x%g", ism_multiplier) else "", as.integer(seed))
} else opt$prefix
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

params <- list(sample_size = sample_size, age = age, ism = ism,
               ism_multiplier = if (ism) ism_multiplier else NA_real_,
               seed = seed, pop_size = pop_size,
               dfe_k = if (ism) dfe_k else NA_integer_,
               driver_prop = driver_prop, mask = mask_path)
message("[", prefix, "] parameters: ",
        paste(names(params), vapply(params, function(x) paste(format(x), collapse = ","), ""),
              sep = "=", collapse = "; "))

## ---- 1. population simulation (rsimpop) ----------------------------------
## Without ISM violations: continuous gamma DFE (Mitchell et al. 2022).
## With ISM violations ("Alternative ISM-violation strategy"): draw K deviates
## from that gamma once, then sample uniformly among them per driver event, so
## independent driver events can share an identical fitness coefficient.
if (ism) {
  set.seed(seed)
  fitness_sample <- rgamma(n = dfe_k, shape = 0.47, rate = 34)
  fitnessFn <- function() sample(fitness_sample, size = 1)
} else {
  fitnessFn <- genGammaFitness(shape = 0.47, rate = 34, seed = seed)
}

## rsimpop prints progress to stdout; keep logs readable when looping this
## script over many replicates.
quiet <- function(expr) {
  invisible(utils::capture.output(res <- force(expr)))
  res
}

initSimPop(seed, bForce = TRUE)
sim <- quiet(run_driver_process_sim(
  simpop                   = NULL,
  initial_division_rate    = 0.1,
  final_division_rate      = 1 / 365,
  target_pop_size          = pop_size,
  nyears                   = age,
  fitness                  = fitnessFn,
  drivers_per_cell_per_day = 2.0e-3 / 365
))

st     <- quiet(get_subsampled_tree(sim, sample_size))
st_mut <- quiet(get_elapsed_time_tree(st, mutrateperdivision = 0, backgroundrate = 16.8 / 365))
tr     <- get_phylo_object(st_mut)

## ---- 2. ground-truth molecular phylogeny ---------------------------------
tr_mol <- scale_branches_to_subs_per_site(tr, mask = mask_path)

## ---- 3. mutation history truth (+ ISM violations) ------------------------
mut_mat <- create_mut_df(tr)

n_ism <- 0L
if (ism) {
  expected_ism <- estimate_ism_violations(tr)
  n_ism <- as.integer(ceiling(ism_multiplier * expected_ism))
  viols <- data.frame(n_origins = 2, count = n_ism)
  mut_mat <- introduce_ism_violations(mut_df = mut_mat,
                                      violation_spec = viols,
                                      seed = seed,
                                      driver_prop = driver_prop)
}

## ---- 4. sequence reads, genotypes, VCF -----------------------------------
sim_reads <- simulate_DP_and_AD(G = mut_mat,
                                D_bar           = 30,
                                min_cov         = 10,
                                max_cov         = 250,
                                site_sdlog      = 0.1,
                                sample_sdlog    = 0.1,
                                theta           = 8,
                                rho_max         = 0.1,
                                site_conc_shape = 4,
                                site_conc_rate  = 0.3,
                                error_rate      = 0,
                                dropout_rate    = 0,
                                seed            = seed)
sim_genotypes <- simulate_GQ_PL_GT(AD = sim_reads$AD, DP = sim_reads$DP)
mut_table     <- sample_mutation_loci(nrow(sim_genotypes$GT), seed = seed)

vcf_df <- build_vcf_df(G = mut_mat,
                       DP = sim_reads$DP,
                       AD = sim_reads$AD,
                       gl = sim_genotypes,
                       edges = mut_mat$edge,
                       mut_table = mut_table)

## Records with no observed alt allele in any sample carry no information
## for SCM or CellPhy, so they are dropped (as in the notebook).
vcf_out <- vcf_df[!grepl("AC=0", vcf_df$INFO, fixed = TRUE), ]

## ---- 5. write outputs ----------------------------------------------------
vcf_file <- file.path(out_dir, paste0(prefix, ".vcf"))
nwk_file <- file.path(out_dir, paste0(prefix, ".nwk"))
rds_file <- file.path(out_dir, paste0(prefix, ".truth.rds"))

write_vcf_df(vcf_out, file = vcf_file, tree = tr)
ape::write.tree(tr_mol, file = nwk_file)

summary <- list(n_tips = length(tr$tip.label),
                n_mutation_records = nrow(mut_mat),
                n_vcf_records = nrow(vcf_out),
                n_ism_violations = sum(mut_mat$is_ism_violation),
                n_ism_violations_requested = n_ism,
                n_driver_events_in_tree = if (is.null(tr$driver_info)) 0L else nrow(tr$driver_info),
                n_distinct_driver_fitness = if (is.null(tr$driver_info)) 0L else length(unique(tr$driver_info$fitness)),
                total_tree_length_mutations = sum(tr$edge.length))
saveRDS(list(tr = tr, tr_mol = tr_mol, params = params, summary = summary), rds_file)

message("[", prefix, "] done: ",
        summary$n_tips, " tips, ", summary$n_vcf_records, " VCF records, ",
        summary$n_ism_violations, " ISM violations -> ", out_dir)
message("CellPhy input: --msa ", vcf_file, " --tree ", nwk_file)
