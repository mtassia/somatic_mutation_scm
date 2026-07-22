# Simulate a phylogeny that reflects the processes underlying hematopoietic evolution.
# From phylogeny, simulate mutations with a known mutation rate and substitution matrix
# Test performance of SCM on this ground truth dataset

# Load libraries
library(phytools)
library(tidyverse)
library(rsimpop)
library(truncdist)
library(BSgenome.Hsapiens.UCSC.hg38)  
library(GenomicRanges)

#### SIMULATE HEMATOPOIETIC HISTORY ####

SEED <- 1212121
initSimPop(SEED, bForce = TRUE)

## --- Selection coefficient distribution -----------------------------------
## rsimpop's internal "fitness" is a per-division coefficient s; for adult HSCs
## dividing ~once/year, exp(s)-1 approximates the annual clonal growth rate,
## which is how CH driver fitness is usually reported in the literature
## (e.g. Watson et al. 2020 eLife; Fabre et al. 2022 Nature — UK Biobank CH
## drivers show a roughly exponential-tailed distribution of annual growth
## rates, with a detection floor around 0.05-0.08/yr and a long tail out to
## ~0.3-0.5/yr for the strongest DNMT3A/TET2/SF3B1-class drivers).
## Treat these threshold/rate values as a starting point to recalibrate
## against whichever published fitness table you're targeting.
genGammaFitness <- function(shape, rate) {
  function() rgamma(n = 1, shape = shape, rate = rate)
}
fitnessFn <- genGammaFitness(shape = 0.47, rate = 34) # parameters from Mitchell et al. 2022

## sanity-check the implied annual growth-rate distribution before committing
## compute to a full run
hist(sapply(1:1e4, function(x) exp(fitnessFn()) - 1),
     breaks = 100, xlim = c(0, 1),
     xlab = "Selection coefficient per year", main = "")

## --- Fixed-population, multi-driver simulation -----------------------------
dps <- run_driver_process_sim(
  simpop                   = NULL,
  initial_division_rate    = 0.1,           # fast early developmental growth
  final_division_rate      = 1/(2*190),     # adult symmetric division rate
  target_pop_size          = 1e5,           # fixed, not time-varying
  nyears                   = 50,            # age at sampling
  fitness                  = fitnessFn,     # called once per driver acquisition
  drivers_per_cell_per_day = 2.0e-3 / 365   # ABC estimate from Mitchell et al. 2022
)

## --- Sample and convert to mutation-scaled branch lengths ------------------
st <- get_subsampled_tree(dps, 100)
st_mut <- get_elapsed_time_tree(st, 
                                mutrateperdivision = 0, 
                                backgroundrate = 16.8/365)

plot_tree_events(st_mut, cex.label = 0, legpos = NULL)

## --- Extract data from simpop objects in prep for simulating variants-------

# Create a phylo object from rsimpop output
get_phylo_object <- function(simpop){
  tr <- list()
  tr$edge.length <- simpop$edge.length
  tr$edge <- simpop$edge
  tr$tip.label <- simpop$tip.label
  tr$Nnode <- simpop$Nnode
  attr(tr, "class") <- "phylo"
  attr(tr, "order") <- "cladewise"
  
  return(tr)
}

# Function to get the node descendant of an edge
get_edge_descendant <- function(phylo, edge_index) {
  edge <- phylo$edge[edge_index, ]
  return(edge[2])
}

# Function to get all tips descendant from a node
get_tips_from_node <- function(phylo, node) {
  descendants <- getDescendants(phylo, node)
  tips <- descendants[descendants <= length(phylo$tip.label)]
  return(tips)
}

# Function to generate a boolean vector of mutation presnce 
#for all tips descendent of a node
create_mut_vector <- function(phylo, node) {
  tips <- get_tips_from_node(phylo, node)
  mut_vector <- rep(0, length(phylo$tip.label))
  mut_vector[tips] <- 1
  return(mut_vector)
}

# function to create a dataframe of mutation presence.
create_mut_df <- function(tree) {
  mat <- matrix(nrow = 0, ncol = length(tree$tip.label) + 1)
  
  for (i in 1:nrow(tree$edge)) {
    if (tree$edge.length[i] > 0) {
      mut_vector <- c(i,create_mut_vector(tree, get_edge_descendant(tree, i)))
      mut_matrix <- matrix(rep(mut_vector, tree$edge.length[i]), 
                           nrow = tree$edge.length[i], 
                           byrow = TRUE)
      mat <- rbind(mat, mut_matrix)
    }
  }  
  
  colnames(mat) <- c("edge",tree$tip.label)
  mat <- as.data.frame(mat)
  return(mat)
}

## --- Create mutation dataframe from phylo object --------------------------------
tr <- get_phylo_object(st_mut) %>% drop.tip("s1")
G <- create_mut_df(tr)
G_edges <- G[, 1]  # store edge column for reference
G <- G[, -1]  # remove edge column, leaving only mutation presence/absence

## --- Depth model, truncated to the post-QC regime -------------------------

n_mut <- nrow(G); n_samp <- ncol(G)

D_bar   <- 20 # mean depth
min_cov <- 6
max_cov <- 250

## site_factor drawn from a log-normal TRUNCATED so implied mean depth
## (site_factor * D_bar, before per-sample variation) stays within
## [min_cov, max_cov] — i.e. this locus would have survived depth filtering
site_factor <- rtrunc(n_mut, "lnorm", 
                      a = min_cov / D_bar, 
                      b = max_cov / D_bar,
                      meanlog = 0, 
                      sdlog = 0.5)

## sample_factor similarly excludes samples that would fail a minimum
## mean-coverage bar (the crude stand-in for the mixmodel clonality check --
## a sample failing on depth this badly wouldn't pass clonality either)
sample_factor <- rtrunc(n_samp, "lnorm", 
                        a = min_cov / D_bar, 
                        b = Inf,
                        meanlog = 0, 
                        sdlog = 0.3)

theta <- 8   # this is "clean" residual noise
mu <- outer(site_factor, sample_factor) * D_bar # Expected 
DP <- matrix(rnbinom(n_mut * n_samp, 
                     mu = mu, 
                     size = theta), 
             nrow = n_mut)

## --- Site-level beta-binomial concentration, truncated to snv_rho <= 0.1 --
## rho = 1/(1+conc), so rho <= 0.1  <=>  conc >= 9
site_conc <- rtrunc(n_mut, "gamma",
                    a = 9, 
                    b = Inf, 
                    shape = 4, 
                    rate = 0.3)
site_rho  <- 1 / (1 + site_conc)   # sanity check: max(site_rho) should be <= 0.1

## --- Allele counts: no germline injection, no non-clonal sample -----------
error_rate   <- 0.002
dropout_rate <- 0.00

AD <- matrix(0L, n_mut, n_samp)
is_het  <- G == 1
dropout <- is_het & matrix(runif(n_mut * n_samp) < dropout_rate, n_mut, n_samp)

AD[!is_het] <- rbinom(sum(!is_het), DP[!is_het], error_rate)

het_ok   <- is_het & !dropout
conc_mat <- matrix(site_conc, n_mut, n_samp)[het_ok]
p_vaf    <- rbeta(sum(het_ok), conc_mat * 0.5, conc_mat * 0.5)
AD[het_ok] <- rbinom(sum(het_ok), DP[het_ok], p_vaf)

AD[is_het & dropout] <- rbinom(sum(is_het & dropout), DP[is_het & dropout], error_rate)

## --- GQ/PL from DP, AD (vectorized over full mutation x sample matrices) ---
compute_GL_PL <- function(DP, AD, error_rate = 0.002) {
  p_hom_ref <- error_rate
  p_het     <- 0.5
  p_hom_alt <- 1 - error_rate
  
  ## log10 genotype likelihoods 
  ll <- function(p) dbinom(AD, DP, p, log = TRUE) / log(10)
  GL0 <- ll(p_hom_ref)   # 0/0
  GL1 <- ll(p_het)       # 0/1
  GL2 <- ll(p_hom_alt)   # 1/1
  
  ## PL: phred-scaled, normalized so the best genotype = 0, capped at 255
  m   <- pmax(GL0, GL1, GL2)
  PL0 <- pmin(round(-10 * (GL0 - m)), 255)
  PL1 <- pmin(round(-10 * (GL1 - m)), 255)
  PL2 <- pmin(round(-10 * (GL2 - m)), 255)
  
  ## GQ: GATK-style genotype quality = difference between the best PL (0)
  ## and the second-best PL, capped at 99
  PL_stack <- array(c(PL0, PL1, PL2), dim = c(dim(DP), 3))
  PL_sorted <- apply(PL_stack, c(1, 2), sort)   # 3 x n_mut x n_samp, ascending
  GQ <- pmin(PL_sorted[2, , ], 99)
  
  list(GQ = GQ, PL0 = PL0, PL1 = PL1, PL2 = PL2)
}

gl <- compute_GL_PL(DP, AD, error_rate = error_rate)


## called GT = argmin PL per site/sample
PL_stack <- array(c(gl$PL0, gl$PL1, gl$PL2), dim = c(n_mut, n_samp, 3))
GT_idx   <- apply(PL_stack, c(1, 2), which.min)   # 1=0/0, 2=0/1, 3=1/1
GT_str   <- matrix(c("0/0", "0/1", "1/1")[GT_idx], n_mut, n_samp)

## Simulate drawing from a genome for mutation loci
genome     <- BSgenome.Hsapiens.UCSC.hg38
chroms     <- paste0("chr", c(1:22, "X"))
chrom_lens <- seqlengths(genome)[chroms]

sample_mutations <- function(n, genome, chroms, chrom_lens,
                             buffer = 1000, oversample = 1.02) {
  
  n_draw <- ceiling(n * oversample)
  
  ## vectorized chrom + position draws (no per-row loop)
  chrom_draw <- sample(chroms, n_draw, replace = TRUE,
                       prob = chrom_lens / sum(chrom_lens))
  pos_draw   <- as.integer(runif(n_draw, buffer + 1,
                                 chrom_lens[chrom_draw] - buffer))
  
  ## single batched genome lookup instead of n_draw separate calls
  gr  <- GRanges(chrom_draw, IRanges(pos_draw, width = 1))
  ref <- as.character(getSeq(genome, gr))
  
  valid <- ref %in% c("A", "C", "G", "T")
  if (sum(valid) < n) {
    ## extremely rare at buffer >= 1000; simple one-shot top-up rather than looping
    extra <- sample_mutations(n - sum(valid), genome, chroms, chrom_lens, buffer, oversample)
    chrom_out <- c(chrom_draw[valid], extra$chrom)[1:n]
    pos_out   <- c(pos_draw[valid],   extra$pos)[1:n]
    ref_out   <- c(ref[valid],        extra$ref)[1:n]
  } else {
    idx <- which(valid)[1:n]
    chrom_out <- chrom_draw[idx]
    pos_out   <- pos_draw[idx]
    ref_out   <- ref[idx]
  }
  
  ## vectorized ALT assignment via lookup table (no per-row sample())
  bases    <- c("A", "C", "G", "T")
  alt_opts <- sapply(bases, function(b) setdiff(bases, b))  # 3 x 4 matrix, columns named by ref base
  alt_out  <- alt_opts[cbind(sample.int(3, n, replace = TRUE),
                             match(ref_out, bases))]
  
  data.frame(chrom = chrom_out, pos = pos_out, ref = ref_out, alt = alt_out,
             stringsAsFactors = FALSE)
}

mut_table <- sample_mutations(nrow(G), genome, chroms, chrom_lens)

## Make VCF like structure

build_vcf_df <- function(mut_table, GT_str, DP, AD, gl, sample_names, edges,
                         id = ".", filter = "PASS") {
  
  n_mut  <- nrow(mut_table)
  n_samp <- ncol(GT_str)
  stopifnot(nrow(GT_str) == n_mut, nrow(DP) == n_mut, nrow(AD) == n_mut,
            length(sample_names) == n_samp, length(edges) == n_mut)
  
  ## fully vectorized per-sample FORMAT string: GT:DP:AD:GQ:PL
  fmt_flat <- sprintf(
    "%s:%d:%d,%d:%d:%d,%d,%d",
    GT_str, DP, DP - AD, AD,
    gl$GQ,
    gl$PL0, gl$PL1, gl$PL2
  )
  fmt_mat <- matrix(fmt_flat, nrow = n_mut, ncol = n_samp)
  colnames(fmt_mat) <- sample_names
  
  vcf_df <- data.frame(
    CHROM  = mut_table$chrom,
    POS    = mut_table$pos,
    ID     = id,
    REF    = mut_table$ref,
    ALT    = mut_table$alt,
    QUAL   = ".",
    FILTER = filter,
    INFO   = sprintf("EDGE=%d", edges),
    FORMAT = "GT:DP:AD:GQ:PL",
    stringsAsFactors = FALSE
  )
  
  cbind(vcf_df, as.data.frame(fmt_mat, stringsAsFactors = FALSE))
}

vcf_df <- build_vcf_df(mut_table, GT_str, DP, AD, gl,
                       sample_names = colnames(G),
                       edges = G_edges)
