# Load libraries
library(phytools)
library(ggtree)
library(tidyverse)
library(rsimpop)
library(truncdist)
library(BSgenome.Hsapiens.UCSC.hg38)  
library(GenomicRanges)

# Gamma distribution for fitness effects of driver mutations, as estimated by Mitchell et al. 2022
genGammaFitness <- function(shape = 0.47, rate = 34) {
  function() rgamma(n = 1, shape = shape, rate = rate)
}


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

# Get the node descendant of an edge
get_edge_descendant <- function(phylo, edge_index) {
  edge <- phylo$edge[edge_index, ]
  return(edge[2])
}

# Get all tips descendant from a node
get_tips_from_node <- function(phylo, node) {
  descendants <- getDescendants(phylo, node)
  tips <- descendants[descendants <= length(phylo$tip.label)]
  return(tips)
}

# Function to generate a boolean vector of mutation presence for all tips descendant of a node
create_mut_vector <- function(phylo, node) {
  tips <- get_tips_from_node(phylo, node)
  mut_vector <- rep(0, length(phylo$tip.label))
  mut_vector[tips] <- 1
  return(mut_vector)
}

# function to create a data frame of mutation presence
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

simulate_depth_and_ad <- function(G,
                                  D_bar          = 20,   # target mean sequencing depth across the whole matrix
                                  min_cov        = 6,    # lower depth bound a site must clear to pass QC (Sequoia-style min_cov)
                                  max_cov        = 250,  # upper depth bound (excludes repeat/mapping-artifact loci)
                                  site_sdlog     = 0.5,  # log-scale spread of per-site coverage bias (GC/mappability)
                                  sample_sdlog   = 0.3,  # log-scale spread of per-sample library depth variation
                                  theta          = 8,    # negative-binomial dispersion for residual (post-QC) depth noise
                                  rho_max        = 0.1,  # beta-binomial overdispersion ceiling a site must clear (Sequoia snv_rho)
                                  site_conc_shape = 4,   # gamma shape for the distribution of per-site VAF concentration
                                  site_conc_rate  = 0.3, # gamma rate for the distribution of per-site VAF concentration
                                  error_rate     = 0.002,# sequencing/mapping error rate (alt reads at true hom-ref sites)
                                  dropout_rate   = 0.00) {# P(complete allelic dropout | truly heterozygous site)
  
  if (colnames(G)[1] == "edge") {
    G <- G[,-1]  # remove the edge column if present
  }
  
  n_mut  <- nrow(G)    # number of mutation/site rows in the ground-truth genotype matrix
  n_samp <- ncol(G)    # number of sample/tip columns
  
  ## --- Depth model, truncated to the post-QC regime ------------------------
  ## Rationale: rather than simulating depth freely and then filtering out
  ## sites/samples that fail QC (which would require discarding rows/columns
  ## after the fact), we draw directly from the region of parameter space
  ## that WOULD survive filtering. This reflects the assumption that this
  ## matrix represents already-QC-passed data.
  
  ## site_factor: multiplicative per-site coverage bias (e.g. GC content,
  ## mappability, replication timing). Truncated so that site_factor * D_bar
  ## stays within [min_cov, max_cov]. This guarantees every simulated site would
  ## have cleared the depth-based QC filter.
  site_factor <- rtrunc(n_mut, "lnorm",
                        a = min_cov / D_bar,   # lower truncation point (in units of D_bar)
                        b = max_cov / D_bar,   # upper truncation point
                        meanlog = 0,           # median multiplier = 1 (no systematic bias)
                        sdlog = site_sdlog)    # controls how much sites vary in coverage
  
  ## sample_factor: multiplicative per-sample depth effect. Only lower-truncated,
  ## as a well-sequenced sample isn't penalized. A sample whose mean depth would 
  ## fall below min_cov is excluded.
  sample_factor <- rtrunc(n_samp, "lnorm",
                          a = min_cov / D_bar,
                          b = Inf,
                          meanlog = 0,
                          sdlog = sample_sdlog)
  
  ## Expected depth per (site, sample) cell = D_bar scaled by both the site's
  ## and the sample's multiplicative factors.
  mu <- outer(site_factor, sample_factor) * D_bar
  
  ## Realized depth drawn from a negative binomial (overdispersed relative to
  ## Poisson). 
  DP <- matrix(rnbinom(n_mut * n_samp,
                       mu = mu,
                       size = theta),
               nrow = n_mut)
  colnames(DP) <- colnames(G)
  
  ## --- Site-level beta-binomial concentration, truncated to snv_rho <= rho_max
  ## Rationale: real heterozygous-site VAFs aren't a clean binomial(0.5) draw
  ## -- there's locus-specific overdispersion (mapping issues, local sequence
  ## context, etc.), captured here with a beta-binomial model. The
  ## overdispersion parameter rho = 1/(1+concentration); Sequoia's snv_rho
  ## filter excludes sites whose VAF is too variable across samples (rho too
  ## high / concentration too low). We truncate the concentration prior so
  ## every generated site would already satisfy rho <= rho_max.
  conc_floor <- (1 - rho_max) / rho_max   # solve rho_max = 1/(1+conc) for conc
  site_conc  <- rtrunc(n_mut, "gamma",
                       a = conc_floor,       # lower bound on concentration (upper bound on rho)
                       b = Inf,
                       shape = site_conc_shape,
                       rate = site_conc_rate)
  site_rho   <- 1 / (1 + site_conc)   # sanity check: max(site_rho) should be <= rho_max
  
  ## dropout: models complete allelic dropout at truly heterozygous sites 
  is_het  <- G == 1   # TRUE where the ground-truth genotype is heterozygous (mutation present)
  dropout <- is_het & matrix(runif(n_mut * n_samp) < dropout_rate, n_mut, n_samp)
  
  ## At true hom-ref sites (not heterozygous), alt reads arise only from
  ## sequencing/mapping error -- binomial draw at the error rate.
  AD <- matrix(0L, n_mut, n_samp)   # alt-read-count matrix, same shape as DP/G
  colnames(AD) <- colnames(G)
  AD[!is_het] <- rbinom(sum(!is_het), DP[!is_het], error_rate)
  
  ## At true het sites WITHOUT dropout: alt-read count is a beta-binomial
  ## draw. The beta distribution's shape parameters (conc*0.5, conc*0.5) are
  ## symmetric around VAF = 0.5 (expected for a true heterozygous variant),
  ## with `conc` (site_conc) controlling how tightly VAF clusters around 0.5
  ## versus how much it's allowed to drift due to locus-specific noise.
  het_ok   <- is_het & !dropout
  conc_mat <- matrix(site_conc, n_mut, n_samp)[het_ok]   # broadcast site_conc across samples, then subset
  p_vaf    <- rbeta(sum(het_ok), conc_mat * 0.5, conc_mat * 0.5)
  AD[het_ok] <- rbinom(sum(het_ok), DP[het_ok], p_vaf)
  
  ## At true het sites WITH dropout: no true alt signal is observable, so
  ## alt reads arise only from background error, same as a hom-ref site --
  ## this is what makes dropout "invisible" to a caller (the site looks
  ## exactly like hom-ref, not like a low-confidence het call).
  AD[is_het & dropout] <- rbinom(sum(is_het & dropout), DP[is_het & dropout], error_rate)
  
  ## Return everything needed downstream (DP/AD for VCF construction) plus
  ## the intermediate quantities (useful for diagnostics/sanity checks, e.g.
  ## confirming max(site_rho) <= rho_max, or inspecting which sites/samples
  ## got the most extreme depth-bias draws).
  list(DP = DP, AD = AD, site_factor = site_factor, sample_factor = sample_factor,
       site_conc = site_conc, site_rho = site_rho, dropout = dropout)
}

## --- GQ/PL from DP, AD (vectorized over full mutation x sample matrices) ---
compute_GQ_PL <- function(DP, AD, error_rate = 0.002) {
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