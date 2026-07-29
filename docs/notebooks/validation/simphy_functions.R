# Load libraries
library(phytools)
library(ggtree)
library(tidyverse)
library(rsimpop)
library(truncdist)
library(BSgenome.Hsapiens.UCSC.hg38)  
library(GenomicRanges)
library(rtracklayer)

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

simulate_DP_and_AD <- function(   G,                     # output from `create_mut_df` (mutation presence matrix)
                                  D_bar          = 20,   # target mean sequencing depth across the whole matrix
                                  min_cov        = 6,    # lower depth bound a site must clear to pass QC (Sequoia-style min_cov)
                                  max_cov        = 250,  # upper depth bound (excludes repeat/mapping-artifact loci)
                                  site_sdlog     = 0.5,  # log-scale spread of per-site coverage bias (GC/mappability)
                                  sample_sdlog   = 0.3,  # log-scale spread of per-sample library depth variation
                                  theta          = 8,    # negative-binomial dispersion for residual (post-QC) depth noise
                                  rho_max        = 0.1,  # beta-binomial overdispersion ceiling a site must clear (Sequoia snv_rho)
                                  site_conc_shape = 4,   # gamma shape for the distribution of per-site VAF concentration
                                  site_conc_rate  = 0.3, # gamma rate for the distribution of per-site VAF concentration
                                  error_rate     = 0.00, # sequencing/mapping error rate (alt reads at true hom-ref sites)
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

## --- GQ/PL/GT from DP, AD (vectorized over full mutation x sample matrices) ---
simulate_GQ_PL_GT <- function(DP, AD, error_rate = 0) {
  p_hom_ref <- error_rate
  p_het     <- 0.5
  p_hom_alt <- 1 - error_rate
  
  n_mut  <- nrow(DP)
  n_samp <- ncol(DP)
  
  ## log10 genotype likelihoods (still needed internally to derive PL/GQ,
  ## not returned directly)
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
  ## and the second-best PL, capped at 99. 
  GQ <- pmin(PL0 + PL1 + PL2 - pmax(PL0, PL1, PL2) - pmin(PL0, PL1, PL2), 99)

  ## called GT = argmin PL per site/sample. Vectorized in place of
  ## apply(..., which.min) for the same reason as GQ above; first-index
  ## tie-break (0/0 over 0/1 over 1/1) matches which.min()'s behavior.
  GT_idx <- ifelse(PL0 <= PL1 & PL0 <= PL2, 1L, ifelse(PL1 <= PL2, 2L, 3L))
  GT_str <- matrix(c("0/0", "0/1", "1/1")[GT_idx], n_mut, n_samp)
  colnames(GT_str) <- colnames(PL0)
  
  list(GQ = GQ, PL0 = PL0, PL1 = PL1, PL2 = PL2, GT = GT_str)
}

## --- Sample mutation genomic loci (locus, ref, alt) ------------------------
sample_mutation_loci <- function(n,
                                 genome_pkg = "BSgenome.Hsapiens.UCSC.hg38",
                                 chroms     = paste0("chr", c(1:22, "X")),
                                 buffer     = 1000,
                                 oversample = 1.02) {
  
  library(genome_pkg, character.only = TRUE)
  genome     <- get(genome_pkg)
  chrom_lens <- seqlengths(genome)[chroms]
  
  ## internal recursive core -- genome/chroms/chrom_lens computed once above,
  ## then just threaded through recursive top-up calls (rare, only fires if
  ## an N/gap region is hit, or a locus collides with one already drawn --
  ## `used` accumulates accepted loci across recursive top-ups so every
  ## returned (chrom, pos) is unique, consistent with an infinite-sites
  ## assumption; duplicate draws otherwise become non-negligible at these
  ## mutation counts via the birthday paradox).
  .sample_core <- function(n, used = character(0)) {
    n_draw <- ceiling(n * oversample)

    ## vectorized chrom + position draws (no per-row loop)
    chrom_draw <- sample(chroms, n_draw, replace = TRUE,
                         prob = chrom_lens / sum(chrom_lens))
    pos_draw   <- as.integer(runif(n_draw, buffer + 1,
                                   chrom_lens[chrom_draw] - buffer))
    locus_draw <- paste(chrom_draw, pos_draw, sep = "_")

    ## single batched genome lookup instead of n_draw separate calls
    gr  <- GRanges(chrom_draw, IRanges(pos_draw, width = 1))
    ref <- as.character(getSeq(genome, gr))

    valid <- ref %in% c("A", "C", "G", "T") &
      !duplicated(locus_draw) & !(locus_draw %in% used)
    if (sum(valid) < n) {
      ## extremely rare at buffer >= 1000; simple one-shot top-up rather than looping
      extra <- .sample_core(n - sum(valid), used = c(used, locus_draw[valid]))
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
  
  .sample_core(n)
}

## Make VCF like structure
build_vcf_df <- function(G, DP, AD, gl, edges,
                         id = ".", filter = "PASS") {

  if (colnames(G)[1] == "edge") {
    G <- G[,-1]  # remove the edge column if present
  }
  
  sample_names <- colnames(G)
  
  n_mut  <- nrow(mut_table)
  n_samp <- ncol(gl$GT)
  stopifnot(nrow(DP) == n_mut, nrow(AD) == n_mut,
            length(sample_names) == n_samp, length(edges) == n_mut)

  fmt_flat <- sprintf(
    "%s:%d:%d,%d:%d:%d,%d,%d",
    gl$GT, DP, DP - AD, AD,
    gl$GQ,
    gl$PL0, gl$PL1, gl$PL2
  )
  fmt_mat <- matrix(fmt_flat, nrow = n_mut, ncol = n_samp)
  colnames(fmt_mat) <- sample_names

  ## AC/AF
  alt_allele_count <- c("0/0" = 0L, "0/1" = 1L, "1/1" = 2L)[gl$GT]
  dim(alt_allele_count) <- dim(gl$GT)
  AC <- rowSums(alt_allele_count)
  AF <- AC / (2 * n_samp)

  vcf_df <- data.frame(
    CHROM  = mut_table$chrom,
    POS    = mut_table$pos,
    ID     = id,
    REF    = mut_table$ref,
    ALT    = mut_table$alt,
    QUAL   = ".",
    FILTER = filter,
    INFO   = sprintf("AC=%d;AF=%.6g;EDGE=%d", AC, AF, edges),
    FORMAT = "GT:DP:AD:GQ:PL",
    stringsAsFactors = FALSE
  )

  cbind(vcf_df, as.data.frame(fmt_mat, stringsAsFactors = FALSE))
}

## Write a build_vcf_df() data frame out to a VCF v4.2 file
write_vcf_df <- function(vcf_df, file,
                         chrom_order = paste0("chr", c(1:22, "X")),
                         contig_lengths = "auto",
                         sort = TRUE,
                         source = "simphy_functions.R") {

  ## Default: declare a ##contig line for every CHROM actually present, sized
  ## from the same reference genome sample_mutation_loci() draws loci from.
  ## htslib-based readers (e.g. the raxml-ng build CellPhy ships) treat any
  ## CHROM missing a ##contig header as a warning-worthy anomaly, so this
  ## should not be left NULL unless you really want an unheadered VCF.
  if (identical(contig_lengths, "auto")) {
    chroms <- intersect(chrom_order, unique(vcf_df$CHROM))
    contig_lengths <- GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38)[chroms]
  }

  meta <- c(
    "##fileformat=VCFv4.2",
    sprintf("##source=%s", source)
  )

  if (!is.null(contig_lengths)) {
    meta <- c(meta, sprintf("##contig=<ID=%s,length=%d>",
                            names(contig_lengths), contig_lengths))
  }

  meta <- c(meta,
    '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele">',
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency, for each ALT allele">',
    '##INFO=<ID=EDGE,Number=1,Type=Integer,Description="Tree edge on which the mutation arose">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">'
  )

  if (sort) {
    ord <- order(factor(vcf_df$CHROM, levels = chrom_order), vcf_df$POS)
    vcf_df <- vcf_df[ord, ]
  }

  fixed_cols  <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
  sample_cols <- setdiff(colnames(vcf_df), fixed_cols)

  header_line <- paste(c("#CHROM", fixed_cols[-1], sample_cols), collapse = "\t")
  body_lines  <- do.call(paste, c(as.list(vcf_df[, c(fixed_cols, sample_cols)]), sep = "\t"))

  writeLines(c(meta, header_line, body_lines), con = file)
  invisible(file)
}

## Rescale a phylo object's branch lengths from raw mutation counts (as held
## by tr <- get_phylo_object(st_mut) %>% drop.tip("s1")) into a molecular
## phylogeny with branch lengths in substitutions/site. NOTE: The genome 
## accessibility mask here is assumed to be the "good sites"
scale_branches_to_subs_per_site <- function(tr, mask,
                                            mask_chroms = paste0("chr", c(1:22, "X"))) {

  stopifnot(inherits(tr, "phylo"))

  ## `mask` may be a path/URL to the accessibility mask (BED/BED.gz, e.g. the
  ## 1000 Genomes "strict" or "pilot" whole-genome accessibility mask), or an
  ## already-imported GRanges of accessible intervals.
  mask_gr <- if (is.character(mask)) rtracklayer::import(mask, format = "bed") else mask

  ## 1000G masks are typically Ensembl-style ("1", "2", ..., "X") rather than
  ## UCSC-style ("chr1", "chr2", ..., "chrX"); harmonise to match mask_chroms.
  if (!any(GenomeInfoDb::seqlevels(mask_gr) %in% mask_chroms)) {
    GenomeInfoDb::seqlevels(mask_gr) <- paste0("chr", GenomeInfoDb::seqlevels(mask_gr))
  }

  mask_gr <- GenomicRanges::reduce(mask_gr)
  mask_gr <- mask_gr[GenomeInfoDb::seqnames(mask_gr) %in% mask_chroms]

  n_accessible <- sum(as.numeric(GenomicRanges::width(mask_gr)))
  stopifnot(n_accessible > 0)

  tr$edge.length <- tr$edge.length / n_accessible
  tr
}

