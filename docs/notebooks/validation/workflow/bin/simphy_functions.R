# Load libraries
library(phytools)
library(ggtree)
library(tidyverse)
library(rsimpop)
library(truncdist)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(rtracklayer)
library(TreeDist)
library(cowplot)
library(ggsci)
library(pbmcapply)

# Gamma distribution for fitness effects of driver mutations, as estimated by Mitchell et al. 2022
genGammaFitness <- function(shape = 0.47, rate = 34, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  function() rgamma(n = 1, shape = shape, rate = rate)
}

# Generate n random hex hash strings (nchar characters each), for use as mutation IDs. 
random_hash <- function(n, nchar = 12) {
  hex <- c(0:9, letters[1:6])
  apply(matrix(sample(hex, n * nchar, replace = TRUE), nrow = n), 1, paste, collapse = "")
}

# Combine a driver-event table (i.e.,), simpop$events) and a
# driver -> fitness lookup (i.e., simpop$cfg$drivers) into a single data
# frame. Returns NULL if either input is NULL (e.g. a tree that didn't
# originate from run_driver_process_sim()). Called once by get_phylo_object()
# and stored as tree$driver_info.
get_driver_info <- function(events, driver_fitness) {
  if (is.null(events) || is.null(driver_fitness)) return(NULL)

  drivers <- events[events$driverid > 0, c("driverid", "node"), drop = FALSE]
  colnames(drivers) <- c("driver_id", "node")
  drivers$fitness <- driver_fitness$fitness[
    match(drivers$driver_id, driver_fitness$driver)]
  drivers$mutation_id <- random_hash(nrow(drivers))

  drivers[, c("driver_id", "fitness", "node", "mutation_id")]
}

# Create a phylo object from rsimpop output.
get_phylo_object <- function(simpop) {
  tr <- list()
  tr$edge.length <- simpop$edge.length
  tr$edge <- simpop$edge
  tr$tip.label <- simpop$tip.label
  tr$Nnode <- simpop$Nnode
  attr(tr, "class") <- "phylo"
  attr(tr, "order") <- "cladewise"

  tr$driver_info <- get_driver_info(simpop$events, simpop$cfg$drivers)

  tr
}

# Get the node descendant of an edge (or of each edge in a vector of them)
get_edge_descendant <- function(phylo, edge_index) {
  phylo$edge[edge_index, 2]
}

# Get all tips descendant from a node
get_tips_from_node <- function(phylo, node) {
  descendants <- getDescendants(phylo, node)
  tips <- descendants[descendants <= length(phylo$tip.label)]
  tips
}

# Map every edge in a tree to a canonical string key identifying its
# descendant tip-set (its bipartition). This is stable across rerooting,
# ladderizing.
edge_bipartitions <- function(tree) {
  n_tip <- length(tree$tip.label)
  keys <- vapply(seq_len(nrow(tree$edge)), function(i) {
    node <- get_edge_descendant(tree, i)
    tips <- if (node <= n_tip) {
      tree$tip.label[node]
    } else {
      tree$tip.label[get_tips_from_node(tree, node)]
    }
    paste(sort(tips), collapse = "|")
  }, character(1))
  names(keys) <- as.character(seq_len(nrow(tree$edge)))
  keys
}

# Translate an edge index from `from_tree`'s numbering into the
# corresponding edge index in `to_tree`'s numbering, by matching bipartitions. 
# Returns NA if no edge in `to_tree` has the same descendant tip-set.
# `from_keys`/`to_keys` are the precomputed output of edge_bipartitions().
translate_edge_index <- function(edge_idx, from_keys, to_keys) {
  target_key <- from_keys[as.character(edge_idx)]
  match_idx  <- names(to_keys)[to_keys == target_key]
  if (length(match_idx) == 0) return(NA_integer_)
  as.integer(match_idx)
}

# Precompute all-pairs node distances (in edge count, not branch length).
node_edge_distances <- function(tree) {
  unit_tree <- tree
  unit_tree$edge.length <- rep(1, nrow(tree$edge))
  ape::dist.nodes(unit_tree)
}

# Topological distance (in number of edges) between two edges of `tree`,
# defined as the edge-count distance between their child nodes. 
# Pass a precomputed node_edge_distances() matrix via `node_dist` to
# avoid recomputing it on every call.
edge_distance <- function(edge_i, edge_j, tree, node_dist = NULL) {
  if (is.null(node_dist)) node_dist <- node_edge_distances(tree)
  node_dist[get_edge_descendant(tree, edge_i), get_edge_descendant(tree, edge_j)]
}

# Line-graph (edge-adjacency) distance between two edges of `tree`. Unlike
# edge_distance() -- which is 2 for sibling edges, since it goes through
# their shared parent's node -- this is 1 for siblings, and also 1 between
# an edge and its parent edge.
line_graph_distance <- function(edge_i, edge_j, tree, node_dist = NULL) {
  if (edge_i == edge_j) return(0L)
  if (is.null(node_dist)) node_dist <- node_edge_distances(tree)
  pe <- tree$edge[edge_i, 1]; ce <- tree$edge[edge_i, 2]
  pf <- tree$edge[edge_j, 1]; cf <- tree$edge[edge_j, 2]
  1 + min(node_dist[pe, pf], node_dist[pe, cf], node_dist[ce, pf], node_dist[ce, cf])
}

# Line-graph diameter of `tree`: the largest line_graph_distance() between
# any two distinct edges.
line_graph_diameter <- function(tree, node_dist = NULL) {
  if (is.null(node_dist)) node_dist <- node_edge_distances(tree)
  pe <- tree$edge[, 1]; ce <- tree$edge[, 2]
  min_d <- pmin(node_dist[pe, pe], node_dist[pe, ce], node_dist[ce, pe], node_dist[ce, ce])
  diag(min_d) <- NA
  max(min_d, na.rm = TRUE) + 1
}

# Exact transition-probability matrix P(v) = exp(Q*v) 
transition_prob_matrix <- function(Q, v) {
  as.matrix(expm::expm(Q * v))
}

# Read every locus's SCM-assigned edge(s) straight out of a multi_scm() h5
# output file.
#
# RETURNS: data.frame, one row per locus found in the h5 file, with
# `scm_assigned_edges` (comma-separated if SCM placed more than one state
# change for that locus) and `n_edges_assigned`.
read_assigned_edges <- function(h5f_path) {

  fid <- rhdf5::H5Fopen(h5f_path, flags = "H5F_ACC_RDONLY")
  on.exit(rhdf5::H5Fclose(fid), add = TRUE)

  shallow   <- rhdf5::h5ls(fid, recursive = 2)
  ae_groups <- shallow %>%
    dplyr::filter(name == "assigned_edges", otype == "H5I_GROUP") %>%
    dplyr::mutate(path = paste0(ifelse(group == "/", "", group), "/", name))

  if (nrow(ae_groups) == 0) {
    warning("No assigned_edges group found in this h5 file.")
    return(data.frame())
  }

  loc_paths <- lapply(ae_groups$path, function(grp_path) {
    grp <- rhdf5::H5Gopen(fid, grp_path)
    on.exit(rhdf5::H5Gclose(grp), add = TRUE)
    loci <- rhdf5::h5ls(grp, recursive = FALSE)$name
    setNames(paste0(grp_path, "/", loci), loci)
  }) %>% unlist()

  lapply(names(loc_paths), function(loc) {
    scm_vec <- rhdf5::h5read(fid, loc_paths[[loc]])
    scm_idx <- which(scm_vec == 1)

    data.frame(
      locus              = loc,
      scm_assigned_edges = paste(scm_idx, collapse = ","),
      n_edges_assigned   = length(scm_idx),
      stringsAsFactors = FALSE
    )
  }) %>% dplyr::bind_rows()
}

# Compare SCM's inferred mutation-to-edge assignment (multi_scm() h5 output,
# via read_assigned_edges()) against the ground-truth edge(s) each mutation
# was actually emitted on during simulation (the EDGE= INFO tag in the VCF
# written by build_vcf_df()). 
#
# ARGUMENTS:
#   - vcf_df:
#       data.frame from build_vcf_df()/write_vcf_df(), with the ground-truth
#       edge(s) in INFO as "EDGE=<n>" or "EDGE=<n1>,<n2>,...", indexed
#       against `sim_tree`
#   - h5f_path:
#       Path to the h5 file
#   - sim_tree:
#       The ground-truth simulation tree
#   - scm_tree:
#       The tree object output by SCM 
#
# RETURNS: data.frame, one row per locus present in both vcf_df and the h5
# file's assigned_edges group(s)
compare_scm_edges <- function(vcf_df, 
                              h5f_path, 
                              sim_tree, 
                              scm_tree, 
                              cores = 1, 
                              full_df = FALSE) {

  sim_keys <- edge_bipartitions(sim_tree)
  scm_keys <- edge_bipartitions(scm_tree)
  scm_node_dist <- node_edge_distances(scm_tree)

  scm_diameter <- line_graph_diameter(scm_tree, scm_node_dist)

  locus <- paste(vcf_df$CHROM, vcf_df$POS, sep = "_")

  truth_edge_str <- sub(".*EDGE=([0-9,]+).*", "\\1", vcf_df$INFO)
  truth_edges    <- strsplit(truth_edge_str, ",")

  assigned <- read_assigned_edges(h5f_path)
  if (nrow(assigned) == 0) return(data.frame())

  shared <- intersect(locus, assigned$locus)
  if (length(shared) == 0) {
    warning("No loci shared between vcf_df and the h5 file's assigned_edges group(s).")
    return(data.frame())
  }

  mclapply(shared, function(loc) {
    truth_idx <- as.integer(truth_edges[[match(loc, locus)]])
    truth_scm_idx <- vapply(truth_idx, translate_edge_index, integer(1),
                            from_keys = sim_keys, to_keys = scm_keys)
    truth_scm_valid <- truth_scm_idx[!is.na(truth_scm_idx)]

    scm_str <- assigned$scm_assigned_edges[assigned$locus == loc]
    scm_idx <- if (nzchar(scm_str)) as.integer(strsplit(scm_str, ",")[[1]]) else integer(0)

    ## An assignment is correct only if the set of inferred edges is exactly
    ## the set of (translatable) true origin edges
    correct <- length(scm_idx) > 0 && length(truth_scm_valid) > 0 &&
      setequal(scm_idx, truth_scm_valid)

    ## Symmetric nearest-neighbor sum under line-graph distance (see
    ## line_graph_distance()): every inferred edge is charged the distance to
    ## its closest truth edge (spurious/misplaced inferences), AND every
    ## truth edge is charged the distance to its closest inferred edge
    ## (missed origins). Two distinct edges are always >= 1 apart, so this is
    ## 0 exactly when `correct` is TRUE and > 0 otherwise -- including when
    ## SCM recovers only some of the true origins, or adds extras. Line-graph
    ## distance is used rather than node-based edge_distance() because SCM
    ## assigns mutations to edges, not nodes -- e.g. a sibling edge is
    ## genuinely "one step" away from the true edge, not two.
    edge_dist <- if (length(truth_scm_valid) == 0 || length(scm_idx) == 0) {
      NA_integer_
    } else {
      nearest <- function(from, to) {
        sum(vapply(from, function(a) {
          min(vapply(to, function(b) line_graph_distance(a, b, scm_tree, scm_node_dist),
                     numeric(1)))
        }, numeric(1)))
      }
      nearest(scm_idx, truth_scm_valid) + nearest(truth_scm_valid, scm_idx)
    }

    ## Average, per edge considered (inferred + true), what fraction of the
    ## tree's diameter it missed by -- 0 means the two edge sets coincide, 1
    ## means every edge is as far as two edges can possibly be on this tree.
    edge_distance_normalized <- edge_dist /
      ((length(scm_idx) + length(truth_scm_valid)) * scm_diameter)

    if (full_df) {
      data.frame(
        locus                     = loc,
        truth_edge_sim_numbering  = paste(truth_idx, collapse = ","),
        truth_edge_scm_numbering  = paste(truth_scm_idx, collapse = ","),
        scm_assigned_edges        = paste(scm_idx, collapse = ","),
        n_edges_assigned          = length(scm_idx),
        topology_mismatch         = length(truth_scm_valid) == 0,
        correct                   = correct,
        edge_distance             = edge_dist,
        edge_distance_normalized  = edge_distance_normalized,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        locus                     = loc,
        true_edge                 = paste(truth_scm_idx, collapse = ","),
        assigned_edge             = paste(scm_idx, collapse = ","),
        n_assigned                = length(scm_idx),
        correct                   = correct,
        d                         = edge_dist,
        dnorm                     = edge_distance_normalized,
        stringsAsFactors          = FALSE
      )
    }

  }, mc.cores = cores) %>% dplyr::bind_rows()
}

# Function to generate a boolean vector of mutation presence for all tips descendant of a node
create_mut_vector <- function(phylo, node) {
  tips <- get_tips_from_node(phylo, node)
  mut_vector <- rep(0, length(phylo$tip.label))
  mut_vector[tips] <- 1
  mut_vector
}

# Per-mutation metadata columns create_mut_df() adds
MUT_META_COLS <- c("edge", 
                   "mutation_id", 
                   "is_driver", 
                   "selection_coefficient", 
                   "is_ism_violation",
                   "n_origins", 
                   "origin_mutation_ids")

# function to create a data frame of mutation presence.
create_mut_df <- function(tree, driver_info = tree$driver_info) {
  ## Per-edge blocks are collected here and bound together once at the end.
  blocks <- list()

  ## Matching a driver's `node` against
  ## tree$edge[,2] identifies which edge it evolved on.
  has_driver_info <- !is.null(driver_info)

  mutation_id <- character(0)
  is_driver <- logical(0)
  selection_coefficient <- numeric(0)

  for (i in 1:nrow(tree$edge)) {
    n_i <- tree$edge.length[i]
    if (n_i > 0) {
      mut_vector <- c(i,create_mut_vector(tree, get_edge_descendant(tree, i)))
      mut_matrix <- matrix(rep(mut_vector, n_i),
                           nrow = n_i,
                           byrow = TRUE)
      blocks[[length(blocks) + 1L]] <- mut_matrix

      ## A driver event marks the branch it arose on (via its child node)
      driver_flags <- rep(FALSE, n_i)
      sel_coefs <- rep(NA_real_, n_i)
      ids <- character(n_i)
      n_hits <- 0
      if (has_driver_info) {
        child_node <- tree$edge[i, 2]
        hits <- driver_info[driver_info$node == child_node, , drop = FALSE]
        if (nrow(hits) > 0) {
          n_hits <- min(nrow(hits), n_i)  # guard: more drivers matched than mutations to tag on this edge
          driver_flags[seq_len(n_hits)] <- TRUE
          sel_coefs[seq_len(n_hits)] <- hits$fitness[seq_len(n_hits)]
          ids[seq_len(n_hits)] <- hits$mutation_id[seq_len(n_hits)]
        }
      }
      if (n_hits < n_i) {
        ids[(n_hits + 1):n_i] <- random_hash(n_i - n_hits)
      }

      is_driver <- c(is_driver, driver_flags)
      selection_coefficient <- c(selection_coefficient, sel_coefs)
      mutation_id <- c(mutation_id, ids)
    }
  }

  mat <- if (length(blocks) > 0) {
    do.call(rbind, blocks)
  } else {
    matrix(nrow = 0, ncol = length(tree$tip.label) + 1)   # no edge has length > 0
  }
  rm(blocks)
  colnames(mat) <- c("edge",tree$tip.label)
  mat <- as.data.frame(mat)
  ## `edge` is character, not numeric, from the start -- one comma-separated
  ## list of tree-edge indices per mutation. 
  mat$edge <- as.character(mat$edge)
  mat <- mat %>%
    mutate(mutation_id = mutation_id,
           is_driver = is_driver,
           selection_coefficient = selection_coefficient,
           is_ism_violation = FALSE,
           n_origins = 1L,
           .before = 1)
  mat
}

# Introduce infinite-sites-model (ISM) violations into a create_mut_df()
# output by merging groups of otherwise-independent mutation rows into a
# single shared-locus record.
#
# ARGUMENTS:
#   - mut_df:
#       data.frame from create_mut_df()
#   - violation_spec:
#       data.frame with columns `n_origins` (rows to merge per violation,
#       each >= 2) and `count` (how many such violations to introduce at
#       that n_origins)
#   - seed:
#       optional RNG seed for reproducible selection
#   - sel_coef_tol:
#       numeric tolerance for treating two drivers' selection coefficients
#       as the same. The default, 0, requires EXACT equality.
#   - driver_prop:
#       NULL (default) rows are drawn uniformly at random, so a driver is only 
#       included by chance. A number in [0, 1] instead forces exactly
#       round(driver_prop * total violations) violations (chosen at random)
#       to include exactly one driver row.
# RETURNS: a create_mut_df()-shaped data.frame, with an origin_mutation_ids
# column added, and sum(violation_spec$count * (violation_spec$n_origins - 1))
# fewer rows than the input.
introduce_ism_violations <- function(mut_df, 
                                     violation_spec, 
                                     seed = NULL, 
                                     sel_coef_tol = 0,
                                     driver_prop = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(all(c("n_origins", "count") %in% colnames(violation_spec)))
  stopifnot(all(violation_spec$n_origins >= 2), all(violation_spec$count >= 1))
  stopifnot(is.null(driver_prop) || (length(driver_prop) == 1 && driver_prop >= 0 && driver_prop <= 1))

  n_needed <- sum(violation_spec$n_origins * violation_spec$count)
  if (n_needed > nrow(mut_df)) {
    stop(sprintf("violation_spec requires merging %d rows but mut_df only has %d.",
                 n_needed, nrow(mut_df)))
  }

  sample_cols <- setdiff(colnames(mut_df), MUT_META_COLS)

  mut_df$origin_mutation_ids <- mut_df$mutation_id

  used <- logical(nrow(mut_df))
  merged_rows <- list()

  spec_seq <- rep(seq_len(nrow(violation_spec)), violation_spec$count)

  force_driver <- rep(FALSE, length(spec_seq))
  if (!is.null(driver_prop)) {
    n_force <- round(driver_prop * length(spec_seq))
    n_driver_rows <- sum(mut_df$is_driver)
    if (n_force > n_driver_rows) {
      warning(sprintf(paste0("driver_prop = %g asks for %d driver-containing violations but mut_df ",
                             "has only %d driver rows; forcing %d."),
                      driver_prop, n_force, n_driver_rows, n_driver_rows))
      n_force <- n_driver_rows
    }
    force_driver[sample.int(length(spec_seq), n_force)] <- TRUE
  }

  for (i in seq_along(spec_seq)) {
    k <- spec_seq[i]
    n <- violation_spec$n_origins[k]
    available <- which(!used)
    if (length(available) < n) {
      stop("Ran out of unused mutation records while introducing ISM violations.")
    }

    if (!is.null(driver_prop)) {
      ## driver_prop mode: a forced violation takes exactly one unused driver
      avail_drv  <- available[mut_df$is_driver[available]]
      avail_pass <- available[!mut_df$is_driver[available]]
      n_pass_needed <- if (force_driver[i]) n - 1 else n
      if (length(avail_pass) < n_pass_needed) {
        stop("Ran out of unused passenger records while introducing ISM violations.")
      }
      picked <- c(if (force_driver[i]) avail_drv[sample.int(length(avail_drv), 1)],
                  avail_pass[sample.int(length(avail_pass), n_pass_needed)])
    } else {
      ## Reject-and-redraw: draw n rows uniformly at random, accept only if
      ## every driver among them shares the same selection_coefficient.
      max_tries <- 1000
      picked <- NULL
      for (attempt in seq_len(max_tries)) {
        candidate_idx <- sample(available, n)
        driver_coefs <- mut_df$selection_coefficient[candidate_idx][mut_df$is_driver[candidate_idx]]
        if (length(driver_coefs) <= 1 || diff(range(driver_coefs)) <= sel_coef_tol) {
          picked <- candidate_idx
          break
        }
      }
      if (is.null(picked)) {
        stop("Could not find ", n, " mutation records with a consistent driver ",
             "selection coefficient after ", max_tries, " attempts.")
      }
    }

    rows <- mut_df[picked, , drop = FALSE]
    is_driver_merge <- any(rows$is_driver)
    driver_row <- if (is_driver_merge) rows[which(rows$is_driver)[1], ] else NULL

    merged <- rows[1, , drop = FALSE]
    merged[sample_cols] <- as.integer(colSums(rows[sample_cols]) > 0)
    merged$n_origins <- n
    merged$is_ism_violation <- TRUE
    merged$is_driver <- is_driver_merge
    merged$selection_coefficient <- if (is_driver_merge) driver_row$selection_coefficient else NA_real_
    merged$mutation_id <- if (is_driver_merge) driver_row$mutation_id else random_hash(1)
    merged$edge <- paste(rows$edge, collapse = ",")
    merged$origin_mutation_ids <- paste(rows$origin_mutation_ids, collapse = ",")

    merged_rows[[length(merged_rows) + 1]] <- merged
    used[picked] <- TRUE
  }

  survivors <- mut_df[!used, , drop = FALSE]
  dplyr::bind_rows(survivors, merged_rows)
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
                                  dropout_rate   = 0.00, # P(complete allelic dropout | truly heterozygous site)
                                  seed           = NULL) {# optional RNG seed for reproducible draws

  if (!is.null(seed)) set.seed(seed)

  G <- G[, !colnames(G) %in% MUT_META_COLS, drop = FALSE]  # keep only sample columns

  n_mut  <- nrow(G)    # number of mutation/site rows in the ground-truth genotype matrix
  n_samp <- ncol(G)    # number of sample/tip columns

  ## --- Depth model, truncated to the post-QC regime ------------------------
  ## Rationale: rather than simulating depth freely and then filtering out
  ## sites/samples that fail QC (which would require discarding rows/columns
  ## after the fact), we draw directly from the region of parameter space
  ## that would survive filtering. This reflects the assumption that this
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
  conc_floor <- (1 - rho_max) / rho_max # solve rho_max = 1/(1+conc) for conc
  site_conc  <- rtrunc(n_mut, "gamma",
                       a = conc_floor, # lower bound on concentration (upper bound on rho)
                       b = Inf,
                       shape = site_conc_shape,
                       rate = site_conc_rate)
  site_rho   <- 1 / (1 + site_conc) # sanity check: max(site_rho) should be <= rho_max

  ## dropout: models complete allelic dropout at truly heterozygous sites
  is_het  <- G == 1   # TRUE where the ground-truth genotype is heterozygous (mutation present)
  dropout <- is_het & matrix(runif(n_mut * n_samp) < dropout_rate, n_mut, n_samp)

  ## At true hom-ref sites (not heterozygous), alt reads arise only from
  ## sequencing/mapping error -- binomial draw at the error rate.
  AD <- matrix(0L, n_mut, n_samp) # alt-read-count matrix, same shape as DP/G
  colnames(AD) <- colnames(G)
  AD[!is_het] <- rbinom(sum(!is_het), DP[!is_het], error_rate)

  ## Without dropout: alt-read count is a beta-binomial
  ## draw. The beta distribution's shape parameters (conc*0.5, conc*0.5) are
  ## symmetric around VAF = 0.5 (expected for a true heterozygous variant),
  ## with `conc` (site_conc) controlling how tightly VAF clusters around 0.5
  ## versus how much it's allowed to drift due to locus-specific noise.
  het_ok   <- is_het & !dropout
  conc_mat <- matrix(site_conc, n_mut, n_samp)[het_ok]  # broadcast site_conc across samples, then subset
  p_vaf    <- rbeta(sum(het_ok), conc_mat * 0.5, conc_mat * 0.5) 
  AD[het_ok] <- rbinom(sum(het_ok), DP[het_ok], p_vaf)

  ## With dropout: no true alt signal is observable, so
  ## alt reads arise only from background error
  AD[is_het & dropout] <- rbinom(sum(is_het & dropout), DP[is_het & dropout], error_rate)

  ## Return everything needed downstream (DP/AD for VCF construction) plus
  ## the intermediate quantities (useful for diagnostics/sanity checks).
  list(DP = DP, AD = AD, site_factor = site_factor, sample_factor = sample_factor,
       site_conc = site_conc, site_rho = site_rho, dropout = dropout)
}

# GQ/PL/GT from DP, AD (vectorized over full mutation x sample matrices)
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

  ## called GT = argmin PL per site/sample
  GT_idx <- ifelse(PL0 <= PL1 & PL0 <= PL2, 1L, ifelse(PL1 <= PL2, 2L, 3L))
  GT_str <- matrix(c("0/0", "0/1", "1/1")[GT_idx], n_mut, n_samp)
  colnames(GT_str) <- colnames(PL0)

  list(GQ = GQ, PL0 = PL0, PL1 = PL1, PL2 = PL2, GT = GT_str)
}

# Sample mutation genomic loci (locus, ref, alt)
sample_mutation_loci <- function(n,
                                 genome_pkg = "BSgenome.Hsapiens.UCSC.hg38",
                                 chroms     = paste0("chr", c(1:22, "X")),
                                 buffer     = 1000,
                                 oversample = 1.02,
                                 seed       = NULL) {

  if (!is.null(seed)) set.seed(seed)

  library(genome_pkg, character.only = TRUE)
  genome     <- get(genome_pkg)
  chrom_lens <- seqlengths(genome)[chroms]

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
      ## extremely rare at buffer >= 1000
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

# Make VCF like structure.
build_vcf_df <- function(G, DP, AD, gl, edges, mut_table,
                         id = ".", filter = "PASS") {

  ## Pull per-mutation metadata off G (see create_mut_df()) before it gets
  ## stripped down to just sample columns below. `id` is only a fallback for
  ## a G that doesn't carry a mutation_id column (e.g. not from
  ## create_mut_df()); driver status defaults to non-driver in that case.
  mutation_id <- if ("mutation_id" %in% colnames(G)) G$mutation_id else id
  is_driver   <- if ("is_driver" %in% colnames(G)) G$is_driver else FALSE
  selection_coefficient <- if ("selection_coefficient" %in% colnames(G)) {
    G$selection_coefficient
  } else {
    NA_real_
  }
  is_ism_violation <- if ("is_ism_violation" %in% colnames(G)) G$is_ism_violation else FALSE
  n_origins <- if ("n_origins" %in% colnames(G)) G$n_origins else 1L

  G <- G[, !colnames(G) %in% MUT_META_COLS, drop = FALSE]  # keep only sample columns

  sample_names <- colnames(G)
  G_truth <- as.matrix(G)

  n_mut  <- nrow(mut_table)
  n_samp <- ncol(gl$GT)
  stopifnot(nrow(DP) == n_mut, nrow(AD) == n_mut,
            length(sample_names) == n_samp, length(edges) == n_mut)

  fmt_flat <- sprintf(
    "%s:%d:%d,%d:%d:%d,%d,%d:%d",
    gl$GT, DP, DP - AD, AD,
    gl$GQ,
    gl$PL0, gl$PL1, gl$PL2,
    G_truth
  )
  fmt_mat <- matrix(fmt_flat, nrow = n_mut, ncol = n_samp)
  colnames(fmt_mat) <- sample_names

  ## INFO flags
  alt_allele_count <- c("0/0" = 0L, "0/1" = 1L, "1/1" = 2L)[gl$GT]
  dim(alt_allele_count) <- dim(gl$GT)
  AC <- rowSums(alt_allele_count)
  AF <- AC / (2 * n_samp)
  driver_flag <- as.integer(is_driver)
  sel_coef_out <- ifelse(is_driver, selection_coefficient, 0)
  ism_viol_flag <- as.integer(is_ism_violation)

  vcf_df <- data.frame(
    CHROM  = mut_table$chrom,
    POS    = mut_table$pos,
    ID     = mutation_id,
    REF    = mut_table$ref,
    ALT    = mut_table$alt,
    QUAL   = ".",
    FILTER = filter,
    INFO   = sprintf("AC=%d;AF=%.6g;EDGE=%s;DRIVER=%d;S=%.6g;ISM_VIOL=%d;N_ORIGINS=%d",
                     AC, AF, edges, driver_flag, sel_coef_out, ism_viol_flag, n_origins),
    FORMAT = "GT:DP:AD:GQ:PL:TG",
    stringsAsFactors = FALSE
  )

  cbind(vcf_df, as.data.frame(fmt_mat, stringsAsFactors = FALSE))
}

# Write a build_vcf_df() data frame out to a VCF v4.2 file
write_vcf_df <- function(vcf_df, file, tree,
                         chrom_order = paste0("chr", c(1:22, "X")),
                         contig_lengths = "auto",
                         sort = TRUE,
                         source = "simphy_functions.R") {

  ## `tree` is the ground-truth simulation tree the VCF's EDGE= INFO values
  ## are indexed against (see build_vcf_df()) -- required so the tree that
  ## defines "ground truth" travels with the VCF itself
  stopifnot(inherits(tree, "phylo"))

  ## Default: declare a ##contig line for every CHROM actually present, sized
  ## from the same reference genome sample_mutation_loci() draws loci from.
  if (identical(contig_lengths, "auto")) {
    chroms <- intersect(chrom_order, unique(vcf_df$CHROM))
    contig_lengths <- GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38)[chroms]
  }

  meta <- c(
    "##fileformat=VCFv4.2",
    sprintf("##source=%s", source),
    sprintf("##tree=%s", ape::write.tree(tree))
  )

  if (!is.null(contig_lengths)) {
    meta <- c(meta, sprintf("##contig=<ID=%s,length=%d>",
                            names(contig_lengths), contig_lengths))
  }

  meta <- c(meta,
    '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele">',
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency, for each ALT allele">',
    '##INFO=<ID=EDGE,Number=.,Type=Integer,Description="Tree edge(s) on which the mutation arose; more than one indicates an ISM-violating recurrent mutation">',
    '##INFO=<ID=DRIVER,Number=1,Type=Integer,Description="1 if mutation is a driver, 0 otherwise">',
    '##INFO=<ID=S,Number=1,Type=Float,Description="Selection coefficient (fitness); 0 for neutral/passenger mutations">',
    '##INFO=<ID=ISM_VIOL,Number=1,Type=Integer,Description="1 if locus is an infinite-sites-model violation, 0 otherwise">',
    '##INFO=<ID=N_ORIGINS,Number=1,Type=Integer,Description="Number of independent origins of this mutation; 1 unless it is a recurrent/ISM-violating site">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">',
    '##FORMAT=<ID=TG,Number=1,Type=Integer,Description="Ground-truth mutation presence (1) or absence (0) for this sample under the simulated tree, independent of the observed/called GT">'
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

# Rescale a phylo object's branch lengths from raw mutation counts. 
# NOTE: The genome accessibility mask here is assumed to be the "good sites"
scale_branches_to_subs_per_site <- function(tr, mask,
                                            mask_chroms = paste0("chr", c(1:22, "X"))) {

  stopifnot(inherits(tr, "phylo"))

  ## `mask` may be a path/URL to the accessibility mask (BED/BED.gz, e.g. the
  ## 1000 Genomes "strict" or "pilot" whole-genome accessibility mask), or an
  ## already-imported GRanges of accessible intervals.
  mask_gr <- if (is.character(mask)) rtracklayer::import(mask, format = "bed") else mask

  ## 1KGP masks are typically Ensembl-style ("1", "2", ..., "X") rather than
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

# Function estimates the count of expected ISM violations given a 
# phylogenies total branch length (measured in mutations) using birthday paradox
estimate_ism_violations <- function(phylo, L = 3.2e9) {
  ceiling((sum(phylo$edge.length)^2)/(2*L))
}
