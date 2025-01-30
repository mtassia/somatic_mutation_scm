####* INSTALL/LOAD LIBRARIES *####
cran_packages <- c("BiocManager", "tidyverse", "data.table",
                    "pbapply", "pbmcapply", "vcfR",
                    "ape", "phytools", "cowplot",
                    "RColorBrewer", "snow", "markophylo",
                    "ggtext")

for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, quietly = TRUE)
  } else {
    library(pkg, character.only = TRUE, quietly = TRUE)
  }
}

biconductor_packages <- c("ggtree")
for (pkg in biconductor_packages) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
    library(pkg, quietly = TRUE)
  } else {
    library(pkg, character.only = TRUE, quietly = TRUE)
  }
}

####* TREE PARSING *####

#* Collapse tree edges with user-defined support thresholds *#
collapse_poor_supported_edges  <-  function(phylo, support_threshold_to_keep) {
  # ARGUMENTS:
  #   - phylo:
  #       "phylo" tree object; e.g., output from read.tree(<nwk>)
  #   - support_threshold_to_keep:
  #       Node support values below this threshold will be cause the
  #       descendent branch to be distribute equally to all daughter lineages.

  nodes <- phylo$node %>% as.integer()
  badnodes <- which(nodes < support_threshold_to_keep) + length(phylo$tip.label)
  collapsed <- CollapseNode(phylo, badnodes)
  collapsed$node.label <- NULL
  return(collapsed)
}

####* DATA PREP *####

#* Compute genotype priors from phred-scaled genotype likelihood (PL->GP)
unphred <- function(phredscore) {
  # ARGUMENTS:
  #   - phredscore:
  #       integer; genotype likelihood (PL)

  prob <- 10^(-phredscore / 10)
  prob_norm <- prob / sum(prob)
  return(prob_norm)
}

#* Read VCF, generate a genotype prior matrix for each SNP, and store all
#* in a named list.
#TODO: Consolidate this and compile_gt_states.indel into a single function
compile_gt_states.snp <- function(vcf, mat_list = TRUE) {
  # ARGUMENTS:
  #   - vcf:
  #       VCF imported as "vcfR" object using vcfR::read.vcfR(<input>).
  #       All entries must be biallelic.
  #   - mat_list:
  #       Boolean; specify whether output is a named list of datatables
  #       per variant locus (TRUE) or a single long datatable
  
  # Read VCF with vcfR
  print("Reading vcf data...")

  # Convert vcfR object to tidy list of tibbles
  tidy_vcf <- vcfR2tidy(vcf,
                        gt_column_prepend = "",
                        allele.sep = ",",
                        format_fields = c("GT", "GQ", "PL"))

  # Join genotype fields (FORMAT) and site-wise fixed fields by ChromKey,
  # an arbitrary chromosome ID that is indexed by order of occurence. 
  df <- left_join(tidy_vcf$gt,
                  tidy_vcf$fix %>%
                    select(ChromKey, CHROM) %>%
                    unique(),
                  by = "ChromKey") %>%
    select(CHROM, POS, Indiv, GT, GQ, PL)

  # Add REF, ALT, AC, and AF fields to working dataframe
  df <- left_join(df,
                  tidy_vcf$fix %>%
                    select(CHROM, POS, REF, ALT, AC, AF),
                  by = c("CHROM", "POS"))

  # Compute genotype priors from genotype likelihood (PL) fields
  print("Computing genotype prior probabilities...")

  # 1. Add field to dataframe that converts genotype code (e.g., 0/1) to
  #     a string reflecting the diploid genotype state (e.g., AG).
  # 2. Reorder dataframe columns
  # 3. Separate PL field into three columns, one for each diploid state
  # 4. Compute genotype priors from each PL field
  #TODO: Genotype prior compute can be accelerated by extracting PL fields
  #TODO: and applying a vector function.
  df <- df %>%
    mutate(locus = paste(CHROM, POS, sep = "_"),
           gt_str = case_when(GT %in% c("0/0", "0|0") ~ paste0(REF, REF),
                            GT %in% c("1/0", "1|0", 
                                      "0/1", "0|1") ~ paste0(REF, ALT),
                            GT %in% c("1/1", "1|1") ~ paste0(ALT, ALT),
                            is.na(GT) ~ NA)) %>%
    select(locus, CHROM, POS, REF, ALT, AC, AF,
            Indiv, GT, gt_str, GQ, PL) %>%
    separate(PL, sep = ",", into=c("PL_ref", "PL_het", "PL_alt")) %>%
    mutate(GQ = as.integer(GQ),
           PL_ref = as.integer(PL_ref),
           PL_het = as.integer(PL_het),
           PL_alt = as.integer(PL_alt),
           AF = as.numeric(AF),
           AC = as.integer(AC)) %>%
    rowwise() %>%
    mutate(Pref = unphred(c(PL_ref, PL_het, PL_alt))[1],
           Phet = unphred(c(PL_ref, PL_het, PL_alt))[2],
           Palt = unphred(c(PL_ref, PL_het, PL_alt))[3])
  
  # Rearrange data frame and convert to data table for faster access
  df <- df %>%
    arrange(CHROM, POS) %>%
    select(locus, REF, ALT, gt_str, 
            Indiv, Pref, Phet, Palt) %>%
    as.data.table()

  # If mat_list argument is set to FALSE, return long data table. Otherwise,
  # return a list of data tables where each is named as <CHROM>_<POS>.
  if (mat_list == FALSE){
    return(df)
  } else {
    print("Assembling genotype prior state matrices for SNPs...")

    # Subset datatable to SNP loci using REF and ALT character string matching
    df.snp <- df %>%
      group_by(locus) %>%
      filter(REF %in% c("A", "C", "G", "T") & ALT %in% c("A", "C", "G", "T"))
    
    # Sort contents of gt_str field such that it's alphabetical and compatible
    # with the SNP Q matrix generated by read_cellphy_model.
    df.snp$gt_str <- str_split(df.snp$gt_str, pattern = "") %>%
      pblapply(., function(x) {
        sort(x) %>% paste0(., collapse = "")
      }) %>%
      unlist()

    # 1. Split datatable by locus (<CHROM_POS>) and retain only
    #     the locus, Indiv, and 3 observed genotype prior columns
    # 2. Add all genotype states that are unobserved and set
    #     genotype prior to 0.
    gt_list.snp <- df.snp %>%
      group_by(locus) %>%
      group_split() %>%
      pblapply(., function(x){
        ref <- head(x, n = 1) %>% pull(REF)
        alt <- head(x, n = 1) %>% pull(ALT)
        colnames(x) <- c("locus","REF","ALT","gt_str","Indiv",
                       paste0(c(ref, ref), collapse = ""),
                       paste0(sort(c(ref, alt)), collapse = ""),
                       paste0(c(alt, alt), collapse = ""))
        return(x %>%
                 select(!c(REF, ALT, gt_str)) %>%
                 as.data.table())
      }) %>%
      pblapply(., function(x) {
        cols_to_add <- setdiff(c("AA", "CC", "GG", "TT", "AC",
                                  "AG", "AT", "CG", "CT", "GT"),
                                colnames(x)[3:5])
        x[, cols_to_add] = 0

        x <- x %>%
          select(locus,
                 Indiv,
                 c("AA", "CC", "GG", "TT", "AC",
                    "AG", "AT", "CG", "CT", "GT"))

        return(x)
      })

    # Name each datatable in list by <CHROM>_<POS>
    names(gt_list.snp) <- df.snp %>% 
                          group_keys() %>%
                          pull()

    # Return list of named genotype prior datatables for SNPs
    return(gt_list.snp)
  }
}

#* Read VCF, generate a genotype prior matrix for each INDEL, and store all
#* in a named list.
#TODO: Consolidate this into the compile_gt_states.snp function
#! See compile_gt_states.snp for notes on function elements; this function
#! will be removed in future updates.
compile_gt_states.indel <- function(vcf, mat_list = TRUE) {
  print("Reading vcf data...")
  tidy_vcf <- vcfR2tidy(vcf,
                        gt_column_prepend = "",
                        allele.sep = ",",
                        format_fields = c("GT", "GQ", "PL"))
  df <- left_join(tidy_vcf$gt,
                  tidy_vcf$fix %>%
                    select(ChromKey, CHROM) %>%
                    unique(),
                  by="ChromKey") %>%
        select(CHROM, POS, Indiv, GT, GQ, PL)
  df <- left_join(df,
                  tidy_vcf$fix %>%
                    select(CHROM, POS, REF, ALT, AC, AF),
                  by = c("CHROM", "POS"))

  print("Computing genotype probabilities...")
  df <- df %>%
    mutate(locus = paste(CHROM, POS, sep = "_"),
           gt_str = case_when(GT %in% c("0/0", "0|0") ~ paste0(REF, REF),
                            GT %in% c("1/0", "1|0", "0/1", "0|1") ~ paste0(REF, ALT),
                            GT %in% c("1/1", "1|1") ~ paste0(ALT, ALT),
                            is.na(GT) ~ NA)) %>%
    select(locus, CHROM, POS, REF, ALT, AC,
            AF, Indiv, GT, gt_str, GQ, PL) %>%
    separate(PL, sep = ",", into = c("PL_ref", "PL_het", "PL_alt")) %>%
    mutate(GQ = as.integer(GQ),
           PL_ref = as.integer(PL_ref),
           PL_het = as.integer(PL_het),
           PL_alt = as.integer(PL_alt),
           AF = as.numeric(AF),
           AC = as.integer(AC)) %>%
    rowwise() %>%
    mutate(Pref = unphred(c(PL_ref, PL_het, PL_alt))[1],
           Phet = unphred(c(PL_ref, PL_het, PL_alt))[2],
           Palt = unphred(c(PL_ref, PL_het, PL_alt))[3])
  df <- df %>%
    arrange(CHROM, POS) %>%
    select(locus ,REF, ALT, gt_str,
            Indiv, Pref, Phet, Palt) %>%
    as.data.table()
  if (mat_list == FALSE) {
    return(df) #Return a long data frame with genotype probabilities for each sample per variant
  } else {
    print("Assembling genotype state matrices for indels...")
    df.indel <- df %>%
      group_by(locus) %>%
      filter(nchar(REF) > 1 | nchar(ALT) > 1)

    gt_list.indel <- df.indel %>%
      group_by(locus) %>%
      group_split() %>%
      pblapply(., function(x){
        ref <- head(x, n = 1) %>% pull(REF)
        alt <- head(x, n = 1) %>% pull(ALT)
        x <- x %>%
          select(!c(REF,ALT,gt_str))
        colnames(x) <- c("locus", "Indiv", "REF", "HET", "ALT")
        return(as.data.table(x))
      })
    names(gt_list.indel) <- df.indel %>% group_keys() %>% pull()
    return(gt_list.indel)
  }
}

#* Read GT10 substitution model from `CellPhy` *.bestModel output and
#* store as a matrix compatible with, e.g., phytools::plot.Qmatrix()
read_cellphy_model <- function(bestModel_path) {
  # ARGUMENTS:
  #   - bestModel_path:
  #       bestModel text output from `cellphy` run with
  #       best-fit GT10 substitution model.

  RateCats <- c("Zero", "AC", "AG", "AT", "CG", "CT", "GT")
  States <- c("AA", "CC", "GG", "TT", "AC", "AG", "AT", "CG", "CT", "GT")

  #Read bestmodel file
  txt <- read_file(file = as.character(bestModel_path))

  #Digest rates matrix from bestmodel file
  Rates <- txt %>%
    str_split(pattern = "\\+") %>%
    purrr::map(., 1) %>%
    unlist() %>%
    str_remove_all(pattern = "GT10") %>%
    str_remove_all(pattern = "[{}]") %>%
    str_split(pattern = "/") %>%
    unlist() %>%
    as.numeric()
  names(Rates) <- RateCats

  alpha <- Rates["AC"]
  beta <- Rates["AG"]
  gamma <- Rates["AT"]
  kappa <- Rates["CG"]
  lambda <- Rates["CT"]
  mu <- Rates["GT"]

  #Digest state frequencies matrix from bestmodel file
  Freq <- txt %>%
    str_split(pattern = "\\+") %>%
    purrr::map(., 2) %>%
    unlist() %>%
    str_remove_all(pattern = "FU") %>%
    gsub("\\}.*", "", .) %>%
    str_remove_all(pattern = "\\{") %>%
    str_split(pattern = "/") %>%
    unlist() %>%
    as.numeric()
  names(Freq) <- States

  #Populate Q matrix columns in the following order:
  #"AA","CC","GG","TT","AC","AG","AT","CG","CT","GT"
  AAvec <- c(NA,
          Rates["Zero"],
          Rates["Zero"],
          Rates["Zero"],
          alpha * Freq["AA"],
          beta * Freq["AA"],
          gamma * Freq["AA"],
          Rates["Zero"],
          Rates["Zero"],
          Rates["Zero"])
  AAvec <- unname(AAvec)

  CCvec <- c(Rates["Zero"],
          NA,
          Rates["Zero"],
          Rates["Zero"],
          alpha * Freq["CC"],
          Rates["Zero"],
          Rates["Zero"],
          kappa * Freq["CC"],
          lambda * Freq["CC"],
          Rates["Zero"])
  CCvec <- unname(CCvec)

  GGvec <- c(Rates["Zero"],
          Rates["Zero"],
          NA,
          Rates["Zero"],
          Rates["Zero"],
          beta * Freq["GG"],
          Rates["Zero"],
          kappa * Freq["GG"],
          Rates["Zero"],
          mu * Freq["GG"])
  GGvec <- unname(GGvec)

  TTvec <- c(Rates["Zero"],
          Rates["Zero"],
          Rates["Zero"],
          NA,
          Rates["Zero"],
          Rates["Zero"],
          gamma * Freq["TT"],
          Rates["Zero"],
          lambda * Freq["TT"],
          mu * Freq["TT"])
  TTvec <- unname(TTvec)

  ACvec <- c(alpha * Freq["AC"],
          alpha * Freq["AC"],
          Rates["Zero"],
          Rates["Zero"],
          NA,
          kappa * Freq["AC"],
          lambda * Freq["AC"],
          beta * Freq["AC"],
          gamma * Freq["AC"],
          Rates["Zero"])
  ACvec <- unname(ACvec)

  AGvec <- c(beta * Freq["AG"],
          Rates["Zero"],
          beta * Freq["AG"],
          Rates["Zero"],
          kappa * Freq["AG"],
          NA,
          mu * Freq["AG"],
          alpha * Freq["AG"],
          Rates["Zero"],
          gamma * Freq["AG"])
  AGvec <- unname(AGvec)

  ATvec <- c(gamma * Freq["AT"],
          Rates["Zero"],
          Rates["Zero"],
          gamma * Freq["AT"],
          lambda * Freq["AT"],
          mu * Freq["AT"],
          NA,
          Rates["Zero"],
          alpha * Freq["AT"],
          beta * Freq["AT"])
  ATvec <- unname(ATvec)
  CGvec <- c(Rates["Zero"],
          kappa * Freq["CG"],
          kappa * Freq["CG"],
          Rates["Zero"],
          beta * Freq["CG"],
          alpha * Freq["CG"],
          Rates["Zero"],
          NA,
          mu * Freq["CG"],
          lambda * Freq["CG"])
  CGvec <- unname(CGvec)

  CTvec <- c(Rates["Zero"],
          lambda * Freq["CT"],
          Rates["Zero"],
          lambda * Freq["CT"],
          gamma * Freq["CT"],
          Rates["Zero"],
          alpha * Freq["CT"],
          mu * Freq["CT"],
          NA,
          kappa * Freq["CT"])
  CTvec <- unname(CTvec)

  GTvec <- c(Rates["Zero"],
          Rates["Zero"],
          mu * Freq["GT"],
          mu * Freq["GT"],
          Rates["Zero"],
          gamma * Freq["GT"],
          beta * Freq["GT"],
          lambda * Freq["GT"],
          kappa * Freq["GT"],
          NA)
  GTvec <- unname(GTvec)

  #Assemble Q matrix and compute diagonal
  Q <- matrix(c(AAvec, CCvec, GGvec, TTvec, ACvec,
                AGvec, ATvec, CGvec, CTvec, GTvec),
              nrow = 10, ncol = 10, byrow = FALSE)
  Q[1, 1] <- sum(Q[1, ], na.rm = TRUE) * -1
  Q[2, 2] <- sum(Q[2, ], na.rm = TRUE) * -1
  Q[3, 3] <- sum(Q[3, ], na.rm = TRUE) * -1
  Q[4, 4] <- sum(Q[4, ], na.rm = TRUE) * -1
  Q[5, 5] <- sum(Q[5, ], na.rm = TRUE) * -1
  Q[6, 6] <- sum(Q[6, ], na.rm = TRUE) * -1
  Q[7, 7] <- sum(Q[7, ], na.rm = TRUE) * -1
  Q[8, 8] <- sum(Q[8, ], na.rm = TRUE) * -1
  Q[9, 9] <- sum(Q[9, ], na.rm = TRUE) * -1
  Q[10, 10] <- sum(Q[10, ], na.rm = TRUE) * -1

  #Name dimensions
  colnames(Q) <- States
  rownames(Q) <- States
  
  #Return Q matrix
  return(Q)
}

#* Generate INDEL substitution model using observed INDEL genotype states
#* and somatic phylogeny.
generate_indel_model <- function(indel_state_list, tree, cores = 4){
  # ARGUMENTS:
  #   - indel_state_list:
  #       Output from compile_gt_states.indel(mat_list = TRUE)
  #   - tree:
  #       Phylogenetic tree caturing somatic relationships among all
  #       samples represented in indel_state_list
  #   - cores:
  #       Number of parallel computes for `pbmclapply`
  
  print("Constructing indel character-state matrix...")

  # For each locus in indel_state_list:
  #   1. Subset datatable to sample,P(REF),P(HET),P(ALT)
  #   2. Convert datatable to dataframe
  #   3. Assign each sample a genotype state (i.e., REF, HET, ALT)
  #       based on the maximum prior
  #   4. Output a named vector containing the character states per sample
  #   5. Compile all named vectors into a single dataframe where
  #       column names are samples, row names are loci, and values are 
  #       genotype states
  indel.states <- pbmclapply(indel_state_list, mc.cores = cores, function(x) {
    # subset columns
    x <- x %>% 
      select(2:5) %>%
      as.data.frame()
    
    # assign each sample a genotype state (REF, HET, ALT) based one
    # the maximum prior probability. If all priors are equal (e.g., 
    # missing data), to NA. If REF and HET are equal, set state to
    # REF.
    x <- x %>%
      rowwise() %>%
      mutate(
        max_value = max(REF, HET, ALT),
        state = case_when(
          sum(REF == max_value, HET == max_value, ALT == max_value) > 1 ~ NA_character_,
          sum(REF == max_value, HET == max_value) > 1 ~ "REF",
          REF == max_value ~ "REF",
          HET == max_value ~ "HET",
          ALT == max_value ~ "ALT"
        )
      ) %>%
      ungroup() %>%
      select(-max_value) %>% 
      as.data.frame() %>%
      column_to_rownames(var = "Indiv")
    
    # Create a named character state vector, where each element is
    # the indel genotype state, and the name is the sample
    char.c <- x$state
    names(char.c) <- rownames(x)
    return(char.c)
  }) %>% 
    as.data.frame() %>%
    .[ , colSums(is.na(.)) == 0] %>%
    t()
  
  print("Fitting parameters...")

  # Use markophylo::estimaterates to infer a maximum likelihood
  # indel state substitution matrix (REF<->HET<->ALT) given the
  # observed state matrix generated above and the phylogeny. 
  rate <- estimaterates(usertree = tree,
                      userphyl = indel.states,
                      matchtipstodata = TRUE,
                      alphabet = c("REF", "HET", "ALT"),
                      modelmat = matrix(c(NA, 1, 0,
                                          1, NA, 1,
                                          0, 1, NA),
                                        nrow = 3,
                                        ncol = 3),
                      rootprob = "maxlik",
                      numhessian = F)
  
  # Generate INDEL Q matrix using maximum likelihood substitution
  # rate parameter.
  #! Assumes the transition probability between REF<->HET == ALT<->HET.
  #! Does not include character state frequencies, as states are not
  #!    comparable between loci.
  indel.q <- matrix(c((-1 * rate$results$wop$par[1]), rate$results$wop$par[1], 0.001,
                        rate$results$wop$par[1], (-2 * (rate$results$wop$par[1])), rate$results$wop$par[1],
                        0.001, rate$results$wop$par[1], (-1 * rate$results$wop$par[1])),
                    nrow = 3,
                    ncol = 3)
  rownames(indel.q) <- c("REF", "HET", "ALT")
  colnames(indel.q) <- c("REF", "HET", "ALT")

  # Return INDEL Q matrix
  return(indel.q)
}

####* STOCHASTIC CHARACTER MAPPING *####

# Run phytools::make.simmap() for a single locus
runSCM_single <- function(x, tree, gt_state_list, 
                          Qmat, reduced = F, reps = 100,
                          root_state = "equal", cores = 1) {
  # ARGUMENTS:
  #   - x:
  #       Name of target locus in gt_state_list
  #   - tree:
  #       Phylo object with tip for each of the samples in gt_state_list$<x>
  #   - gt_state_list
  #       List of genotype state priors generated by
  #       compile_gt_states.snp(mat_list = T) or
  #       compile_gt_states.indel(mat_list = T)
  #   - Qmat:
  #       Substitution matrix generated by read_cellphy_model or
  #       generate_indel_model
  #   - reduced:
  #       Boolean; subset genotype state matrix (gt_state_list$<x>) to
  #       only the observed genotype states
  #   - reps:
  #       Number of SIMMAP simulations
  #   - root_state:
  #       Prior distribution for the root node of tree
  #   - cores:
  #       Number of cores to use for make.simmap()

  # Report error if more than one locus string is supplied to x
  if (length(x) > 1){
    return(print("ERROR: >1 variant used for input"))
  }
  
  # Grab focal genotype-prior matrix from gt_state_list and name
  # rows by sample.
  gt_matrix <- gt_state_list[[x]] %>% 
    select(!locus) %>%
    column_to_rownames(var = "Indiv")
  
  # If reduced == TRUE, subset gt_matrix and Qto only the observed states
  if (reduced == T) {
    gt_matrix  <-  gt_matrix %>% 
      select(where(~sum(.) != 0))
    Qmat <- Qmat[colnames(gt_matrix), colnames(gt_matrix)]
  }
  
  # Prep cores and print run info
  cl <- makeSOCKcluster(rep("localhost", cores))
  print(paste0("Performing stochastic character mapping on ", x, "."))
  print(paste0("Cores: ", cores))
  print(paste0("Replicates: ", reps))
  print(paste0("Reduced Qmatrix: ", reduced))
  
  # Run phytools::make.simmap() using all cores specified in <cl>
  scm <- clusterApply(cl,
                      x = replicate(cores,
                                    as.matrix(gt_matrix),
                                    simplify = FALSE),
                      fun = make.simmap,
                      tree = tree,
                      Q = Qmat,
                      pi = root_state,
                      nsim = as.integer(round(reps / cores)))
  scm <- do.call("c", scm)
  
  # Reclass simmap if not properly classed by output above
  if(!("multiSimmap" %in% class(scm))) {
    class(scm) <- c("multiSimmap",class(scm))
  }
  
  # Clear cluster
  stopCluster(cl)

  # Return multiimmap object
  return(scm)
}

#* Output a summary for SNP SCM, optionally plot
summarise_scm.snp <- function(multiSimmap, PPthreshold = 0.95, plot = FALSE,
                              legend = TRUE , title = NA) {
  # ARGUMENTS:
  #   - multiSimmap:
  #       Output from runSCM_single()
  #   - PPthreshold:
  #       Genotype state posterior threshold for confident state call
  #   - plot:
  #       Boolean; set TRUE to generate plot summary of multiSimmap
  #   - legend:
  #       Boolean; set TRUE to add genotype state legend to plot
  #   - title:
  #       String; title for plot.

  print("Summarising SCMs...")

  # Load data from multiSimmap (ancestral state estimates and tree)
  scm_summary <- summary(multiSimmap) #! Longest step
  scm_summary <- scm_summary$ace %>% as.data.frame()
  tree <- multiSimmap[[1]]
  
  # Assign each sample it's posterior consensus state using
  # the maximum posterior estimate. If all states < PPthreshold,
  # set constate to NA.
  scm_summary$constate <- apply(scm_summary,
                                1,
                                pull_consensus_state,
                                threshold = PPthreshold)
  # assign rowname to the node column
  scm_summary$node <- rownames(scm_summary)

  # Assign tip labels their node numbers according to the node numbering
  # conventions from ape & phytools. Made consistent with the input tree.
  scm_summary[(tree$Nnode + 1) : nrow(scm_summary), ]$node <- 1:(length(tree$tip.label))
  
  # Order dataframe by node number and pull vector or consensus states
  genotype_states <- arrange(scm_summary, as.integer(node)) %>% 
    pull(constate)
  
  print("Preparing output...")

  # Prepare genotype posteriors per node/tip data frame:
  #   - Ordered sequentially by node/tip number
  #   - All genotype states are present as columns
  #   - Values are genotype posteriors for each node/tip
  out <- scm_summary %>% 
    arrange(as.integer(node)) %>%
    select(!c(node, constate))
  rownames(out) <- NULL
  cols_to_add <- setdiff(c("AA", "CC", "GG", "TT", "AC",
                            "AG", "AT", "CG", "CT", "GT"), 
                          colnames(out))
  out[,cols_to_add] <- 0
  out <- out %>% 
    select("AA", "CC", "GG", "TT", "AC",
            "AG", "AT", "CG", "CT", "GT")

  # Create a vector of consensus state posteriors (in node order)
  constate_posteriors <- apply(out, 1, max)

  # Create an output dataframe for the 95% HPD character state transitions
  # for each permissible genotype substitution
  ## Define all possible unphased genotype states
  states <- c("AA", "CC", "GG", "TT", "AC",
            "AG", "AT", "CG", "CT", "GT")

  ## Create empty dataframe to store the count of
  ## 95% HPD character state transitions
  hpd.df <- expand.grid(states,states) %>%
              filter(Var1 != Var2) %>%
              apply(., 1, function(x){
                from <- strsplit(x[1], "")[[1]]
                to <- strsplit(x[2], "")[[1]]
                if (from[1] == from[2] & from[1] %in% to) {
                  return(x)
                }
                if (sum(from %in% to) == 1) {
                  return(x)
                }
              }) %>%
              compact() %>%
              bind_rows() %>%
              arrange(Var1) %>% 
              rename(from = Var1, to = Var2) %>%
              mutate(lower_95hpd = 0,
                    upper_95hpd = 0)

  ## Create density object for SCM posteriors                  
  dens.data <- density(multiSimmap)

  ## Obtain the 95% HPD character state transitions and 
  ## populate hpd.df
  hpd.df <- lapply(1:nrow(hpd.df), function(x) {
              # Get row from hpd.df
              row <- hpd.df[x,]
              # Create string for genotype transition
              trans.str <- paste0(row[1, 1], "->", row[1, 2])

              # If transition is present in multiSimmap density data
              # populate hpd.df with 95% HPD transition counts
              if (trans.str %in% dens.data$trans) {
                i <- which(dens.data$trans == trans.str)
                row$lower_95hpd <- dens.data$hpd[[i]][1,1]
                row$upper_95hpd <- dens.data$hpd[[i]][1,2]
                return(row)
              } else {
                return(row)
              }
            }) %>%
            bind_rows()

  # Plotting
  if (plot == TRUE) {
    # Set colors for all 10 SNP genotype state
    cols <- setNames(brewer.pal(n = 10, name = "Set3"),
                  c("AA", "CC", "GG", "TT", "AC",
                    "AG", "AT", "CG", "CT", "GT"))
    # Use `phytools` default plotting method for multiSimmap
    plot(summary(multiSimmap),
         type = "phylogram",
         direction = "downwards",
         colors = cols,
         fsize = 1,
         ftype = "off",
         lwd = 1,
         cex = c(0.25, 0.25),
         cex = c(0.5, 0.5),
         mar = c(2, 0.1, 2, 0.1),
         outline = FALSE)
    
    # Add title (if specified)
    title(main = title)

    # Add legend (if specified)
    if (legend == TRUE) {
      legend("bottom",
             horiz = TRUE,
             xpd = TRUE,
             legend = c("AA", "CC", "GG", "TT", "AC",
                        "AG", "AT", "CG", "CT", "GT"),
             xjust = 0.5,
             yjust = 0.5,
             pch = 22,
             pt.cex = 2,
             pt.bg = cols,
             inset = c(0, -0.02),
             bty = "n",
             cex = 1)
    }

    # Add scale bar
    add.scale.bar(length = 0.01)

    # If state change detected, plot on tree as black diamond
    if (sum(detect_state_changes(tree$edge, genotype_states)) > 0) {
      edgelabels(text = "",
                 edge = which(detect_state_changes(tree$edge, genotype_states) == 1),
                 col = "black",
                 frame = "none",
                 pch = 18,
                 cex = 1.5,
                 adj = c(0.5, 0.5))
    }
  }

  # Return a list containing:
  #   1. scm_summary:
  #       Dataframe of genotype state posteriors
  #   2. assigned:
  #       Binary vector of whether a character state change occurred
  #       and reported in same order as <phylo>$edge.length
  #   3. constate_PPs:
  #       Posterior estimate for each node/tip (in node order)
  #   4. QloGL:
  #       Log likelihood of Q matrix given observations of tree and 
  #       tip states.
  #       ! Note that QlogL will have to be revised if Q is input
  #       ! as a prior distribution rather than fixed.
  #   5. Number of SCM replicates
  #   6. Count of 95% HPD character state transitions
  return(list("gt_posteriors" = out,
              "assigned_edges" = detect_state_changes(tree$edge, genotype_states),
              "consensus_posteriors" = constate_posteriors,
              "QlogL" = multiSimmap[[1]]$logL[1],
              "scm_reps" = length(multiSimmap),
              "hpd_counts" = hpd.df))
}

#* Output a summary for indel SCM, optionally plot
#! See summarise_scm.snp for notes on function elements; this function
#! will be removed in future updates when consolidated with summarise_scm.snp.
summarise_scm.indel <- function(multiSimmap, PPthreshold=0.95, plot = FALSE,
                                legend = TRUE, title = NA) {
  #load primary df from input
  scm_summary <- summary(multiSimmap)
  scm_summary <- scm_summary$ace %>% as.data.frame()
  tree <- multiSimmap[[1]]
  
  #annotate posteriors with consensus state and node numbers (renumbering tip names to node number)
  scm_summary$constate <- apply(scm_summary, 1, pull_consensus_state, threshold = PPthreshold)
  scm_summary$node <- rownames(scm_summary)
  scm_summary[(tree$Nnode+1):nrow(scm_summary), ]$node <- 1:(length(tree$tip.label))
  
  #pull vector of node consensus states
  genotype_states <- arrange(scm_summary, as.integer(node)) %>% 
    pull(constate)
  
  #format state posteriors df output 
  out <- scm_summary %>% 
    arrange(as.integer(node)) %>% 
    select(!c(node, constate))
  rownames(out) <- NULL
  cols_to_add <- setdiff(c("REF", "HET", "ALT"), colnames(out))
  out[,cols_to_add] = 0
  out <- out %>% 
    select("REF", "HET", "ALT")
  
  #pull consensus state posteriors (in same order as out)
  constate_posteriors <- apply(out, 1, max)
  
  #plotting
  if ( plot == TRUE ){
    cols <- setNames(c("gray", "royalblue", "tomato"),
                   c("REF", "HET", "ALT"))
    plot(summary(multiSimmap),
         type = "phylogram",
         direction = "downwards",
         colors = cols,
         fsize = 1,
         ftype = "off",
         lwd = 1,
         #offset=0.4,
         #ylim=c(-1,Ntip(tree)),
         cex = c(0.25, 0.25),
         cex = c(0.5, 0.5),
         mar = c(2, 0.1, 2, 0.1),
         outline = FALSE)
    title(main = title)
    if ( legend == TRUE ){
      legend("bottom",
             horiz = TRUE,
             xpd = TRUE,
             legend = c("0/0", "0/1", "1/1"),
             xjust = 0.5,
             yjust = 0.5,
             pch = 22,
             pt.cex = 2,
             pt.bg = cols,
             inset = c(0, -0.02),
             #title.adj = 0,
             bty = "n",
             cex = 1)
    }
    #add.scale.bar(length = 0.01, y = 0.1, x = 0)
    add.scale.bar(length = 0.01)
    if ( sum(detect_state_changes(tree$edge, genotype_states)) > 0) {
      edgelabels(text = "",
                 edge = which(detect_state_changes(tree$edge, genotype_states) == 1),
                 #col = "white",
                 #bg = "black",
                 col = "black",
                 frame = "none",
                 pch = 18,
                 cex = 1.5,
                 #frame = "rect",
                 adj = c(0.5, 0.5))
    }
  }
  
  #output
  return(list("scm_summary" = out,
              "assigned" = detect_state_changes(tree$edge, genotype_states),
              "constate_PPs" = constate_posteriors,
              "QlogL" = multiSimmap[[1]]$logL[1]))
}

#* Return a vector of consensus states (>= threshold) all rows in summary(multiSimmap)$ace
pull_consensus_state  <-  function(row, threshold) {
  # ARGUMENTS:
  #   - row:
  #       row (df) from summary(multiSimmap)$ace
  #   - threshold:
  #       lower-bound posterior cutoff for confident genotype state call
  
  if (max(row) >= threshold) {
    names(row)[which.max(row)]
  } else {
    # Return NA if all state posteriors are < threshold
    NA
  }
} 

#* Returns a vector (in same order as <phylo>$edge) 
#* where state changes have occurred
detect_state_changes <- function(tree_edges, states){
  # ARGUMENTS:
  #   - tree_edges:
  #       matrix from <phylo>$edge, where:
  #         1. tree_edges[,1] is the parent node index
  #         2. tree_edges[,2] is the daughter node index
  #         3. tree_edges[n,] are edge indices
  #   - states:
  #       Consensus state in node order

  # If root state is NA, return 0's for all branches and exit
  if ( is.na(states[tree_edges[1, ][1]]) ){
    return(rep(0, length(states)))
  }

  # Pull root genotype from first index of states
  root_state <- states[1]
  
  # For each edge of tree (row in tree_edges):
  # vector[1] = parent node
  # vector[2] = daughter node
  apply(tree_edges, 1, function(x){

    # If parent state uncertain, assign parent state to most recent high confidence state;
    # this will always resolve to a non-NA state given the conditional test at the root.
    while ( is.na(states[x[1]]) ) {
      x[1]  <-  tree_edges[tree_edges[, 2] == x[1]][1]
    }
    
    # If daughter state is uncertain, no state change on edge
    if ( is.na(states[x[2]]) ) {
      return(0)
    } 
    
    # If parent state == daughter state, no state change on edge
    if ( states[x[1]] == states[x[2]] ){
      return(0)
    } 
    
    # Else, state change present
    else {
      return(1)
    }
  })
} 

#* Plot scm density histogram for INDELS. Reproduces the 
#* output from plot(density(multiSimmap)) using ggplot.
plot_density_hist.indel <- function(multiSimmap_density){
  # ARGUMENTS:
  #   - multiSimmap_density:
  #       Output from running density() on multiSimmap
  #       object. If multiSimmap object input, function
  #       will run density(multiSimmap). 

  # If input hasn't been analyzed with density(), do so
  if (any(class(multiSimmap_density) == "multiSimmap")) {
    multiSimmap_density <- density(scm_example)
  }
  
  # Create vector of permitted indel state changes
  trans <- c("REF->HET", "ALT->HET",
           "HET->REF", "HET->ALT")
  
  # Density data
  df.dens <- lapply(trans, function(x) {
    if (! x %in% multiSimmap_density$trans) {
      # If state change not measured in multiSimmap_density, 
      # store empty observations to df
      df <- data.frame(state_trans = x,
                     num_trans = 0,
                     count = NA,
                     density = 1.0)
    } else {
      # If state change measured in multiSimmap_density,
      # store observations in df 
      df <- data.frame(state_trans = x,
                     num_trans = multiSimmap_density$p[[x]]$mids,
                     count = multiSimmap_density$p[[x]]$counts,
                     density = multiSimmap_density$p[[x]]$density)
    }
    return(df)
  }) %>% bind_rows()

  # Convert state_trans field to factor
  df.dens$state_trans <- factor(df.dens$state_trans, levels = trans)

  # Cateogrize allele state and the zygosity gain/loss
  df.dens <- df.dens %>% 
    mutate(allele_cat = case_when(state_trans %in% c("REF->HET", "HET->REF") ~ as.factor("Homozygous ref."),
                                state_trans %in% c("ALT->HET", "HET->ALT") ~ as.factor("Homozygous alt.")),
           heterozygosity = case_when(state_trans %in% c("REF->HET", "ALT->HET") ~ as.factor("gain"),
                                    state_trans %in% c("HET->REF", "HET->ALT") ~ as.factor("loss")))
  
  # Extract high probability density (HPD) intervals
  df.hpd <- lapply(trans, function(x) {
    if (! x %in% multiSimmap_density$trans) {
      df <- data.frame(state_trans = x,
                     low = 0,
                     high = 0)
    } else {
      df <- data.frame(state_trans = x,
                     low = multiSimmap_density$hpd[[x]][1],
                     high = multiSimmap_density$hpd[[x]][2])
    }
    return(df)
  }) %>% bind_rows()

  # Convert HPD states to factors
  df.hpd$state_trans <- factor(df.hpd$state_trans, levels = trans)
  
  # Cateogrize allele state and the zygosity gain/loss
  df.hpd <- df.hpd %>% 
    mutate(allele_cat = case_when(state_trans %in% c("REF->HET", "HET->REF") ~ as.factor("Homozygous ref."),
                                state_trans %in% c("ALT->HET", "HET->ALT") ~ as.factor("Homozygous alt.")),
           heterozygosity = case_when(state_trans %in% c("REF->HET", "ALT->HET") ~ as.factor("gain"),
                                    state_trans %in% c("HET->REF", "HET->ALT") ~ as.factor("loss")),
           med=((low+high)/2))
  
  # Join HPD intervals to density data
  df.hpd <- df.dens %>% 
    group_by(state_trans) %>% 
    select(state_trans, density) %>% 
    summarise(max_den = max(density)) %>%
    mutate(max_den = max_den + 0.05) %>%
    right_join(df.hpd, by = "state_trans")
  
  # Structure additional annotations for plotting
  annotations <- data.frame(lab = c("gain of het.", "loss of het."),
                            allele_cat = factor("Homozygous alt.",
                            levels = c("Homozygous ref.", "Homozygous alt.")),
                            count = max(df.dens$num_trans),
                            density = c(0.5, -0.5))
  annotations2 <- data.frame(lab = c("Homozygous ref.", "Homozygous alt."),
                            allele_cat = factor(c("Homozygous ref.", "Homozygous alt."),
                            levels = c("Homozygous ref.", "Homozygous alt.")),
                            y = 1.2,
                            x = df.dens %>% pull(num_trans) %>% max() / 2)
  
  #Plot
  par(mfrow=c(1,1), ask = FALSE)
  g <- ggplot(data = df.dens %>% filter(heterozygosity == "gain"),
              mapping = aes(x = num_trans, y = density))+
        geom_col(mapping = aes(fill = state_trans),
                 position = "identity",
                 width = 0.95)+
        geom_col(data = df.dens %>% filter(heterozygosity == "loss"),
                 mapping = aes(x = num_trans, y = -1 * density, fill = state_trans),
                 position = "identity",
                 width = 0.95)+
        geom_hline(yintercept = 0, linewidth = 0.5, color = "white")+
        geom_text(data = annotations,
                  mapping = aes(label = lab, x = count, y = density, color = lab),
                  fontface = "bold",
                  angle = -90,
                  hjust = 0.5,
                  nudge_x = 1,
                  color = c("gray70", "gray30"),
                  size = 5)+
        geom_text(data = annotations2,
                  mapping = aes(x = x, y = y, label = lab),
                  hjust = 0.5,
                  vjust = 0.5,
                  fontface = "bold",
                  color = c("royalblue", "salmon"),
                  size = 5)+
        geom_segment(data = df.hpd %>% filter(heterozygosity == "gain"),
                     mapping = aes(y = max_den, yend = max_den, 
                                   x = low - 0.5, xend = high + 0.5),
                     color = "black")+
        geom_segment(data = df.hpd %>% filter(heterozygosity == "loss"),
                     mapping = aes(y = -max_den, yend = -max_den,
                                   x = low - 0.5, xend = high + 0.5),
                     color = "black")+
        geom_text(data = df.hpd %>% filter(heterozygosity == "gain"),
                  mapping = aes(y = max_den + 0.075, x = med),
                  color = "black",
                  hjust = 0.5,
                  label = "95% HPD")+
        geom_text(data = df.hpd %>% filter(heterozygosity == "loss"),
                  mapping = aes(y = -1 * (max_den + 0.075), x = med),
                  color = "black",
                  hjust = 0.5,
                  label = "95% HPD")+
        facet_wrap(~allele_cat, ncol=2) +
        scale_x_continuous(limits = c(min(df.dens$num_trans), max(df.dens$num_trans) + 1),
                           breaks = seq(0, max(df.dens$num_trans), 2)) +
        scale_y_continuous(limits = c(-1.2, 1.2), breaks=seq(-1, 1, 0.2),
                           labels = c(rev(seq(0, 1, 0.2)), seq(0.2, 1, 0.2))) +
        scale_fill_manual(values = c("royalblue", "salmon",
                                   "royalblue4", "salmon4"),
                          guide="none")+
        theme_cowplot() +
        theme(strip.text = element_blank(),
              strip.background = element_blank(),
              plot.title = element_text(face = "bold", hjust = 0.5),
              plot.subtitle = element_text(face = "italic", hjust = 0.5),
              axis.line.y = element_line()) +
        labs(x="mutations (count)", y="density")
  suppressWarnings(print(g))

  # Return list of dataframes
  return(list("densities" = df.dens %>% select(-allele_cat, -heterozygosity),
              "hpd95" = df.hpd %>% select(-allele_cat, -heterozygosity, -med, -max_den)))
}