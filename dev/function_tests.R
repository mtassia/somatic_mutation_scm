source("dev/bin/SCM_functions.R")
options(vsc.dev.args = list(width = 1200,height=1000))
par(mar = c(2,2,2,2))

##### Data import #####
#Read tree
tr<-read.tree(file=paste0(getwd(), "/dev/data/test.vcf.cellphy.raxml.supportFBP")) %>%
  root(outgroup="outgroup", resolve.root = T) %>%
  drop.tip("outgroup") %>%
  ladderize()
#tr<-collapse_poor_supported_edges(tr,90)

#Read vcf
vcf<-read.vcfR(paste0(getwd(), "/dev/data/test.vcf.gz"))
gt_list.snp<-compile_gt_states.snp(vcf)
# gt_list.indel<-compile_gt_states.indel(vcf)

#Prep Q matrices
snp.q <- read_cellphy_model(bestModel_path = paste0(getwd(), "/dev/data/test.vcf.cellphy.raxml.bestModel"))
# indel.q<-generate_indel_model(indel_state_list = gt_list.indel, tree = tr)

##### SCM SNPS #####
#Focal variant
# var <- "chr5_1295113"
# var <- "chr22_50057483"
# var <- "chr9_5073770"
var <- sample(names(gt_list.snp), size = 1)
scm_example <- runSCM_single(x = var,
                             tree = tr,
                             gt_state_list = gt_list.snp,
                             Qmat = snp.q,
                             reduced = FALSE,
                             reps = 1000,
                             cores = 6)
par(mfrow = c(1, 1), ask = FALSE)
summarise_scm.snp(multiSimmap = scm_example,
                  locus = var,
                  plot = TRUE)

par(mar = c(5,5,5,5))
plot(density(scm_example))

#Plot state priors and posteriors w.r.t. tree
plot_state_probs_on_tree(locus_str = var,
                         gt_state_list = gt_list.snp,
                         scm_summary = summarise_scm.snp(multiSimmap = scm_example,
                                                         locus = var,
                                                         plot = FALSE),
                         tree = tr)

#Visualize individual SCM generations
par(mfrow = c(5, 5))
null <- sapply(scm_example[sample(seq(1, 1000), size = 25, replace = FALSE)],
             plot,
             ftype = "off",
             lwd = 1.5)

#Visualize random sample of SCM summaries
par(mfrow = c(4, 5))
for (i in sample(names(gt_list.snp), size = 20, replace = FALSE)) {
  scm <- runSCM_single(
    x = i,
    tree = tr,
    gt_state_list = gt_list.snp,
    Qmat = snp.q,
    reps = 100,
    cores = 6,
    reduced = TRUE)
  
  print(paste0("Q logL: ", scm[[1]]$logL[1]))
  summarise_scm.snp(scm,
                    plot = TRUE,
                    title = i,
                    legend = FALSE,
                    locus = i)
  lapply(scm, function(x) {x$logL[1]}) %>% unlist() %>% summary()
  tmp <- summary(scm)
}

##### SCM INDELS #####
site = sample(names(gt_list.indel),size=1)
scm_example <- runSCM_single(
  x = site,
  tree = tr,
  gt_state_list = gt_list.indel,
  Qmat = indel.q,
  reduced = FALSE,
  reps = 1000,
  cores = 6)
par(mfrow = c(1, 1), ask = FALSE)
plot_density_hist.indel(scm_example)
summarise_scm.indel(scm_example, plot = TRUE, title = site)

par(mfrow = c(5, 5))
null <- sapply(scm_example[sample(seq(1, 1000), size = 25, replace = FALSE)],
            plot,
            ftype = "off",
            lwd = 1.5,
            colors = c("REF" = "gray", "HET" = "royalblue", "ALT" = "tomato"))

par(mfrow = c(4, 5))
for (i in sample(names(gt_list.indel), size = 20, replace = FALSE)){
  scm < -runSCM_single(x = i,
    tree = tr,
    gt_state_list = gt_list.indel,
    Qmat = indel.q,
    reps = 1000,
    cores = 6,
    reduced = TRUE)
  print(paste0("Q logL: ",scm[[1]]$logL[1]))
  summarise_scm.indel(scm, plot = TRUE, title = i, legend = FALSE)
  lapply(scm, function(x) {x$logL[1]}) %>% unlist() %>% summary()
  tmp <- summary(scm)
}

scm_summary <- summary(scm_example)
Prior_post <- data.frame(gt_list.indel[[site]]) %>%
  select(Indiv, REF, HET, ALT) %>%
  pivot_longer(cols=2:4, names_to = "genotype", values_to = "prior") %>%
  left_join(as.data.frame(scm_summary$tips) %>%
              rownames_to_column(var = "Indiv") %>%
              pivot_longer(cols = 2:ncol(.),
                           names_to = "genotype",
                           values_to = "posterior")) %>%
  mutate(posterior = case_when(is.na(posterior) ~ 0,
                               TRUE ~ posterior))
Prior_post$genotype <- factor(Prior_post$genotype,
                              levels = c("REF", "HET", "ALT"))

ggplot(Prior_post, aes(x = prior, y = posterior, fill = genotype)) +
  geom_point(size = 4, shape = 21) +
  geom_abline(slope = 1, linetype = 2, color = "gray") +
  facet_wrap(~genotype) +
  scale_fill_manual(values = c("REF" = "gray",
                               "HET" = "royalblue",
                               "ALT" = "tomato"),
                    guide = "none")+
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_cowplot() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"))

##### HDF5 TESTS #####
source("dev/bin/SCM_functions.R")
# file.remove("dev/data/test.h5")
multi_scm(gt_state_list = gt_list.snp,
          tree = tr,
          Q = snp.q,
          scm_its = 100,
          cores = 6,
          h5f_path = "dev/data/test.h5",
          chr = "all",
          overwrite = FALSE,
          dryrun = TRUE)

h5 <- H5Fopen("dev/data/test_chr21.h5")
h5_summary <- bind_rows(h5$`summary`) %>%
              tibble()
events <- lapply(h5$`assigned_edges`, sum) %>%
  unlist() %>%
  stack()
h5closeAll()

ggplot(h5_summary %>% rownames_to_column(var = "index")) +
  geom_point(aes(y = runtime / iterations, x = index, color = QlogL), alpha = 0.5) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_color_viridis(option = "A", direction = -1) +
  theme_cowplot() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 24, face = "bold")) +
  labs(y = "Mean runtime per SCM (sec)", x = "Locus index")

ggplot(h5_summary %>%
         group_by(muts95low,muts95high) %>%
         summarise(count = n()),
       aes(x = as.factor(muts95low), y = as.factor(muts95high), fill = log10(count))) +
  geom_tile() +
  geom_label(aes(label = count), fill = "white", size = 4) +
  coord_equal() +
  scale_fill_viridis(option = "A") +
  theme_cowplot() +
  theme(panel.background = element_rect(fill = "black"),
        axis.title = element_text(size = 24, face = "bold")) +
  labs(x = "95% HPD lower bound",
       y = "95% HPD upper bound",
       fill = "log10(count)")
ggplot(events, aes(x = values)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0,6,1)) +
  theme_cowplot() +
  labs(x = "Mutations assigned to tree (PP ≥ 95%)",
       y = "Count")
h5_summary %>%
  mutate(HPDrange = (muts95high - muts95low) + 1) %>%
  group_by(HPDrange,muts_assigned) %>%
  summarise(count = n()) %>%
  ggplot(aes(y = muts_assigned, x = HPDrange, fill = log10(count))) +
    geom_tile() +
    geom_label(aes(label = count), fill = "white", size = 6) +
    coord_equal() +
    scale_fill_viridis(option = "A", begin = 0.25, guide = "none") +
    theme_cowplot() +
    theme(panel.background = element_rect(fill = "black"),
          axis.title = element_text(size = 24, face = "bold")) +
    labs(x = "Mutation events range (95% HPD, inclusive)",
         y = "Assigned mutations (PP ≥ 95%)")

ggplot(h5_summary, 
       aes(x = as.factor(muts_assigned),
           y = QlogL)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray80") +
  geom_quasirandom() +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  theme_cowplot() +
  labs(x = "Assigned mutations (PP ≥ 95%)",
       y = "Q matrix log-likelihood")

