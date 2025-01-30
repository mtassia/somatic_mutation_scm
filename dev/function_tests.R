source("dev/bin/SCM_functions.R")
options(vsc.dev.args = list(width = 1200,height=1000))

##### Data import #####
#Read tree
tr<-read.tree(file="dev/data/test.vcf.cellphy.raxml.supportFBP") %>%
  root(outgroup="outgroup", resolve.root = T) %>%
  drop.tip("outgroup") %>%
  ladderize()
#tr<-collapse_poor_supported_edges(tr,90)

#Read vcf
vcf<-read.vcfR("dev/data/test.vcf.gz")
gt_list.snp<-compile_gt_states.snp(vcf)
# gt_list.indel<-compile_gt_states.indel(vcf)

#Prep Q matrices
snp.q<-read_cellphy_model(bestModel_path = "dev/data/test.vcf.cellphy.raxml.bestModel")
# indel.q<-generate_indel_model(indel_state_list = gt_list.indel, tree = tr)

##### SCM SNPS #####
#Focal variant
scm_example<-runSCM_single(x = "chr5_1295113",
                           tree = tr,
                           gt_state_list = gt_list.snp,
                           Qmat = snp.q,
                           reduced = FALSE,
                           reps = 1000,
                           cores = 6)
par(mfrow=c(1,1),ask=F)
summarise_scm.snp(multiSimmap = scm_example,
                  locus = "chr5:1295113",
                  plot = FALSE)

#Visualize individual SCM generations
par(mfrow=c(5,5))
null<-sapply(scm_example[sample(seq(1,1000),size = 25,replace = F)],
             plot,
             ftype="off",
             lwd=1.5)

#Visualize random sample of SCM summaries
par(mfrow=c(4,5))
for ( i in sample(names(gt_list.snp),size = 20,replace = F) ) {
  scm<-runSCM_single(
    x=i,
    tree=tr,
    gt_state_list = gt_list.snp,
    Qmat = snp.q,
    reps = 100,
    cores = 6,
    reduced = T)
  
  print(paste0("Q logL: ",scm[[1]]$logL[1]))
  summarise_scm.snp(scm,plot = T,title=i,legend = F)
  lapply(scm,function(x) {x$logL[1]}) %>% unlist() %>% summary()
  tmp<-summary(scm)
}

##### SCM INDELS #####
site=sample(names(gt_list.indel),size=1)
scm_example <- runSCM_single(
  x = site,
  tree = tr,
  gt_state_list = gt_list.indel,
  Qmat = indel.q,
  reduced = F,
  reps = 1000,
  cores=6)
par(mfrow=c(1,1),ask=F)
plot_density_hist.indel(scm_example)
summarise_scm.indel(scm_example,plot=T,title=site)

par(mfrow=c(5,5))
null<-sapply(scm_example[sample(seq(1,1000),size = 25,replace = F)],
             plot,
             ftype="off",
             lwd=1.5,
             colors=c("REF"="gray","HET"="royalblue","ALT"="tomato"))

par(mfrow=c(4,5))
for ( i in sample(names(gt_list.indel),size = 20,replace = F) ){
  scm<-runSCM_single(x=i,
    tree=tr,
    gt_state_list = gt_list.indel,
    Qmat = indel.q,
    reps = 1000,
    cores = 6,
    reduced = T)
  print(paste0("Q logL: ",scm[[1]]$logL[1]))
  summarise_scm.indel(scm,plot = T,title=i,legend = F)
  lapply(scm,function(x) {x$logL[1]}) %>% unlist() %>% summary()
  tmp<-summary(scm)
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
# Create h5 file with groups
h5createFile("dev/data/test.h5")
h5createGroup("dev/data/test.h5", "gt_posteriors")
h5createGroup("dev/data/test.h5", "assigned_edges")
h5createGroup("dev/data/test.h5", "consensus_posteriors")
h5createGroup("dev/data/test.h5", "QlogL")
h5createGroup("dev/data/test.h5", "scm_reps")
h5createGroup("dev/data/test.h5", "hpd_counts")
h5ls("dev/data/test.h5")

loc.len <- length(grep(pattern = "chr9",
                  x = names(gt_list.snp),
                  value = TRUE))
i <- 1
for (loc in grep(pattern = "chr9",
                 x = names(gt_list.snp),
                 value = TRUE)) {
  print(paste0("Processing locus ", i, " of ", loc.len))
  start_time <- Sys.time()

  scm <- runSCM_single(x = loc,
                       tree = tr,
                       gt_state_list = gt_list.snp,
                       Qmat = snp.q,
                       reduced = FALSE,
                       reps = 100,
                       cores = 6)
  scm <- summarise_scm.snp(multiSimmap = scm,
                           locus = loc,
                           plot = FALSE)

  h5write(obj = scm$gt_posteriors,
          file = "dev/data/test.h5",
          name = paste0("gt_posteriors/", loc))
  h5write(obj = scm$assigned_edges,
          file = "dev/data/test.h5",
          name = paste0("assigned_edges/", loc))
  h5write(obj = scm$consensus_posteriors,
          file = "dev/data/test.h5",
          name = paste0("consensus_posteriors/", loc))
  h5write(obj = scm$QlogL,
          file = "dev/data/test.h5",
          name = paste0("QlogL/", loc))
  h5write(obj = scm$scm_reps,
          file = "dev/data/test.h5",
          name = paste0("scm_reps/", loc))
  h5write(obj = scm$hpd_counts,
          file = "dev/data/test.h5",
          name = paste0("hpd_counts/", loc))
  
  end_time <- Sys.time()
  print(paste0("Time elapsed: ", end_time - start_time))
  i <- i + 1
}
h5ls("dev/data/test.h5")

# Read h5 file into object
h5f <- H5Fopen("dev/data/test.h5")

h5closeAll()
