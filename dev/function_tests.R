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
for ( i in sample(names(gt_list.snp), size = 20, replace = FALSE ) ) {
  scm<-runSCM_single(
    x=i,
    tree=tr,
    gt_state_list = gt_list.snp,
    Qmat = snp.q,
    reps = 100,
    cores = 6,
    reduced = T)
  
  print(paste0("Q logL: ",scm[[1]]$logL[1]))
  summarise_scm.snp(scm,
                    plot = T,
                    title = i,
                    legend = FALSE,
                    locus = i)
  lapply(scm,function(x) {x$logL[1]}) %>% unlist() %>% summary()
  tmp <- summary(scm)
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
# TODO: Add ML tree as static group
# TODO: Add metadata as static group (chr, records, runtime, date, etc.)
# TODO: Add scm_reps as static group
# TODO: Add Q matrix as static group

source("dev/bin/SCM_functions.R")
df <- multi_scm(gt_state_list = gt_list.snp,
          tree = tr,
          Q = snp.q,
          scm_its = 100,
          cores = 6,
          h5f_path = "dev/data/test.h5",
          chr = "chr7",
          overwrite = TRUE)
