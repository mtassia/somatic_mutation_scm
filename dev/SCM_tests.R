source("scripts/SCM_functions.R")

##### Data import #####
#Read tree
tr<-read.tree(file="trees/DH.phylogenetic_sites.synthetic_reference.vcf.cellphy.raxml.supportFBP") %>%
  root(outgroup="DH",resolve.root = T) %>% 
  drop.tip("DH") %>%
  ladderize()
#tree<-collapse_poor_supported_edges(tree,90)

#Read vcf
vcf<-read.vcfR("vcfs/DH.filtered.06282024.curated.funcotated.vcf.gz")
gt_list.snp<-compile_gt_states.snp(vcf)
gt_list.indel<-compile_gt_states.indel(vcf)

#Prep Q matrices
snp.q<-read_cellphy_model(bestModel_path = "models/DH.phylogenetic_sites.synthetic_reference.vcf.cellphy.raxml.bestModel")
indel.q<-generate_indel_model(indel_state_list = gt_list.indel,tree=tr)

##### SCM SNPS #####
scm_example<-runSCM_single(#x = sample(names(gt_list.snp),size = 1),
                           x = "chr5_1295113",
                           tree = tr,
                           gt_state_list = gt_list.snp,
                           Qmat = snp.q,
                           reduced = F,
                           reps = 1000,
                           cores=6)

par(mfrow=c(5,5))
null<-sapply(scm_example[sample(seq(1,1000),size = 25,replace = F)],
             plot,
             ftype="off",
             lwd=1.5)
par(mfrow=c(1,1),ask=F)
summarise_scm.snp(scm_example,plot=T)

par(mfrow=c(4,5))
for ( i in sample(names(gt_list.snp),size = 20,replace = F) ){
  scm<-runSCM_single(#x="chr7_124853116",
    #x=sample(names(gt_list),size = 1,replace = F),
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
scm_example<-runSCM_single(x = site,
  #x="chr3_195717011",
  #x="chr4_10145447",
  tree = tr,
  gt_state_list = gt_list.indel,
  Qmat = indel.q,
  reduced = F,
  reps = 1000,
  cores=6)
plot_density_hist.indel(scm_example)
summarise_scm.indel(scm_example,plot=T,title=site)

par(mfrow=c(5,5))
null<-sapply(scm_example[sample(seq(1,1000),size = 25,replace = F)],
             plot,
             ftype="off",
             lwd=1.5,
             colors=c("REF"="gray","HET"="royalblue","ALT"="tomato"))

par(mfrow=c(1,1),ask=F)
summarise_scm.indel(scm_example,plot=T)

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

scm_summary<-summary(scm_example)
Prior_post<-data.frame(gt_list.indel$chr3_195717011) %>%
  select(Indiv,REF,HET,ALT) %>%
  pivot_longer(cols=2:4,names_to="genotype",values_to="prior") %>%
  left_join(as.data.frame(scm_summary$tips) %>%
              rownames_to_column(var="Indiv") %>%
              pivot_longer(cols=2:4,names_to="genotype",values_to="posterior"))
Prior_post$genotype<-factor(Prior_post$genotype,levels=c("REF","HET","ALT"))

ggplot(Prior_post,aes(x=prior,y=posterior,fill=genotype)) +
  geom_point(size=4,shape=21) +
  geom_abline(slope=1,linetype=2,color="gray") +
  facet_wrap(~genotype) +
  scale_fill_manual(values=c("REF"="gray","HET"="royalblue","ALT"="tomato"),guide="none")+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1)) +
  theme_cowplot() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face="bold"),
        plot.title=element_text(hjust=0.5,face="bold"))+
  labs(title="INDEL - chr3:195717011")

#### TODO ####
#Write method to compose results into a summary data frame and the posterior estimates into an HDF5-compatible array.
#Summary output should have the following data:
#c("locus","root","derived","n_state_changes","class","logL_Q","edge_profile")
