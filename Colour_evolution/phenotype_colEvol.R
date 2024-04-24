## colour evoltion
library(tidyverse)
library(SNPRelate)
library(genedroppeR)
library(sequoia)

setwd("/Users/ahewett1/Documents")

mc1r=read.csv("Gene_drop/MC1R.csv", header = T)%>%
  rename(RingId=Ring, MC1R=V126I)

mc1r[mc1r=="VV"]=2
mc1r[mc1r=="VI"]=1
mc1r[mc1r=="II"]=0


## pedigree 

ped_birth_yr2=ped_birth_yr%>%
  unique()%>%
  mutate(year=as.numeric(year))%>%
   filter(year>=1997)%>%
   filter(year<2018)%>%
  mutate(cohort=(year-min(year))+1)%>%
  filter(RingId!="889913") ##sort this out later

cohort_yr=ped_birth_yr2%>%
  select(year, cohort)%>%
  unique()


owl_col_gens=mc1r%>%
  right_join(ped_birth_yr2)

owl_col_gens$MC1R=as.numeric(owl_col_gens$MC1R)

owl.summ <- summary_cohort(id = owl_col_gens$RingId,
                               mother = owl_col_gens$dadid,
                               father = owl_col_gens$momid,
                               cohort = owl_col_gens$cohort,
                               genotype = owl_col_gens$MC1R)

## Locus summary AA, AB, BB corresponds to 0 (II), 1, 2 (VV)

owl.summ


owl.summ%>%left_join(cohort_yr)%>%
ggplot(aes(year, A)) +
  geom_line() +
  stat_smooth(method = "lm") +
  ggtitle("Temporal dynamics of 'I' allele")

owl.summ%>%left_join(cohort_yr)%>%
ggplot(aes(cohort, PropGenotyped)) +
  geom_line()  +
  ggtitle("Proportion of Genotyped IDs per year")



owl.mc1r <- genedrop_snp(owl_col_gens$RingId,
                         mother = owl_col_gens$momid,
                         father = owl_col_gens$dadid,
                         cohort = owl_col_gens$cohort,
                         genotype = owl_col_gens$MC1R,
                          nsim = 1000,
                          n_founder_cohorts = 5,
                          fix_founders = F,
                          verbose = T,
                          interval = 200)


owl.mc1r.summ <- summary_genedrop(owl.mc1r)
plot_genedrop_results(owl.mc1r.summ)


plot_genedrop_cumulative_change(owl.mc1r.summ)
