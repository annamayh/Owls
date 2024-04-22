## colour evoltion
library(tidyverse)
library(SNPRelate)
library(genedroppeR)
library(sequoia)

setwd("/Users/ahewett1/Documents")


# load the data with many phenotypes, not all are useful
NewPhenoColor_Order3K = read.table("Gene_drop/NewPhenoColor_Order3K.txt", header = TRUE)
# load the table with the genotypes, you will meed the package "SNPRelate"
gen_mat <- snpgdsGetGeno('Gene_drop/All3102_NewNames_RP502SNPs.variantsList.gds')
# change the genotypes so 0 mean homozygote white
gen_mat[,1] = abs(gen_mat[,1]-2)
gen_mat[,2] = abs(gen_mat[,2]-2)
# make females haploid
gen_mat[NewPhenoColor_Order3K$Sex==2 & gen_mat[,3]==2 & (! is.na(NewPhenoColor_Order3K$Sex)),3] = 1

# format the table, I now realise I could have just send this one ... 
TabSummary = data.frame(NewPhenoColor_Order3K$ID, NewPhenoColor_Order3K$MeanColor, NewPhenoColor_Order3K$SurfaceBlackSpotsCm, NewPhenoColor_Order3K$SexWGS, NewPhenoColor_Order3K$StageBin, gen_mat[,1], gen_mat[,2], gen_mat[,3])
colnames(TabSummary) = c("ID", "MeanColor", "SurfaceBlackSpotsCm", "SexWGS", "StageBin", "GenoCorin", "GenoMC1R", "GenoZ")
TabSummary = TabSummary[! is.na(TabSummary$ID),]
rownames(TabSummary) = TabSummary$ID

# Now you can play with the table :)


genotypes=TabSummary%>%
  select(ID, GenoCorin, GenoMC1R, GenoZ)%>%
  rename(RingId=ID)


## pedigree 

ped_birth_yr2=ped_birth_yr%>%
  unique()%>%
  mutate(year=as.numeric(year))%>%
  filter(year>=1995)%>%
  filter(year<2020)%>%
  mutate(cohort=(year-min(year))+1)%>%
  filter(RingId!="889913") ##sort this out later

cohort_yr=ped_birth_yr2%>%
  select(year, cohort)%>%
  unique()


owl_col_gens=genotypes%>%
  right_join(ped_birth_yr2)


owl.summ <- summary_cohort(id = owl_col_gens$RingId,
                               mother = owl_col_gens$mumID,
                               father = owl_col_gens$dadID,
                               cohort = owl_col_gens$cohort,
                               genotype = owl_col_gens$GenoMC1R)

owl.summ


owl.summ%>%left_join(cohort_yr)%>%
ggplot(aes(year, B)) +
  geom_line() +
  stat_smooth(method = "lm") +
  ggtitle("Temporal dynamics of red allele")

owl.summ%>%left_join(cohort_yr)%>%
ggplot(aes(year, PropGenotyped)) +
  geom_line()  +
  ggtitle("Proportion of Genotyped IDs per year")



owl.mc1r <- genedrop_snp(owl_col_gens$RingId,
                         mother = owl_col_gens$mumID,
                         father = owl_col_gens$dadID,
                         cohort = owl_col_gens$cohort,
                         genotype = owl_col_gens$GenoMC1R
                          nsim = 1000,
                           n_founder_cohorts = 5,
                           fix_founders = F,
                           verbose = T,
                           interval = 200)
