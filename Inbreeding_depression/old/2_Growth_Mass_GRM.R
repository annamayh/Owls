library(tidyverse)
library(brms)
library(corpcor)

##################################################################################
########################### ~~ MASS ~~ #############################################
#####################################################################################

mass_df=read.table("./input_dfs/mass_all_pheno_df.txt",sep=",", header=T)

mass_df$clutch_merge=as.factor(mass_df$clutch_merge)
mass_df$GeneticSex=as.factor(mass_df$GeneticSex)
mass_df$RingId=as.factor(mass_df$RingId)
mass_df$year=as.factor(mass_df$year)
mass_df$Observer=as.factor(mass_df$Observer)
mass_df$nestboxID=as.factor(mass_df$nestboxID)
mass_df$rank=as.numeric(mass_df$rank)

## read in GRM
grm=read_rds("./input_dfs/All3085_AUTOSAUMES_RP502SNPs.RDS")
grm_filt=grm[rownames(grm)%in%mass_df[["RingId"]],colnames(grm)%in%mass_df[["RingId"]]]##filtering grm for ids only in df

mass_df[,'RingId_pe']=mass_df[,'RingId'] ##add permanent env variable to get h2 estimate

grm_filt_pd <- make.positive.definite(grm_filt)
GRM <- as(grm_filt_pd, "dgCMatrix")

##priors
prior_mass<- c(
  prior(normal(3.5, 2), nlpar = "asym",  class="b"),##
  prior(normal(4.5, 1), nlpar = "b",  class="b"), ## 
  prior(normal(1, 0.5), nlpar = "c",  class="b"), ## 
  
  prior(student_t(3, 0, 2.5), class = "sigma", lb=0),
  prior(student_t(3, 0, 2.5), class="sd",nlpar = "asym", lb=0), # 
  prior(student_t(3, 0, 2.5),  class="sd", nlpar = "b", lb=0),
  prior(student_t(3, 0, 2.5),  class="sd", nlpar = "c", lb=0), 
  
  prior(cauchy(0, 1), class="sd", group="RingId_pe", nlpar = "asym", lb=0)#, #
  # prior(cauchy(0, 0.5),  class="sd", group="RingId_pe", nlpar = "b", lb=0),
  # prior(cauchy(0, 0.2),  class="sd", group="RingId_pe", nlpar = "c", lb=0)
  
)

#for now, only fit GRM for asymtote 
growth_mass.mod.GRM=brm(
  ## model 
  bf(mass_scale ~ asym * exp(-b*(c)^age_days),
     asym ~ 1 + FuniWE + rank + GeneticSex + (1|clutch_merge) + (1|Observer) + (1|year) + (1|nestboxID) + (1|gr(RingId, cov=GRM))+(1|RingId_pe),
     b ~ 1 + FuniWE + rank + GeneticSex + (1|clutch_merge) + (1|year) +(1|RingId_pe),
     c ~ 1 + FuniWE + rank + GeneticSex +  (1|year) +(1|RingId_pe), 
     nl=TRUE),
  data=mass_df, 
  data2 = list(GRM = GRM),
  family = gaussian(),
  chains = 4,
  prior = prior_mass,
  control = list(adapt_delta = 0.85),
  init = 0, 
  cores = 4,
  iter = 60000, 
  warmup = 10000
  
)

summary(growth_mass.mod.GRM)


saveRDS(growth_mass.mod.GRM,file="./outputs/growth_mass_GRM_Funi.RDS")
