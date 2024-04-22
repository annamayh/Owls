
## using pedigree RM to see effect of region etc.
## quicker than running GRM but will use this as proxy for running the big boy on the cluster later 

library(tidyverse)
library(brms)


setwd("/Users/ahewett1/Documents")

#######################################################################
##################### ~~ TARSUS ~~ #################################### 
#####################################################################

tarsus_df=read.table("Inbreeding_depression_owls/pheno_df/tarsus_all_pheno_df.txt",sep=",", header=T)%>%
  mutate(RingId_pe=RingId)


grm=readRDS("sequioa/All3085_AUTOSAUMES_RP502SNPs.RDS")
grm_filt=grm[rownames(grm)%in%tarsus_df[["RingId"]],colnames(grm)%in%tarsus_df[["RingId"]]]##filtering grm for ids only in df

grm_filt_pd <- make.positive.definite(grm_filt)
GRM <- as(grm_filt_pd, "dgCMatrix")


head(tarsus_df)
n_distinct(tarsus_df$RingId) ## 

tarsus_df$clutch_merge=as.factor(tarsus_df$clutch_merge)
tarsus_df$GeneticSex=as.factor(tarsus_df$GeneticSex)
tarsus_df$RingId=as.factor(tarsus_df$RingId)
tarsus_df$RingId_pe=as.factor(tarsus_df$RingId_pe)
tarsus_df$year=as.factor(tarsus_df$year)
tarsus_df$Observer=as.factor(tarsus_df$Observer)
tarsus_df$nestboxID=as.factor(tarsus_df$nestboxID)

tarsus_df$month=as.factor(tarsus_df$month)
tarsus_df$rank=as.numeric(tarsus_df$rank)

## priors from get_priors func .. are the same as when using GRM model
prior=c(prior(student_t(3,6.9,2.5), class = "Intercept"), ## for fixed effect of intercept setting
        prior(student_t(3,0,2.5), class = "sd", lb=0),
        prior(cauchy(0, 0.25), class = "sd", group="RingId_pe", lb=0))## for random effect of pe expect very low variance centered on 0



## running in brms 
## converges pretty well even with default number of itts and default priors
mod_tarsus.Funi <- brm(tarsus_scale ~  1 + FuniWE+GeneticSex+rank+SSgompertz(age_days, 7.17, 4.29, 0.88)+
                               (1|gr(RingId, cov=GRM))+ (1|RingId_pe)+(1|Observer)+(1|clutch_merge)+(1|year)+(1|nestboxID) + (1|month),
                             data = tarsus_df,
                             data2 = list(GRM = GRM),
                             prior=prior,
                             chains = 3,
                             cores=3 
                             # iter = 5000, 
                             # warmup = 1000
                             
)


summary(mod_tarsus.Funi)
plot(mod_tarsus.Funi)
