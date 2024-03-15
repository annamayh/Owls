library(tidyverse)
library(brms)
library(corpcor)

##################################################################################
########################### ~~ MASS ~~ #############################################
#####################################################################################

mass_df=read.table("./input_dfs/mass_all_pheno_df.txt",sep=",", header=T)%>%
  select(-min_mes)

mass_df$clutch_merge=as.factor(mass_df$clutch_merge)
mass_df$GeneticSex=as.factor(mass_df$GeneticSex)
mass_df$RingId=as.factor(mass_df$RingId)
mass_df$year=as.factor(mass_df$year)
mass_df$Observer=as.factor(mass_df$Observer)
mass_df$rank=as.numeric(mass_df$rank)

## read in GRM
grm=read_rds("./input_dfs/All3085_AUTOSAUMES_RP502SNPs.RDS")
grm_filt=grm[rownames(grm)%in%tarsus_df[["RingId"]],colnames(grm)%in%tarsus_df[["RingId"]]]##filtering grm for ids only in df

tarsus_df[,'RingId_pe']=tarsus_df[,'RingId'] ##add permanent env variable to get h2 estimate

grm_filt_pd <- make.positive.definite(grm_filt)
GRM <- as(grm_filt_pd, "dgCMatrix")

## same as priors for simple model + prior for pe ~= 0 because we expect little var 
prior_mass=c(prior(student_t(3, 334,65), class = "Intercept"), ## 
             prior(student_t(3,0,65), class = "sd"),
             prior(student_t(3,0,65), class = "sigma"),
             prior(cauchy(0, 2), class = "sd", group="RingId_pe"))


## slight trouble converging when using default number of itts so increased and using priors
mod_mass_all_GRM.1.Funi <- brm(Mass ~  1 + FuniWE+GeneticSex+rank+SSgompertz(age_days, 352, 4.29, 0.88)+
                             (1|gr(RingId, cov=Amat))+(1|RingId_pe)+(1|Observer)+(1|clutch_merge)+(1|year),
                           data = mass_df,
                           control=list(adapt_delta=0.95),
                           data2 = list(Amat = GRM),
                           chains = 4,
                           cores=4,
                           prior=prior_mass, 
                           iter = 100000,
                           warmup = 25000,
                           thin=20
)


saveRDS(mod_mass_all_GRM.1.Funi,file="./outputs/mass_all_GRM_Funi.RDS")
