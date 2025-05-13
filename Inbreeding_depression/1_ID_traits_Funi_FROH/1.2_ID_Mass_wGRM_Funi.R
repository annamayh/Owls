.libPaths(c("/work/FAC/FBM/DEE/jgoudet/barn_owl/ahewett/R", .libPaths())) #specify library cluster path

library(brms)
library(corpcor)
library(readr)
library(dplyr)

##################################################################################
########################### ~~ MASS ~~ #############################################
#####################################################################################

mass_df=read.table("./input_dfs/mass_all_pheno_df.txt",sep=",", header=T)

length(unique(mass_df$RingId))


mass_df$clutch_merge=as.factor(mass_df$clutch_merge)
mass_df$sex=as.factor(mass_df$sex)
mass_df$RingId=as.factor(mass_df$RingId)
mass_df$year=as.factor(mass_df$year)
mass_df$Observer=as.factor(mass_df$Observer)
mass_df$nestboxID=as.factor(mass_df$nestboxID)
mass_df$rank=as.numeric(mass_df$rank)

# mean centering for age so estimates are for an individual of an ave age (not 0 days old)
mass_df$mc_age_acc <- mass_df$age_acc - mean(mass_df$age_acc)


## read in GRM
grm=read_rds("./input_dfs/All3085_AUTOSAUMES_RP502SNPs.RDS")
grm_filt=grm[rownames(grm)%in%mass_df[["RingId"]],colnames(grm)%in%mass_df[["RingId"]]]##filtering grm for ids only in df
grm_filt_pd <- make.positive.definite(grm_filt)
GRM <- as(grm_filt_pd, "dgCMatrix")

## same as priors for simple model + prior for pe ~= 0 because we expect little var 
prior_mass=c(prior(student_t(3, 330, 60), class = "Intercept"), ## 
             prior(student_t(3, 0, 60), class = "sd"),
             prior(student_t(3, 0, 60), class = "sigma"),
             prior(cauchy(0, 5), class = "sd", group="RingId_pe"))


## slight trouble converging when using default number of itts so increased and using priors
mod_mass_GRM.Funi <- brm(Mass ~  1 + FuniW + sex + rank + mc_age_acc +
                           (1|gr(RingId, cov=Amat)) + (1|RingId_pe) + (1|Observer) + (1|clutch_merge) +
                           (1|year) + (1|month) + (1|nestboxID),                         
                         data = mass_df,
                         control=list(adapt_delta=0.95),
                         data2 = list(Amat = GRM),
                         chains = 4,
                         cores=4,
                         prior=prior_mass, ##
                         iter = 50000,
                         warmup = 20000,
                         thin=5
)


saveRDS(mod_mass_GRM.Funi,file="./outputs/1_traitID_subset/1.1.ID_mass_FuniW.RDS") ##

summary(mod_mass_GRM.Funi) ###


