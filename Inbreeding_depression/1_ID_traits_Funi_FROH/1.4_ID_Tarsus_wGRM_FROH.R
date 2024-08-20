library(brms,  lib = "/users/ahewett1/R")
library(corpcor,  lib = "/users/ahewett1/R")
library(readr,  lib = "/users/ahewett1/R")

##################################################################################
########################### ~~ Tarsus LENGTH ~~ #############################################
#####################################################################################

tarsus_df=read.table("./input_dfs/tarsus_all_pheno_df.txt",sep=",", header=T)

tarsus_df$clutch_merge=as.factor(tarsus_df$clutch_merge)
tarsus_df$sex=as.factor(tarsus_df$sex)
tarsus_df$RingId=as.factor(tarsus_df$RingId)
tarsus_df$year=as.factor(tarsus_df$year)
tarsus_df$Observer=as.factor(tarsus_df$Observer)
tarsus_df$nestboxID=as.factor(tarsus_df$nestboxID)

tarsus_df$rank=as.numeric(tarsus_df$rank)

#mean centring age_acc so intercept isnt when an individual is like -23 days old
tarsus_df$mc_age_acc <- tarsus_df$age_acc - mean(tarsus_df$age_acc)


## read in GRM
grm=read_rds("./input_dfs/All3085_AUTOSAUMES_RP502SNPs.RDS")

grm_filt=grm[rownames(grm)%in%tarsus_df[["RingId"]],colnames(grm)%in%tarsus_df[["RingId"]]]##filtering grm for ids only in df

grm_filt_pd <- make.positive.definite(grm_filt)
GRM <- as(grm_filt_pd, "dgCMatrix")

## same as priors for simple model + prior for pe ~= 0 because we expect little var 
prior_tarsus=c(prior(student_t(3, 650, 30), class = "Intercept"),
             prior(student_t(3, 0, 90), class = "sd"),
             prior(student_t(3, 0, 90), class = "sigma"),
             prior(cauchy(0, 5), class = "sd", group="RingId_pe"))


## slight trouble converging when using default number of itts so increased and using priors
mod_tarsus_GRM.Froh <- brm(LeftTarsus ~  1 + FHBD512gen+sex+rank+mc_age_acc+
                           (1|gr(RingId, cov=Amat)) + (1|RingId_pe) + (1|Observer) + (1|clutch_merge) +
                           (1|year) + (1|month) + (1|nestboxID),
                         data = tarsus_df,
                         control=list(adapt_delta=0.98),
                         data2 = list(Amat = GRM),
                         chains = 4,
                         cores=4,
                         prior=prior_tarsus, 
                         iter = 35000,
                         warmup = 15000,
                         thin=5
)

summary(mod_tarsus_GRM.Froh)

saveRDS(mod_tarsus_GRM.Froh,file="./outputs/1_unscaled_traits/1.4.ID_tarsus_GRM_FROH.RDS")


