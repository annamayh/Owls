library(brms,  lib = "/users/ahewett1/R")
library(corpcor,  lib = "/users/ahewett1/R")
library(readr,  lib = "/users/ahewett1/R")

##################################################################################
########################### ~~ Tarsus LENGTH ~~ #############################################
#####################################################################################

wing_df=read.table("./input_dfs/wing_all_pheno_df.txt",sep=",", header=T)

wing_df$clutch_merge=as.factor(wing_df$clutch_merge)
wing_df$sex=as.factor(wing_df$sex)
wing_df$RingId=as.factor(wing_df$RingId)
wing_df$year=as.factor(wing_df$year)
wing_df$Observer=as.factor(wing_df$Observer)
wing_df$nestboxID=as.factor(wing_df$nestboxID)

wing_df$rank=as.numeric(wing_df$rank)


## read in GRM
grm=read_rds("./input_dfs/All3085_AUTOSAUMES_RP502SNPs.RDS")

grm_filt=grm[rownames(grm)%in%wing_df[["RingId"]],colnames(grm)%in%wing_df[["RingId"]]]##filtering grm for ids only in df

grm_filt_pd <- make.positive.definite(grm_filt)
GRM <- as(grm_filt_pd, "dgCMatrix")

## same as priors for simple model + prior for pe ~= 0 because we expect little var 
prior_wing=c(prior(student_t(3, 200, 90), class = "Intercept"), ## 
             prior(student_t(3, 0, 90), class = "sd"),
             prior(student_t(3, 0, 90), class = "sigma"),
             prior(cauchy(0, 5), class = "sd", group="RingId_pe"))


## slight trouble converging when using default number of itts so increased and using priors
mod_wing_GRM.Funi <- brm(LeftWing ~  1 + FuniWE+sex+rank+age_acc+
                           (1|gr(RingId, cov=Amat)) + (1|RingId_pe) + (1|Observer) + (1|clutch_merge) +
                           (1|year) + (1|month) + (1|nestboxID),
                         data = wing_df,
                         control=list(adapt_delta=0.85),
                         data2 = list(Amat = GRM),
                         chains = 4,
                         cores=4,
                         prior=prior_wing, 
                         iter = 50000,
                         warmup = 10000,
                         thin=5
)

summary(mod_wing_GRM.Funi)

saveRDS(mod_wing_GRM.Funi,file="./outputs/1.3.ID_wing_GRM_Funi_unscaled.RDS")


