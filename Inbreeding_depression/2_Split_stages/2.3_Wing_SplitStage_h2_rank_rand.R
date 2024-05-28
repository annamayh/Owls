library(brms,  lib = "/users/ahewett1/R")
library(corpcor,  lib = "/users/ahewett1/R")
library(readr,  lib = "/users/ahewett1/R")

##################################################################################
########################### ~~ wing ~~ #############################################
#####################################################################################

args = commandArgs(trailingOnly = TRUE)
scratch = as.character(args[1]) 

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
prior_wing=c(prior(student_t(3, 200,90), class = "Intercept"), ## 
             prior(student_t(3,0,90), class = "sd"),
             prior(cauchy(0, 10), class = "sd", group="RingId_pe"))


mod_wing_GRM.split_stage_h2 <- brm(
                      
                      bf(LeftWing ~  1 +sex+age_acc+ # F removed as a fixed effect 
                           (0 + gr_stage||gr(RingId, cov=Amat)) + (0 + gr_stage||RingId_pe) + (0 + gr_stage||Observer) + (0 + gr_stage||clutch_merge) +
                           (0 + gr_stage||year) + (0 + gr_stage||month) + (0 + gr_stage||nestboxID) + (0 + gr_stage||rank),
                         sigma ~ gr_stage-1), 
                      
                         data = wing_df,
                         control=list(adapt_delta=0.97),
                         data2 = list(Amat = GRM),
                         chains = 4,
                         cores=4,
                         prior=prior_wing, #
                         iter = 25000,
                         warmup = 8000,
                         thin=5
)

summary(mod_wing_GRM.split_stage_h2) ###

saveRDS(mod_wing_GRM.split_stage_h2,file=paste0(scratch,"2.3.Wing_split_stage_h2_rank_rand.RDS")) ##