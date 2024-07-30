library(brms,  lib = "/users/ahewett1/R")
library(corpcor,  lib = "/users/ahewett1/R")
library(readr,  lib = "/users/ahewett1/R")

##################################################################################
########################### ~~ Tarsus ~~ #############################################
#####################################################################################

args = commandArgs(trailingOnly = TRUE)
scratch = as.character(args[1]) 

tarsus_df=read.table("./input_dfs/tarsus_all_pheno_df.txt",sep=",", header=T)

tarsus_df$clutch_merge=as.factor(tarsus_df$clutch_merge)
tarsus_df$sex=as.factor(tarsus_df$sex)
tarsus_df$RingId=as.factor(tarsus_df$RingId)
tarsus_df$year=as.factor(tarsus_df$year)
tarsus_df$Observer=as.factor(tarsus_df$Observer)
tarsus_df$nestboxID=as.factor(tarsus_df$nestboxID)

tarsus_df$rank=as.numeric(tarsus_df$rank)


## read in GRM
grm=read_rds("./input_dfs/All3085_AUTOSAUMES_RP502SNPs.RDS")
grm_filt=grm[rownames(grm)%in%tarsus_df[["RingId"]],colnames(grm)%in%tarsus_df[["RingId"]]]##filtering grm for ids only in df

grm_filt_pd <- make.positive.definite(grm_filt)
GRM <- as(grm_filt_pd, "dgCMatrix")

## same as priors for simple model + prior for pe ~= 0 because we expect little var 
prior_tarsus=c( 
             prior(student_t(3,0,90), class = "sd"),
             prior(cauchy(0, 10), class = "sd", group="RingId_pe")
             )


mod_tarsus_GRM.split_stage_h2 <- brm(
                      
                      bf(LeftTarsus ~  1 +sex+age_acc+ # F removed as a fixed effect 
                           (0 + gr_stage||gr(RingId, cov=Amat)) + (0 + gr_stage||RingId_pe) + (0 + gr_stage||Observer) + (0 + gr_stage||clutch_merge) +
                           (0 + gr_stage||year) + (0 + gr_stage||month) + (0 + gr_stage||nestboxID) + (0 + gr_stage||rank),
                         sigma ~ gr_stage-1), 
                      
                         data = tarsus_df,
                         control=list(adapt_delta=0.97),
                         data2 = list(Amat = GRM),
                         chains = 4,
                         cores=4,
                         prior=prior_tarsus, #
                         iter = 25000,
                         warmup = 8000,
                         thin=5
)

summary(mod_tarsus_GRM.split_stage_h2) ###

saveRDS(mod_tarsus_GRM.split_stage_h2,file=paste0(scratch,"2.4.Tarsus_split_stage_h2_rank_rand.RDS")) ##