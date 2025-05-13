.libPaths(c("/work/FAC/FBM/DEE/jgoudet/barn_owl/ahewett/R", .libPaths())) #specify library cluster path1-11

library(tidyverse)
library(brms)
library(corpcor)

##################################################################################
########################### ~~ bill ~~ #############################################
#####################################################################################

args = commandArgs(trailingOnly = TRUE)
scratch = as.character(args[1]) 

bill_df=read.table("./input_dfs/bill_all_pheno_df.txt",sep=",", header=T) ##
  
bill_df$clutch_merge=as.factor(bill_df$clutch_merge)
bill_df$sex=as.factor(bill_df$sex)
bill_df$RingId=as.factor(bill_df$RingId)
bill_df$year=as.factor(bill_df$year)
bill_df$Observer=as.factor(bill_df$Observer)
bill_df$nestboxID=as.factor(bill_df$nestboxID)

bill_df$rank=as.numeric(bill_df$rank)


## read in GRM
grm=read_rds("./input_dfs/All3085_AUTOSAUMES_RP502SNPs.RDS")
grm_filt=grm[rownames(grm)%in%bill_df[["RingId"]],colnames(grm)%in%bill_df[["RingId"]]]##filtering grm for ids only in df

grm_filt_pd <- make.positive.definite(grm_filt)
GRM <- as(grm_filt_pd, "dgCMatrix")

## same as priors for simple model + prior for pe ~= 0 because we expect little var 
prior_bill=c(prior(student_t(3, 180,20), class = "Intercept"), ## 
             prior(student_t(3,0,20), class = "sd"),
             prior(normal(0, 5), class = "sd", group="RingId_pe"))

mod_bill_GRM.split_stage_h2 <- brm(
                      
                      bf(BillLength ~  1 +sex+age_acc+ # F removed as a fixed effect 
                           (0 + gr_stage||gr(RingId, cov=Amat)) + (0 + gr_stage||RingId_pe) + (0 + gr_stage||Observer) + (0 + gr_stage||clutch_merge) +
                           (0 + gr_stage||year) + (0 + gr_stage||month) + (0 + gr_stage||nestboxID) + (0 + gr_stage||rank),
                         sigma ~ gr_stage-1), 
                      
                         data = bill_df,
                         control=list(adapt_delta=0.95),
                         data2 = list(Amat = GRM),
                         chains = 4,
                         cores=4,
                         prior=prior_bill, ##
                         iter = 35000,
                         warmup = 5000,
                         thin=5
)

summary(mod_bill_GRM.split_stage_h2) ###

saveRDS(mod_bill_GRM.split_stage_h2,file=paste0(scratch,"2.1.bill_split_stage_subset.RDS")) ##