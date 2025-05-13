.libPaths(c("/work/FAC/FBM/DEE/jgoudet/barn_owl/ahewett/R", .libPaths())) #specify library cluster path

library(brms)
library(corpcor)
library(readr)
library(dplyr)

##################################################################################
########################### ~~ Bill ~~ #############################################
#####################################################################################
## inbreeding depression estimates on bill length ##
# using FuniW and fitting GRM as random
# also mean centering age (so estimates are for an ave aged indiv)  

bill_df=read.table("./input_dfs/bill_all_pheno_df.txt",sep=",", header=T)

length(unique(bill_df$RingId))

# sorting out variables 
bill_df$clutch_merge=as.factor(bill_df$clutch_merge)
bill_df$sex=as.factor(bill_df$sex)
bill_df$RingId=as.factor(bill_df$RingId)
bill_df$year=as.factor(bill_df$year)
bill_df$Observer=as.factor(bill_df$Observer)
bill_df$nestboxID=as.factor(bill_df$nestboxID)
bill_df$gr_stage=as.factor(bill_df$gr_stage)
bill_df$rank=as.numeric(bill_df$rank)

# mean centering for age so estimates are for an individual of an ave age (not 0 days old)
bill_df$mc_age_acc <- bill_df$age_acc - mean(bill_df$age_acc)

## read in GRM
grm=read_rds("./input_dfs/All3085_AUTOSAUMES_RP502SNPs.RDS")
grm_filt=grm[rownames(grm)%in%bill_df[["RingId"]],colnames(grm)%in%bill_df[["RingId"]]]##filtering grm for ids only in df

grm_filt_pd <- make.positive.definite(grm_filt)
GRM <- as(grm_filt_pd, "dgCMatrix")

## same as priors for simple model + prior for pe ~= 0 because we expect little var 
prior_bill=c(prior(student_t(3, 180,20), class = "Intercept"), ## 
             prior(student_t(3,0,20), class = "sd"),
             prior(student_t(3,0,20), class = "sigma"),
             prior(cauchy(0, 5), class = "sd", group="RingId_pe"))


## slight trouble converging when using default number of itts so increased and using priors
mod_bill_GRM.Funi <- brm(BillLength ~  1 + FuniW + sex + mc_age_acc + rank +
                           (1|gr(RingId, cov=Amat)) + (1|RingId_pe) + (1|Observer) + (1|clutch_merge) +
                           (1|year) + (1|month) + (1|nestboxID),
                         data = bill_df,
                         control=list(adapt_delta=0.95),
                         data2 = list(Amat = GRM),
                         chains = 4,
                         cores=4,
                         prior=prior_bill, ##
                         iter = 25000,
                         warmup = 7500,
                         thin=5
)

saveRDS(mod_bill_GRM.Funi,file="./outputs/1_traitID_subset/1.1.ID_bill_FuniW.RDS") ##

summary(mod_bill_GRM.Funi) ###
