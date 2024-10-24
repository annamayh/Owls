library(brms,  lib = "/users/ahewett1/R")
library(corpcor,  lib = "/users/ahewett1/R")
library(readr,  lib = "/users/ahewett1/R")

##################################################################################
########################### ~~ bill ~~ #############################################
#####################################################################################
args = commandArgs(trailingOnly = TRUE)
scratch = as.character(args[1]) 

bill_df=read.table("./input_dfs/bill_all_pheno_df.txt",sep=",", header=T)

bill_df$clutch_merge=as.factor(bill_df$clutch_merge)
bill_df$sex=as.factor(bill_df$sex)
bill_df$RingId=as.factor(bill_df$RingId)
bill_df$year=as.factor(bill_df$year)
bill_df$Observer=as.factor(bill_df$Observer)
bill_df$nestboxID=as.factor(bill_df$nestboxID)

bill_df$rank=as.numeric(bill_df$rank)


grm=read_rds("./input_dfs/All3085_AUTOSAUMES_RP502SNPs.RDS")
grm_filt=grm[rownames(grm)%in%bill_df[["RingId"]],colnames(grm)%in%bill_df[["RingId"]]]##filtering grm for ids only in df

grm_filt_pd <- make.positive.definite(grm_filt)
GRM <- as(grm_filt_pd, "dgCMatrix")

## same as priors for simple model + prior for pe ~= 0 because we expect little var 
prior_bill=c(prior(student_t(3, 180,20), class = "Intercept"), ## 
             prior(student_t(3,0,20), class = "sd"),
             prior(student_t(3,0,20), class = "sigma"),
             prior(cauchy(0, 5), class = "sd", group="RingId_pe"))


mod_bill_GRM.FROH_devstage <- brm(BillLength ~  1 + FHBD512gen*dev_stage+rank+sex+
                           (1|gr(RingId, cov=Amat)) + (1|RingId_pe) + (1|Observer) + (1|clutch_merge) +
                           (1|year) + (1|month) + (1|nestboxID),
                         data = bill_df,
                         control=list(adapt_delta=0.9),
                         data2 = list(Amat = GRM),
                         chains = 4,
                         cores=4,
                         prior=prior_bill, ##
                         iter = 35000,
                         warmup = 10000,
                         thin=5
)

summary(mod_bill_GRM.FROH_devstage) ###


saveRDS(mod_bill_GRM.FROH_devstage,file=paste0(scratch,"5.1.ID_bill_GRM_FROH_inter_DevStage.RDS") )##

