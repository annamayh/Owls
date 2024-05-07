library(tidyverse)
library(brms)
library(corpcor)

##################################################################################
########################### ~~ bill ~~ #############################################
#####################################################################################

bill_df=read.table("./input_dfs/bill_all_pheno_df.txt",sep=",", header=T)%>% ##
  mutate(age_acc=184*exp(-0.99*0.932^age_days)) %>%## account for age using gompertz growth
  mutate(RingId_pe=RingId) # add permanent environment for repeated measures 

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
             prior(student_t(3,0,20), class = "sigma"),
             prior(cauchy(0, 5), class = "sd", group="RingId_pe"))


mod_bill_GRM.split_stage_h2_inter <- brm(
  
  bf(BillLength ~  1 +sex+age_acc+rank+(FuniWE*gr_stage)+# interaction of F with stage to see if we can run it all in one big model 
       (0 + gr_stage||gr(RingId, cov=Amat)) + (0 + gr_stage||RingId_pe) + (0 + gr_stage||Observer) + (0 + gr_stage||clutch_merge) +
       (0 + gr_stage||year) + (0 + gr_stage||month) + (0 + gr_stage||nestboxID),
     sigma ~ gr_stage-1), 
  
  data = bill_df,
  control=list(adapt_delta=0.95),
  data2 = list(Amat = GRM),
  chains = 4,
  cores=4,
  prior=prior_bill, ##
  iter = 50000,
  warmup = 10000,
  thin=5
)

summary(mod_bill_GRM.split_stage_h2_inter) ###

saveRDS(mod_bill_GRM.split_stage_h2_inter,file="./outputs/2_split_stages/2.1.Bill_split_stage_h2_PLUSinter.RDS") ##