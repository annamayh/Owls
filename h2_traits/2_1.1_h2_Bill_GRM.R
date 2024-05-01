library(tidyverse)
library(brms)
library(corpcor)

##################################################################################
########################### ~~ bill ~~ #############################################
#####################################################################################
setwd("/Users/ahewett1/Documents")

bill_df=read.table("Gen_arch_owls/pheno_df/bill_all_pheno_df.txt",sep=",", header=T)%>% ##
  mutate(age_acc=184*exp(-0.99*0.932^age_days)) %>% ## account for age using gompertz growth
  mutate(RingId_pe=RingId)

bill_df$clutch_merge=as.factor(bill_df$clutch_merge)
bill_df$sex=as.factor(bill_df$sex)
bill_df$RingId=as.factor(bill_df$RingId)
bill_df$year=as.factor(bill_df$year)
bill_df$Observer=as.factor(bill_df$Observer)
bill_df$nestboxID=as.factor(bill_df$nestboxID)


bill_df$julian_hatchdate=as.numeric(bill_df$julian_hatchdate)
bill_df$rank=as.numeric(bill_df$rank)
bill_df$month=as.numeric(bill_df$month)


## read in GRM
grm=read_rds("sequioa/All3085_AUTOSAUMES_RP502SNPs.RDS")

grm_filt=grm[rownames(grm)%in%bill_df[["RingId"]],colnames(grm)%in%bill_df[["RingId"]]]##filtering grm for ids only in df

grm_filt_pd <- make.positive.definite(grm_filt)
GRM <- as(grm_filt_pd, "dgCMatrix")


## same as priors for simple model + prior for pe ~= 0 because we expect little var 
prior_bill=c(prior(student_t(3, 180,20), class = "Intercept"), ## 
             prior(student_t(3,0,20), class = "sd"),
             prior(student_t(3,0,20), class = "sigma"),
             prior(cauchy(0, 5), class = "sd", group="RingId_pe"))


h2_mod_bill_GRM <- get_prior(BillLength ~  1 +sex+age_acc+
                         (1|RingId_pe)+(1|gr(RingId, cov=Amat))+(1|Observer)+(1|rank)+(1|nestboxID)+
                         (1|clutch_merge)+(1|year)+(1|month)+(1|MumId)+(1|DadId),
                         data = bill_df,
                         control=list(adapt_delta=0.85),
                         data2 = list(Amat = GRM),
                         chains = 2,
                         cores=2,
                         prior=prior_bill ##
                         # iter = 65000,
                         # warmup = 15000,
                         # thin=5
)

summary(h2_mod_bill_GRM) ###

saveRDS(h2_mod_bill_GRM,file="./outputs/1.1.ID_bill_GRM.RDS") ##