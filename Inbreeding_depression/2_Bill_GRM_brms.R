library(tidyverse)
library(brms)
library(corpcor)


############################################################
################# ~~ BILL ~~ ################################
############################################################


bill_df=read.table("./input_dfs/bill_all_pheno_df.txt",sep=",", header=T)%>%
  select(-min_mes)


bill_df$clutch_merge=as.factor(bill_df$clutch_merge)
bill_df$GeneticSex=as.factor(bill_df$GeneticSex)
bill_df$RingId=as.factor(bill_df$RingId)
bill_df$year=as.factor(bill_df$year)
bill_df$Observer=as.factor(bill_df$Observer)
bill_df$rank=as.numeric(bill_df$rank)

## read in GRM
grm=read_rds("./input_dfs/All3085_AUTOSAUMES_RP502SNPs.RDS")
grm_filt=grm[rownames(grm)%in%bill_df[["RingId"]],colnames(grm)%in%bill_df[["RingId"]]]##filtering grm for ids only in df

bill_df[,'RingId_pe']=bill_df[,'RingId'] ##add permanent env variable to get h2 estimate

grm_filt_pd <- make.positive.definite(grm_filt)
GRM <- as(grm_filt_pd, "dgCMatrix")



prior_bill=c(prior(student_t(3, 170,18), class = "Intercept"), ## for fixed effect of intercept setting mean as the rough weight
             prior(student_t(3,0,18), class = "sd"),
             prior(student_t(3,0,18), class = "sigma"),
             prior(cauchy(0, 2), class = "sd", group="RingId_pe"))




mod_bill_all_GRM_Funi <- brm(BillLength ~  1 + FuniWE+GeneticSex+rank+SSgompertz(age_days, 184, 1.01, 0.93)+
                               (1|gr(RingId, cov=Amat))+(1|RingId_pe)+(1|Observer)+(1|clutch_merge)+(1|year),
                           data = bill_df,
                           prior=prior_bill,
                           control=list(adapt_delta=0.98),
                           data2 = list(Amat = GRM),
                           chains = 4,
                           cores=4,
                           iter = 100000,
                           warmup = 25000,
                           thin=20
                           )


saveRDS(mod_bill_all_GRM_Funi,file="./outputs/bill_all_GRM_Funi.RDS")

