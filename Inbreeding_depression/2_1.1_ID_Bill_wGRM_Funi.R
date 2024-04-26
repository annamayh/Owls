library(tidyverse)
library(brms)
library(corpcor)

##################################################################################
########################### ~~ bill ~~ #############################################
#####################################################################################

bill_df=read.table("./input_dfs/bill_all_pheno_df.txt",sep=",", header=T)%>% ##
  mutate(age_acc=1.84*exp(-0.99*0.932^age_days)) ## account for age using gompertz growth

bill_df$clutch_merge=as.factor(bill_df$clutch_merge)
bill_df$GeneticSex=as.factor(bill_df$GeneticSex)
bill_df$RingId=as.factor(bill_df$RingId)
bill_df$year=as.factor(bill_df$year)
bill_df$Observer=as.factor(bill_df$Observer)
bill_df$nestboxID=as.factor(bill_df$nestboxID)


bill_df$rank=as.numeric(bill_df$rank)
bill_df$CH1903X=as.numeric(bill_df$CH1903X) # locations of the nestboxes
bill_df$CH1903Y=as.numeric(bill_df$CH1903Y)


## read in GRM
grm=read_rds("./input_dfs/All3085_AUTOSAUMES_RP502SNPs.RDS")

grm_filt=grm[rownames(grm)%in%bill_df[["RingId"]],colnames(grm)%in%bill_df[["RingId"]]]##filtering grm for ids only in df

bill_df[,'RingId_pe']=bill_df[,'RingId'] ##add permanent env variable to get h2 estimate

grm_filt_pd <- make.positive.definite(grm_filt)
GRM <- as(grm_filt_pd, "dgCMatrix")

## same as priors for simple model + prior for pe ~= 0 because we expect little var 
prior_bill=c(prior(student_t(3, 1.8,2.5), class = "Intercept"), ## 
             prior(student_t(3,0,2.5), class = "sd"),
             prior(student_t(3,0,2.5), class = "sigma"),
             prior(cauchy(0, 0.5), class = "sd", group="RingId_pe"))


## slight trouble converging when using default number of itts so increased and using priors
mod_bill_GRM.Funi <- brm(bill_scale ~  1 + FuniWE+GeneticSex+rank+age_acc+CH1903X+CH1903Y+
                           (1|gr(RingId, cov=Amat))+(1|RingId_pe)+(1|Observer)+(1|clutch_merge)+(1|year)+(1|nestboxID),
                         data = bill_df,
                         control=list(adapt_delta=0.85),
                         data2 = list(Amat = GRM),
                         chains = 4,
                         cores=4,
                         prior=prior_bill, ##
                         iter = 65000,
                         warmup = 15000,
                         thin=5
)

summary(mod_bill_GRM.Funi) ###

saveRDS(mod_bill_GRM.Funi,file="./outputs/1.1.ID_bill_GRM_Funi.RDS") ##