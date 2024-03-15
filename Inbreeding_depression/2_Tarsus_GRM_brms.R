library(tidyverse)
library(brms)
library(corpcor)

#### TARSUS ####
## read in df of tarsus
tarsus_df=read.table("./input_dfs/tarsus_all_pheno_df.txt",sep=",", header=T)

tarsus_df[,'clutch_merge']=as.factor(tarsus_df[,'clutch_merge'])
tarsus_df[,'GeneticSex']=as.factor(tarsus_df[,'GeneticSex'])
tarsus_df[,'RingId']=as.factor(tarsus_df[,'RingId'])
tarsus_df[,'year']=as.factor(tarsus_df[,'year'])
tarsus_df[,'Observer']=as.factor(tarsus_df[,'Observer'])
tarsus_df[,'rank']=as.numeric(tarsus_df[,'rank'])


grm=read_rds("./input_dfs/All3085_AUTOSAUMES_RP502SNPs.RDS")
grm_filt=grm[rownames(grm)%in%tarsus_df[["RingId"]],colnames(grm)%in%tarsus_df[["RingId"]]]##filtering grm for ids only in df

tarsus_df[,'RingId_pe']=tarsus_df[,'RingId'] ##add permanent env variable to get h2 estimate

grm_filt_pd <- make.positive.definite(grm_filt)
GRM <- as(grm_filt_pd, "dgCMatrix")

## first 3 priors are from get_priors() function in brms 
prior=c(prior(student_t(3, 694,40), class = "Intercept"), ## for fixed effect of intercept setting mean as the rough weight
        prior(student_t(3,0,40), class = "sd"),
        prior(student_t(3,0,40), class = "sigma"),
        prior(cauchy(0, 2), class = "sd", group="RingId_pe"))## for random effect of pe expect very low variance centered on 0




mod_tarsus_all.3.grm <- brm(LeftTarsus ~  1 + FuniWE+GeneticSex+rank+SSgompertz(age_days, 717, 2.27, 0.89)+
                              (1|gr(RingId, cov=Amat))+(1|RingId_pe)+(1|Observer)+(1|clutch_merge)+(1|year),
                            prior = prior,
                            data = tarsus_df,
                            control=list(adapt_delta=0.95),
                            data2 = list(Amat = GRM),
                            chains = 4,
                            cores=4,
                            iter = 100000,
                            warmup = 25000,
                            thin=20
)



saveRDS(mod_tarsus_all.3.grm,file="./outputs/tarsus_all_GRM_Funi.RDS")

