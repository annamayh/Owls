
library(tidyverse)
library(brms)
library(corpcor)



#### TARSUS ####
## read in df of tarsus 
fledge_tarsus_df=read.table("./input_dfs/fledge_pheno_tarsus_df.txt",sep=",", header=T)

fledge_tarsus_df$clutch_merge=as.factor(fledge_tarsus_df$clutch_merge)
fledge_tarsus_df$GeneticSex=as.factor(fledge_tarsus_df$GeneticSex)
fledge_tarsus_df$RingId=as.factor(fledge_tarsus_df$RingId)
fledge_tarsus_df$year=as.factor(fledge_tarsus_df$year)
fledge_tarsus_df$Observer=as.factor(fledge_tarsus_df$Observer)
fledge_tarsus_df$rank=as.numeric(fledge_tarsus_df$rank)


grm=read_rds("./input_dfs/All3085_AUTOSAUMES_RP502SNPs.RDS")
grm_filt=grm[rownames(grm)%in%fledge_tarsus_df[["RingId"]],colnames(grm)%in%fledge_tarsus_df[["RingId"]]]##filtering grm for ids only in df

fledge_tarsus_df$RingId_pe=fledge_tarsus_df$RingId ##add permanent env variable to get h2 estimate

grm_filt_pd <- make.positive.definite(grm_filt) 
GRM <- as(grm_filt_pd, "dgCMatrix")


c(prior(normal(1,2), class = b, coef = FuniWE),
    prior(normal(1,2), class = b, coef = GeneticSex),
    prior(normal(1,2), class = b, coef = rank),
    prior(normal(1,2), class = b, coef = treat),
    prior(normal(1,2), class = b, coef = treat))


mod_tarsus.3.grm <- brm(LeftTarsus ~  1 + FuniWE+GeneticSex+rank+SSgompertz(age_days, 719.5, 2.274, 0.89)+
                          (1|gr(RingId, cov=Amat))+(1|RingId_pe)+(1|Observer)+(1|clutch_merge)+(1|year),
                        data = fledge_tarsus_df,
                        data2 = list(Amat = GRM),
                        chains = 4, 
                        cores=4,
                        iter = 10000, 
                        warmup = 2500, 
                        thin=1
)



save(mod_tarsus.3.grm,file="./outputs/tarsus_GRM_Funi.RData")
