library(tidyverse)
library(brms)
library(corpcor)

#### TARSUS ####
## read in df of tarsus
tarsus_df=read.table("./input_dfs/tarsus_all_pheno_df.txt",sep=",", header=T)


tarsus_df$clutch_merge=as.factor(tarsus_df$clutch_merge)
tarsus_df$GeneticSex=as.factor(tarsus_df$GeneticSex)
tarsus_df$RingId=as.factor(tarsus_df$RingId)
tarsus_df$year=as.factor(tarsus_df$year)
tarsus_df$Observer=as.factor(tarsus_df$Observer)
tarsus_df$nestboxID=as.factor(tarsus_df$nestboxID)
tarsus_df$stage=as.factor(tarsus_df$stage)


tarsus_df$rank=as.numeric(tarsus_df$rank)
tarsus_df$tarsus_scale=as.numeric(tarsus_df$tarsus_scale)




grm=read_rds("./input_dfs/All3085_AUTOSAUMES_RP502SNPs.RDS")
grm_filt=grm[rownames(grm)%in%tarsus_df[["RingId"]],colnames(grm)%in%tarsus_df[["RingId"]]]##filtering grm for ids only in df

tarsus_df[,'RingId_pe']=tarsus_df[,'RingId'] ##add permanent env variable to get h2 estimate

grm_filt_pd <- make.positive.definite(grm_filt)
GRM <- as(grm_filt_pd, "dgCMatrix")

## first 3 priors are from get_priors() function in brms 
prior=c(prior(student_t(3,6.9,2.5), class = "Intercept"), ## for fixed effect of intercept setting mean as the rough weight
        prior(student_t(3,0,2.5), class = "sd", lb=0),
        prior(cauchy(0, 0.5), class = "sd", group="RingId_pe", lb=0))## for random effect of pe expect very low variance centered on 0




mod_tarsus_all.3.grm.split.stages <- brm(
  
    bf(tarsus_scale ~  1 + FuniWE+GeneticSex+rank+SSgompertz(age_days, 7.17, 2.27, 0.89)+stage+
    (0 + stage||gr(RingId, cov=Amat))+(0 + stage||RingId_pe)+(0 + stage||Observer)+
      (0 + stage||clutch_merge) + (0 + stage||year) + (0 + stage||nestboxID),
    sigma ~ stage-1), 
  
    prior = prior,
    data = tarsus_df,
    control=list(adapt_delta=0.95),
    data2 = list(Amat = GRM),
    chains = 4,
    cores=4,
    iter = 50000,
    warmup = 10000,
    thin=10
)



saveRDS(mod_tarsus_all.3.grm.split.stages,file="./outputs/tarsus_GRM_Funi_split_stages.RDS")

