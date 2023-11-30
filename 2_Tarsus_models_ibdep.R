## inbreeding depression in tarsus length ###

library(MCMCglmm)
library(tidyverse)
library(brms)
library(pedigree)

setwd("D:/")
## read in df of tarsus 

head(fledge_tarsus_df)

n_distinct(fledge_tarsus_df$RingId) ## 

fledge_tarsus_df$clutch_merge=as.factor(fledge_tarsus_df$clutch_merge)
fledge_tarsus_df$GeneticSex=as.factor(fledge_tarsus_df$GeneticSex)
fledge_tarsus_df$RingId=as.factor(fledge_tarsus_df$RingId)
fledge_tarsus_df$year=as.factor(fledge_tarsus_df$year)
fledge_tarsus_df$Observer=as.factor(fledge_tarsus_df$Observer)
fledge_tarsus_df$rank=as.numeric(fledge_tarsus_df$rank)

## find rough growth curve for all data
fm1 <- nls(LeftTarsus ~ SSgompertz(age_days, asym, b2, b3),
           data = fledge_tarsus_df)
summary(fm1)
# Estimate Std. Error t value Pr(>|t|)    
# asym 7.195e+02  7.005e-01 1027.11   <2e-16 ***
#   b2   2.274e+00  2.635e-02   86.29   <2e-16 ***
#   b3   8.897e-01  7.375e-04 1206.40   <2e-16 ***

plot(fledge_tarsus_df$age_days, fledge_tarsus_df$LeftTarsus)
curve(7.195e+02*exp(-2.274*8.897e-01^x), from = 0, to=80, add = TRUE, col="red", lwd=2)


## first simple glmm using fROH as inbreeding coeff
prior1.2 <- list(G = list(G1 = list(V=1,nu=1,alpha.mu=0,alpha.V=1000), ##after accounting for age etc then the variance bt ids should be v small 
                          G2 = list(V=1,nu=0.02),
                          G3 = list(V=1,nu=0.02),
                          G4 = list(V=1,nu=0.02)), 
                 R = list(V=1,nu=0.02))

## FHBD512gen == FROH
mod_tarsus.1 <- MCMCglmm(LeftTarsus ~  1 + FHBD512gen+GeneticSex+rank+SSgompertz(age_days, 719.5, 2.274, 0.89),
                   random = ~RingId+Observer+clutch_merge+year, 
                   data = fledge_tarsus_df,
                   prior = prior1.2, 
                   nitt = 110000, ##100k iters probably a bit of overkill tbh, can probs get by on 55k and 5k burnin
                   # thin = 50, 
                   burnin = 10000)


plot(mod_tarsus.1)
summary(mod_tarsus.1)

save(mod_tarsus.1, file="Owls_temp/Outputs_and_df_inbdep/tarsus_output_model.RData")

################################
## ~~  now with pedigree ~~  ###
################################

ped=read.table("Owls_temp/Outputs_and_df_inbdep/pedigreeCORRECTED.tab", header = T)%>%
  select(-sex)%>%
  rename(RingId=id)

ids_in<-as.matrix(fledge_tarsus_df%>%select(RingId))
pruned<-MCMCglmm::prunePed(ped, ids_in)
ord<-orderPed(pruned)
ped_ordered=pruned[order(ord),]

ped_ordered$RingId <- as.factor(ped_ordered$RingId)
ped_ordered$dadid <- as.factor(ped_ordered$dadid)
ped_ordered$momid <- as.factor(ped_ordered$momid)

Ainv<-inverseA(ped_ordered, nodes="ALL")$Ainv

fledge_tarsus_df$RingId_pe=fledge_tarsus_df$RingId ##add permanent env variable to get h2 estimate

prior1.3 <- list(G = list(G1 = list(V=1,nu=1,alpha.mu=0,alpha.V=1000), ##after accounting for age etc then the variance bt ids should be v small 
                          G2 = list(V=1,nu=1,aplha.mu=0,alpha.V=1000),
                          G3 = list(V=1,nu=0.02),
                          G4 = list(V=1,nu=0.02),
                          G5 = list(V=1,nu=0.02)), 
                 R = list(V=1,nu=0.02))


mod_tarsus.2.ped <- MCMCglmm(LeftTarsus ~  1 + FHBD512gen+GeneticSex+rank+SSgompertz(age_days, 719.5, 2.274, 0.89),
                         random = ~RingId+RingId_pe+Observer+clutch_merge+year, 
                         data = fledge_tarsus_df,
                         ginv = list(RingId = Ainv),
                         prior = prior1.3, 
                         nitt = 110000, 
                         ##thin = 50, ##next time thinning at 50 
                         burnin = 10000)

plot(mod_tarsus.2.ped)
summary(mod_tarsus.2.ped)

save(mod_tarsus.2.ped, file="Owls_temp/Outputs_and_df_inbdep/tarsus_output_modelPED.RData")


## ~~ and finally with GRM fitted in brms ~~ ###

grm=read_rds("swisstransfer_9dccae57-94d3-4b48-a0b7-ecf0b75014b3/All3085_AUTOSAUMES_RP502SNPs.RDS")
grm_filt=grm[rownames(grm)%in%fledge_mes_rep[["RingId"]],colnames(grm)%in%fledge_mes_rep[["RingId"]]]##filtering grm for ids only in df

library(corpcor)
grm_filt_pd <- make.positive.definite(grm_filt) 
GRM <- as(grm_filt_pd, "dgCMatrix")

mod_tarsus.3.grm <- brm(LeftTarsus ~  1 + FuniWE+GeneticSex+rank+SSgompertz(age_days, 719.5, 2.274, 0.89)+
                    (1|gr(RingId, cov=Amat))+(1|RingId_pe)+(1|Observer)+(1|clutch_merge)+(1|year),
                  data = fledge_tarsus_df,
                  data2 = list(Amat = GRM),
                  chains = 4, 
                  cores=2,
                  iter = 1500
)


##### maybe compare Funi and FROH ?? ####

