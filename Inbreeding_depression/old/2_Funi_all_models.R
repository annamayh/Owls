
## Main models for Elo's chapter ####
## Decided on Funi as inbreeding coefficient 
## one simple model and one including the pedigree (to estimate h2)
## will also run same 2 models for FROH but only for supp

library(MCMCglmm)
library(tidyverse)
library(brms)
library(pedigree)

<<<<<<< HEAD
setwd("/Users/ahewett1/Documents")

#### TARSUS ####

## read in df of tarsus 
fledge_tarsus_df=read.table("Inbreeding_depression_owls/pheno_df/tarsus_fledge_pheno_df.txt",sep=",", header=T)
=======
setwd("D:/")
setwd("/Volumes/Seagate Portable Drive")
#### TARSUS ####

## read in df of tarsus 
fledge_tarsus_df=read.table("Owls_temp/fledge_pheno_tarsus_df.txt",sep=",", header=T)%>%
  unique()

>>>>>>> a4d6f8cb4427ae105de4d7196205bb961ceb8f0d
head(fledge_tarsus_df)

n_distinct(fledge_tarsus_df$RingId) ## 

fledge_tarsus_df$clutch_merge=as.factor(fledge_tarsus_df$clutch_merge)
fledge_tarsus_df$GeneticSex=as.factor(fledge_tarsus_df$GeneticSex)
fledge_tarsus_df$RingId=as.factor(fledge_tarsus_df$RingId)
fledge_tarsus_df$year=as.factor(fledge_tarsus_df$year)
fledge_tarsus_df$Observer=as.factor(fledge_tarsus_df$Observer)
fledge_tarsus_df$rank=as.numeric(fledge_tarsus_df$rank)

## find rough growth curve for all data
gomp_gr <- nls(LeftTarsus ~ SSgompertz(age_days, asym, b2, b3),
           data = fledge_tarsus_df)
summary(gomp_gr)
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



#################################
## ### Fgrm simple model ### ###
#################################
mod_tarsus.1.2 <- MCMCglmm(LeftTarsus ~  1 + FuniWE +GeneticSex+rank+SSgompertz(age_days, 719.5, 2.274, 0.89),
                           random = ~RingId+Observer+clutch_merge+year, 
                           data = fledge_tarsus_df,
                           prior = prior1.2, 
                           nitt = 80000, 
                           burnin = 10000,
                           #thin=20
                           )


plot(mod_tarsus.1.2)
summary(mod_tarsus.1.2)

growth_log=function(x){
  
  720-(720*exp(-0.08*x))
  
}


mod_tarsus.1.2.2 <- MCMCglmm(LeftTarsus ~  1 + FuniWE +GeneticSex+rank+growth_log(age_days),
                           random = ~RingId+Observer+clutch_merge+year, 
                           data = fledge_tarsus_df,
                           prior = prior1.2, 
                           nitt = 80000, 
                           burnin = 10000,
                           #thin=20
)


plot(mod_tarsus.1.2.2)
summary(mod_tarsus.1.2.2)

### DIC lower using gompertz ###
# DIC: 61794.25 
# DIC: 63896.57 


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


mod_tarsus.2.ped <- MCMCglmm(LeftTarsus ~  1 + FuniWE+GeneticSex+rank+SSgompertz(age_days, 719.5, 2.274, 0.89),
                             random = ~RingId+RingId_pe+Observer+clutch_merge+year, 
                             data = fledge_tarsus_df,
                             ginv = list(RingId = Ainv),
                             prior = prior1.3, 
                             nitt = 250000, 
                             thin = 20,  
                             burnin = 50000)

plot(mod_tarsus.2.ped)
summary(mod_tarsus.2.ped)

autocorr.diag(mod_tarsus.2.ped$VCV)

##################################################################################
### MASS ###

fledge_mass_df=read.table("Owls_temp/fledge_pheno_mass_df.txt",sep=",", header=T)

head(fledge_mass_df)
n_distinct(fledge_mass_df$RingId) ## 

fledge_mass_df$clutch_merge=as.factor(fledge_mass_df$clutch_merge)
fledge_mass_df$GeneticSex=as.factor(fledge_mass_df$GeneticSex)
fledge_mass_df$RingId=as.factor(fledge_mass_df$RingId)
fledge_mass_df$year=as.factor(fledge_mass_df$year)
fledge_mass_df$Observer=as.factor(fledge_mass_df$Observer)
fledge_mass_df$rank=as.numeric(fledge_mass_df$rank)

## find rough growth curve for all data
fm2 <- nls(Mass ~ SSgompertz(age_days, asym, b2, b3),
           data = fledge_mass_df)
summary(fm2)


# Formula: Mass ~ SSgompertz(age_days, asym, b2, b3)
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# asym 3.695e+02  7.555e-01  489.02   <2e-16 ***
#   b2   4.425e+00  1.174e-01   37.68   <2e-16 ***
#   b3   8.906e-01  1.348e-03  660.85   <2e-16 ***


#################################
## ### Fgrm simple model ### ###
#################################
mod_mass.1.2 <- MCMCglmm(Mass ~  1 + FuniWE +GeneticSex+rank+SSgompertz(age_days, 3.695e+02, 4.425,8.906e-01),
                           random = ~RingId+Observer+clutch_merge+year, 
                           data = fledge_mass_df,
                           prior = prior1.2, 
                           nitt = 80000, 
                           burnin = 10000,
                          # thin=50
                         )


plot(mod_mass.1.2)
summary(mod_mass.1.2)


ids_in<-as.matrix(fledge_mass_df%>%select(RingId))
pruned<-MCMCglmm::prunePed(ped, ids_in)
ord<-orderPed(pruned)
ped_ordered=pruned[order(ord),]

ped_ordered$RingId <- as.factor(ped_ordered$RingId)
ped_ordered$dadid <- as.factor(ped_ordered$dadid)
ped_ordered$momid <- as.factor(ped_ordered$momid)

Ainv<-inverseA(ped_ordered, nodes="ALL")$Ainv

fledge_mass_df$RingId_pe=fledge_mass_df$RingId ##add permanent env variable to get h2 estimate

mod_mass.2.ped <- MCMCglmm(Mass ~  1 + FuniWE+GeneticSex+rank+SSgompertz(age_days, 3.695e+02, 4.425,8.906e-01),
                             random = ~RingId+RingId_pe+Observer+clutch_merge+year, 
                             data = fledge_mass_df,
                             ginv = list(RingId = Ainv),
                             prior = prior1.3, 
                             nitt = 250000, 
                             thin = 20,  
                             burnin = 50000)

plot(mod_mass.2.ped)
summary(mod_mass.2.ped)

###############################################################################
### Bill Length ########

fledge_bill_df=read.table("Owls_temp/fledge_pheno_bill_df.txt",sep=",", header=T)

head(fledge_bill_df)
n_distinct(fledge_bill_df$RingId) ## 

fledge_bill_df$clutch_merge=as.factor(fledge_bill_df$clutch_merge)
fledge_bill_df$GeneticSex=as.factor(fledge_bill_df$GeneticSex)
fledge_bill_df$RingId=as.factor(fledge_bill_df$RingId)
fledge_bill_df$year=as.factor(fledge_bill_df$year)
fledge_bill_df$Observer=as.factor(fledge_bill_df$Observer)
fledge_bill_df$rank=as.numeric(fledge_bill_df$rank)

## find rough growth curve for all data
fm2 <- nls(BillLength ~ SSgompertz(age_days, asym, b2, b3),
           data = fledge_bill_df)
summary(fm2)

plot(fledge_bill_df$age_days, fledge_bill_df$BillLength)
curve(1.798e+02*exp(-1.127* 9.205e-01^x), from = 0, to=80, add = TRUE, col="red", lwd=2)


mod_bill.1.2 <- MCMCglmm(BillLength ~  1 + FuniWE +GeneticSex+rank+SSgompertz(age_days, 1.798e+02, 1.127, 9.205e-01),
                           random = ~RingId+Observer+clutch_merge+year, 
                           data = fledge_bill_df,
                           prior = prior1.2, 
                           nitt = 80000, 
                           burnin = 10000,
                           #thin=20
)


plot(mod_bill.1.2)
summary(mod_bill.1.2)



ids_in<-as.matrix(fledge_bill_df%>%select(RingId))
pruned<-MCMCglmm::prunePed(ped, ids_in)
ord<-orderPed(pruned)
ped_ordered=pruned[order(ord),]

ped_ordered$RingId <- as.factor(ped_ordered$RingId)
ped_ordered$dadid <- as.factor(ped_ordered$dadid)
ped_ordered$momid <- as.factor(ped_ordered$momid)

Ainv<-inverseA(ped_ordered, nodes="ALL")$Ainv

fledge_bill_df$RingId_pe=fledge_bill_df$RingId ##add permanent env variable to get h2 estimate

mod_bill.2.ped <- MCMCglmm(BillLength ~  1 + FuniWE+GeneticSex+rank+SSgompertz(age_days, 1.798e+02, 1.127, 9.205e-01),
                           random = ~RingId+RingId_pe+Observer+clutch_merge+year, 
                           data = fledge_bill_df,
                           ginv = list(RingId = Ainv),
                           prior = prior1.3, 
                           nitt = 250000, 
                           thin = 20,  
                           burnin = 50000)

plot(mod_bill.2.ped)
summary(mod_bill.2.ped)



save(mod_tarsus.1.2, mod_tarsus.2.ped,
     mod_mass.1.2, mod_mass.2.ped,
     mod_bill.1.2,mod_bill.2.ped,
     file="Owls_temp/Outputs_and_df_inbdep/Funi_all_model_outputs.RData")

