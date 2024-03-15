## inbreeding depression models using all available data ##
## not including GRM ##

library(tidyverse)
library(brms)

setwd("/Users/ahewett1/Documents")

#######################################################################
##################### ~~ TARSUS ~~ #################################### 
#####################################################################

tarsus_df=read.table("Inbreeding_depression_owls/pheno_df/tarsus_all_pheno_df.txt",sep=",", header=T)%>%
  select(-min_mes)

head(tarsus_df)
n_distinct(tarsus_df$RingId) ## 

tarsus_df$clutch_merge=as.factor(tarsus_df$clutch_merge)
tarsus_df$GeneticSex=as.factor(tarsus_df$GeneticSex)
tarsus_df$RingId=as.factor(tarsus_df$RingId)
tarsus_df$year=as.factor(tarsus_df$year)
tarsus_df$Observer=as.factor(tarsus_df$Observer)
tarsus_df$rank=as.numeric(tarsus_df$rank)

## priors from get_priors func .. are the same as when using GRM model
prior_tarsus=c(prior(student_t(3, 694,40), class = "Intercept"), ## for fixed effect of intercept setting mean as the rough weight
             prior(student_t(3,0,40), class = "sd"),
             prior(student_t(3,0,40), class = "sigma"))



## running in brms 
## converges pretty well even with default number of itts and default priors
mod_tarsus_all.1.Funi <- brm(LeftTarsus ~  1 + FuniWE+GeneticSex+rank+SSgompertz(age_days, 717, 2.27, 0.89)+
                             (1|RingId)+(1|Observer)+(1|clutch_merge)+(1|year),
                            data = tarsus_df,
                            prior=prior_tarsus,
                            chains = 4,
                            cores=4, 
                            iter = 5000, 
                            warmup = 1000
                            
)

plot(mod_tarsus_all.1.Funi)
summary(mod_tarsus_all.1.Funi)
pp_check(mod_tarsus_all.1.Funi)  # posterior predictive checks

################################################################################
## running using Froh too 
mod_tarsus_all.1.FROH <- brm(LeftTarsus ~  1 + FHBD512gen +GeneticSex+rank+SSgompertz(age_days, 717, 2.27, 0.89)+
                               (1|RingId)+(1|Observer)+(1|clutch_merge)+(1|year),
                             data = tarsus_df,
                             prior=prior_tarsus,
                             chains = 4,
                             cores=4, 
                             iter = 5000, 
                             warmup = 1000
                             
)

plot(mod_tarsus_all.1.FROH)
summary(mod_tarsus_all.1.FROH)
pp_check(mod_tarsus_all.1.FROH)  # posterior predictive checks


##################################################################################
########################### ~~ MASS ~~ #############################################
#####################################################################################
mass_df=read.table("Inbreeding_depression_owls/pheno_df/mass_all_pheno_df.txt",sep=",", header=T)%>%
  select(-min_mes)

head(mass_df)

n_distinct(mass_df$RingId) ## 

mass_df$clutch_merge=as.factor(mass_df$clutch_merge)
mass_df$GeneticSex=as.factor(mass_df$GeneticSex)
mass_df$RingId=as.factor(mass_df$RingId)
mass_df$year=as.factor(mass_df$year)
mass_df$Observer=as.factor(mass_df$Observer)
mass_df$rank=as.numeric(mass_df$rank)

prior_mass=c(prior(student_t(3, 334,65), class = "Intercept"), ## for fixed effect of intercept setting mean as the rough weight
             prior(student_t(3,0,65), class = "sd"),
             prior(student_t(3,0,65), class = "sigma"))


## slight trouble converging when using default number of itts so increased and using priors
mod_mass_all.1.Funi <- brm(Mass ~  1 + FuniWE+GeneticSex+rank+SSgompertz(age_days, 352, 4.29, 0.88)+
                               (1|RingId)+(1|Observer)+(1|clutch_merge)+(1|year),
                             data = mass_df,
                             chains = 4,
                             cores=4,
                             prior=prior_mass, 
                            iter = 15000, 
                            warmup = 5000
)


summary(mod_mass_all.1.Funi)
plot(mod_mass_all.1.Funi)

#################
### froh ###

mod_mass_all.1.FROH <- brm(Mass ~  1 + FHBD512gen+GeneticSex+rank+SSgompertz(age_days, 352, 4.29, 0.88)+
                             (1|RingId)+(1|Observer)+(1|clutch_merge)+(1|year),
                           data = mass_df,
                           chains = 4,
                           cores=4,
                           prior=prior_mass, 
                           iter = 15000, 
                           warmup = 5000
)


summary(mod_mass_all.1.FROH)
plot(mod_mass_all.1.FROH)

############################################################
################# ~~ BILL ~~ ################################
############################################################
bill_df=read.table("Inbreeding_depression_owls/pheno_df/bill_all_pheno_df.txt",sep=",", header=T)%>%
  select(-min_mes)

head(bill_df)

n_distinct(bill_df$RingId) ## 

bill_df$clutch_merge=as.factor(bill_df$clutch_merge)
bill_df$GeneticSex=as.factor(bill_df$GeneticSex)
bill_df$RingId=as.factor(bill_df$RingId)
bill_df$year=as.factor(bill_df$year)
bill_df$Observer=as.factor(bill_df$Observer)
bill_df$rank=as.numeric(bill_df$rank)


prior_bill=c(prior(student_t(3, 170,18), class = "Intercept"), ## for fixed effect of intercept setting mean as the rough weight
             prior(student_t(3,0,18), class = "sd"),
             prior(student_t(3,0,18), class = "sigma"))

## running in brms using default priors
## converges pretty well even with default number of itts
mod_bill_all.1.Funi <- brm(BillLength ~  1 + FuniWE+GeneticSex+rank+SSgompertz(age_days, 184, 1.01, 0.93)+
                             (1|RingId)+(1|Observer)+(1|clutch_merge)+(1|year),
                           data = bill_df,
                           prior=prior_bill,
                           chains = 4,
                           cores=4,
                           iter = 5000, 
                           warmup = 1000)



summary(mod_bill_all.1.Funi)
plot(mod_bill_all.1.Funi)

#### FROH ######

mod_bill_all.1.FROH <- brm(BillLength ~  1 + FHBD512gen+GeneticSex+rank+SSgompertz(age_days, 184, 1.01, 0.93)+
                             (1|RingId)+(1|Observer)+(1|clutch_merge)+(1|year),
                           data = bill_df,
                           prior=prior_bill,
                           chains = 4,
                           cores=4,
                           iter = 5000, 
                           warmup = 1000)



summary(mod_bill_all.1.FROH)
plot(mod_bill_all.1.FROH)

#################################################################################
##########################
### ~~~ WING LENGTH ~~ ###
##########################
## less acurate than tarsus but has more data 
## maybe need to include 'irritability'?? 

wing_df=read.table("Inbreeding_depression_owls/pheno_df/wing_all_pheno_df.txt",sep=",", header=T)%>%
  select(-min_mes)

head(wing_df)

n_distinct(wing_df$RingId) ## 

wing_df$clutch_merge=as.factor(wing_df$clutch_merge)
wing_df$GeneticSex=as.factor(wing_df$GeneticSex)
wing_df$RingId=as.factor(wing_df$RingId)
wing_df$year=as.factor(wing_df$year)
wing_df$Observer=as.factor(wing_df$Observer)
wing_df$rank=as.numeric(wing_df$rank)

prior_wing=c(prior(student_t(3, 196,94), class = "Intercept"), ## for fixed effect of intercept setting mean as the rough weight
        prior(student_t(3,0,94), class = "sd"),
        prior(student_t(3,0,94), class = "sigma"))

## struggled to converge with default
mod_wing_all.1.Funi <- brm(LeftWing ~  1 + FuniWE+GeneticSex+rank+SSgompertz(age_days, 317, 4.19, 0.95)+
                             (1|RingId)+(1|Observer)+(1|clutch_merge)+(1|year),
                           prior=prior_wing,
                           data = wing_df,
                           chains = 4,
                           cores=4, 
                           iter = 15000, 
                           warmup = 5000
)

summary(mod_wing_all.1.Funi)
plot(mod_wing_all.1.Funi)


###### froh #######

mod_wing_all.1.FROH <- brm(LeftWing ~  1 + FHBD512gen+GeneticSex+rank+SSgompertz(age_days, 317, 4.19, 0.95)+
                             (1|RingId)+(1|Observer)+(1|clutch_merge)+(1|year),
                           prior=prior_wing,
                           data = wing_df,
                           chains = 4,
                           cores=4, 
                           iter = 15000, 
                           warmup = 5000
)

summary(mod_wing_all.1.FROH)
plot(mod_wing_all.1.FROH)


### save all model outputs 

# save(mod_tarsus_all.1.Funi, mod_tarsus_all.1.FROH, tarsus_df, 
#      mod_mass_all.1.Funi, mod_mass_all.1.FROH, mass_df, 
#      mod_bill_all.1.Funi, mod_bill_all.1.FROH, bill_df, 
#      mod_wing_all.1.Funi, mod_wing_all.1.FROH, wing_df,
#      file = "Inbreeding_depression_owls/Model_outputs/brms_funi_froh_simple_all.RData")
# 


