library(tidyverse)
library(brms)
library(corpcor)

##################################################################################
########################### ~~ MASS ~~ #############################################
#####################################################################################

mass_df=read.table("Inbreeding_depression_owls/pheno_df/mass_all_pheno_df.txt",sep=",", header=T)

mass_df$clutch_merge=as.factor(mass_df$clutch_merge)
mass_df$sex=as.factor(mass_df$sex)
mass_df$RingId=as.factor(mass_df$RingId)
mass_df$year=as.factor(mass_df$year)
mass_df$Observer=as.factor(mass_df$Observer)
mass_df$nestboxID=as.factor(mass_df$nestboxID)

mass_df$rank=as.numeric(mass_df$rank)


## read in GRM


## same as priors for simple model + prior for pe ~= 0 because we expect little var 
prior_mass=c(prior(student_t(3, 330, 60), class = "Intercept"), ## 
             prior(student_t(3, 0, 60), class = "sd"),
             prior(student_t(3, 0, 60), class = "sigma"),
             prior(cauchy(0, 5), class = "sd", group="RingId_pe"))


## slight trouble converging when using default number of itts so increased and using priors
mod_mass_GRM.FROH <- brm(Mass ~  1 + FHBD512gen+sex+rank+age_acc+
                           (1|RingId_pe) + (1|Observer) + (1|clutch_merge) +
                           (1|year) + (1|month) + (1|nestboxID),                         
                         data = mass_df,
                         control=list(adapt_delta=0.9),
                         chains = 4,
                         cores=4,
                         prior=prior_mass, ##
                         iter = 50000,
                         warmup = 10000,
                         thin=5
)

summary(mod_mass_GRM.FROH) ###


f16=read.table("Inbreeding_depression_owls/3k_FHBD_16gens.txt", header = T)

