library(tidyverse)
library(brms)
library(viridis)
library(patchwork)
library(ggExtra)
setwd("/Users/ahewett1/Documents")

## assessing models 

ib_growth_bill_fgrm=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/3.1.bill_gr_fixedb.RDS")
summary(ib_growth_bill_fgrm)

ib_growth_bill_froh=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/3.1.bill_gr_FROH_fixedb.RDS")
summary(ib_growth_bill_froh)

## no fixed b = shows id in growth stage but model doesnt match super well
ib_growth_bill=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/superceeded/bill_unscaled_gr_25k_new_pr.RDS")
summary(ib_growth_bill)




### Mass ###
## using Fgrm worked annoyingly well 
ib_growth_mass_fgrm=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/3.2_mass_gr_fixedb.RDS")
summary(ib_growth_mass_fgrm)


## Tarsus ##

## both are pretty simualr between Froh and Fgrm - both estimate a sig effect of F in c but not elsewhere
ib_growth_tarsus_fgrm=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/3.4_tarsus_gr.RDS")
summary(ib_growth_tarsus_fgrm)

ib_growth_tarsus_froh=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/3.4_tarsus_gr_FROH_fixedb.RDS")
summary(ib_growth_tarsus_froh)


