library(viridis)
library(patchwork)
library(ggExtra)
library(brms)
library(tidyverse)

## Basic inbreeding depression models accounting for age ... using 2 types of inbreeding coefficient: FROH and Funi ##
## 3 traits investigated: bill, mass and wing length

setwd("/Users/ahewett1/Documents")
## BILL##
id_bill_funi=readRDS("Inbreeding_depression_owls/Model_outputs/unscaled_ID_wGRM/1.1.ID_bill_GRM_Funi_unscaled.RDS")
#id_bill_froh=readRDS("Inbreeding_depression_owls/Model_outputs/unscaled_ID_wGRM/1.1.ID_bill_GRM_Funi_unscaled.RDS")
## MASS ##
id_mass_funi=readRDS("Inbreeding_depression_owls/Model_outputs/unscaled_ID_wGRM/1.2.ID_mass_GRM_Funi_unscaled.RDS")
id_mass_froh=readRDS("Inbreeding_depression_owls/Model_outputs/unscaled_ID_wGRM/1.2.ID_mass_GRM_FROH_unscaled.RDS")
## WING ##
id_wing_funi=readRDS("Inbreeding_depression_owls/Model_outputs/unscaled_ID_wGRM/1.3.ID_wing_GRM_Funi_unscaled.RDS")
id_wing_froh=readRDS("Inbreeding_depression_owls/Model_outputs/unscaled_ID_wGRM/1.3.ID_wing_GRM_FROH_unscaled.RDS")


summary(id_wing_funi)

Fs_to_predict=data.frame(
  sex = factor(c("1","1","1","1","1", "2","2","2","2","2")),
  age_acc = c(30),
  FuniWE = c(0,0.05,0.1,0.15,0.2),
  rank = 1
  
)


pp=predict(id_wing_funi, 
           newdata = Fs_to_predict,
           re_formula = NULL)%>%
  as.data.frame()

pp$sex=c("1","1","1","1","1", "2","2","2","2","2")
pp$F=c(0,0.05,0.1,0.15,0.2,0,0.05,0.1,0.15,0.2)



ggplot(pp, aes(x=F, y=Estimate, color=sex))+
  geom_line()+
  geom_ribbon(aes(ymin=Q2.5, ymax=Q97.5), alpha = 0.1)+
  theme(text = element_text(size = 18))

