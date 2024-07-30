library(viridis)
library(patchwork)
library(ggExtra)
library(brms)
library(tidyverse)

## Basic inbreeding depression models accounting for age ... using 2 types of inbreeding coefficient: FROH and Funi ##
## 3 traits investigated: bill, mass and wing length

setwd("/Users/ahewett1/Documents")
## BILL##
id_bill_funi=readRDS("Inbreeding_depression_owls/Model_outputs/1_unscaled_ID_wGRM/1.1.ID_bill_GRM_Funi_unscaled.RDS")
#id_bill_froh=readRDS("Inbreeding_depression_owls/Model_outputs/unscaled_ID_wGRM/1.1.ID_bill_GRM_Funi_unscaled.RDS")
## MASS ##
id_mass_funi=readRDS("Inbreeding_depression_owls/Model_outputs/unscaled_ID_wGRM/1.2.ID_mass_GRM_Funi_unscaled.RDS")
id_mass_froh=readRDS("Inbreeding_depression_owls/Model_outputs/unscaled_ID_wGRM/1.2.ID_mass_GRM_FROH_unscaled.RDS")
## WING ##
id_wing_funi=readRDS("Inbreeding_depression_owls/Model_outputs/unscaled_ID_wGRM/1.3.ID_wing_GRM_Funi_unscaled.RDS")
id_wing_froh=readRDS("Inbreeding_depression_owls/Model_outputs/unscaled_ID_wGRM/1.3.ID_wing_GRM_FROH_unscaled.RDS")


summary(id_bill_funi)
bill_data=read.table("Inbreeding_depression_owls/pheno_df/bill_all_pheno_df.txt",sep=",", header=T)

age_30day_bill=184*exp(-0.99*0.932^30) ##

x=184*exp(-0.99*0.932^-23.5) ##


f=0.25

est_bill_outbred=age_30day_bill+5.45-0.19+(0*-15.73)
est_bill_inbred=age_30day_bill+5.45-0.19+(f*-15.73)





wing_data=read.table("Inbreeding_depression_owls/pheno_df/wing_all_pheno_df.txt",sep=",", header=T)
summary(id_wing_funi)

age_30day_wing=298*exp(-4.19*0.945^3)## account for age unsing gompertz growth

f=0.25

est_wing_outbred=age_30day_wing+1.25-1.14+(0*-9.32)
est_wing_inbred=age_30day_wing+1.25-1.14+(f*-9.32)



mass_data=read.table("Inbreeding_depression_owls/pheno_df/mass_all_pheno_df.txt",sep=",", header=T)


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

