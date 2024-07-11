library(tidyverse)
library(brms)
library(viridis)
library(patchwork)
library(ggExtra)
setwd("/Users/ahewett1/Documents")

f=0.25
colinb="goldenrod1"
colave="chocolate"

ib_growth_bill=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/bill_unscaled_gr_25k_new_pr.RDS")
summary(ib_growth_bill)
#plot(ib_growth_bill)


asym=fixef(ib_growth_bill, pars = "asym_Intercept")[,1]
b=fixef(ib_growth_bill, pars = "b_Intercept")[,1]
c=fixef(ib_growth_bill, pars = "c_Intercept")[,1]

fasym=asym+(f*fixef(ib_growth_bill, pars = "asym_FuniWE")[,1])
fb=b+(f*fixef(ib_growth_bill, pars = "b_FuniWE")[,1])
fc=c+(f*fixef(ib_growth_bill, pars = "c_FuniWE")[,1])


id_bill=ggplot(ib_growth_bill$data, aes(x=age_days, y=BillLength))+
  geom_point(alpha=0.1, colour="black", size=2)+
  xlim (0,80)+
  #ylim(25,200)+
  theme_classic()+
  labs(x="Age in days", y="Bill Length")+
  stat_function(fun=~ asym*exp(-b*(c)^.x), colour=colave, size=1.5, xlim=c(0,88), alpha=0.8)+ ##male (intercept)
  
  stat_function(fun=~fasym*exp(-fb*(fc)^.x), colour=colinb, size=1.5, xlim=c(0,88), alpha=0.8)
## inbred at F=0.15 

id_bill

#################
#####  Mass #####
##################
ib_growth_mass=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/3.2_mass_unscaled_gr_45k.RDS")
summary(ib_growth_mass)
#plot(ib_growth_mass)

asym_m=fixef(ib_growth_mass, pars = "asym_Intercept")[,1]
b_m=fixef(ib_growth_mass, pars = "b_Intercept")[,1]
c_m=fixef(ib_growth_mass, pars = "c_Intercept")[,1]

fasym_m=asym_m+(f*fixef(ib_growth_mass, pars = "asym_FuniWE")[,1])
fb_m=b_m+(f*fixef(ib_growth_mass, pars = "b_FuniWE")[,1])
fc_m=c_m+(f*fixef(ib_growth_mass, pars = "c_FuniWE")[,1])



id_mass=ggplot(ib_growth_mass$data, aes(x=age_days, y=Mass))+
  geom_point(alpha=0.1, colour="black", size=2)+
  xlim (0,80)+
  ylim(0,600)+
  theme_classic()+
  labs(x="Age in days", y="Mass")+
  stat_function(fun=~ asym_m*exp(-b_m*(c_m)^.x), colour=colave, size=1.5, xlim=c(0,88), alpha=0.8)+ ##male (intercept)
  
  stat_function(fun=~fasym_m*exp(-fb_m*(fc_m)^.x), colour=colinb, size=1.5, xlim=c(0,88), alpha=0.8)

id_mass



##################
### Wing ########
################


##  wing
ib_growth_wing=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/3.3_wing_unscaled_gr_stronger_pr_30k.RDS")
summary(ib_growth_wing)
#plot(ib_growth_wing)

asym_w=fixef(ib_growth_wing, pars = "asym_Intercept")[,1]
b_w=fixef(ib_growth_wing, pars = "b_Intercept")[,1]
c_w=fixef(ib_growth_wing, pars = "c_Intercept")[,1]

fasym_w=asym_w+(f*fixef(ib_growth_wing, pars = "asym_FuniWE")[,1])
fb_w=b_w+(f*fixef(ib_growth_wing, pars = "b_FuniWE")[,1])
#fc=c+(f*fixef(ib_growth_wing, pars = "c_FuniWE")[,1]) ## c not sig


id_wing=ggplot(ib_growth_wing$data, aes(x=age_days, y=LeftWing))+
  geom_point(alpha=0.1, colour="black", linewidth=2)+
  xlim (0,100)+
  ylim(0,300)+
  theme_classic()+
  labs(x="Age in days", y="Wing length")+
  stat_function(fun=~ asym_w*exp(-b_w*(c_w)^.x), colour=colave, size=1.5, xlim=c(0,90), alpha=0.8)+ ##male (intercept)
  stat_function(fun=~fasym_w*exp(-fb_w*(c_w)^.x), colour=colinb, size=1.5, xlim=c(0,90), alpha=0.8)
## inbred at F=0.15 

id_wing





id_bill/id_mass/id_wing+plot_annotation(tag_levels = 'A')
