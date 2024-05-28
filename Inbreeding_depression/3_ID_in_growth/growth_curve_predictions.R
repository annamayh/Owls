library(tidyverse)
library(brms)
library(viridis)
library(patchwork)
library(ggExtra)


setwd("/Users/ahewett1/Documents")



ib_growth_bill=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/bill_unscaled_gr_25k_new_pr.RDS")
summary(ib_growth_bill)
#plot(ib_growth_bill)


asym=fixef(ib_growth_bill, pars = "asym_Intercept")[,1]
b=fixef(ib_growth_bill, pars = "b_Intercept")[,1]
c=fixef(ib_growth_bill, pars = "c_Intercept")[,1]

f=0.15
fasym=asym+(f*fixef(ib_growth_bill, pars = "asym_FuniWE")[,1])
fb=b+(f*fixef(ib_growth_bill, pars = "b_FuniWE")[,1])
fc=c+(f*fixef(ib_growth_bill, pars = "c_FuniWE")[,1])


colinb="goldenrod1"
colave="chocolate"

id_mass=ggplot(ib_growth_bill$data, aes(x=age_days, y=BillLength))+
  geom_point(alpha=0.1, colour="black", linewidth=2)+
  xlim (0,100)+
  ylim(25,200)+
  theme_classic()+
  labs(x="Age in days", y="Mass")+
  stat_function(fun=~ asym*exp(-b*(c)^.x), colour=colave, size=1.5, xlim=c(0,88), alpha=0.8)+ ##male (intercept)
  
  stat_function(fun=~fasym*exp(-fb*(fc)^.x), colour=colinb, size=1.5, xlim=c(0,88), alpha=0.8)
## inbred at F=0.15 

id_mass




##  wing
ib_growth_wing=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/3.3_wing_unscaled_gr_25k_new_pr.RDS")
summary(ib_growth_wing)
plot(ib_growth_wing)



asym=fixef(ib_growth_wing, pars = "asym_Intercept")[,1]
b=fixef(ib_growth_wing, pars = "b_Intercept")[,1]
c=fixef(ib_growth_wing, pars = "c_Intercept")[,1]

f=0.1
fasym=asym+(f*fixef(ib_growth_wing, pars = "asym_FuniWE")[,1])
fb=b+(f*fixef(ib_growth_wing, pars = "b_FuniWE")[,1])
fc=c+(f*fixef(ib_growth_wing, pars = "c_FuniWE")[,1])


colinb="goldenrod1"
colave="chocolate"

id_mass=ggplot(ib_growth_wing$data, aes(x=age_days, y=LeftWing))+
  geom_point(alpha=0.1, colour="black", linewidth=2)+
  xlim (0,100)+
  ylim(0,500)+
  theme_classic()+
  labs(x="Age in days", y="Mass")+
  stat_function(fun=~ asym*exp(-b*(c)^.x), colour=colave, size=1.5, xlim=c(0,88), alpha=0.8)+ ##male (intercept)
  
  stat_function(fun=~fasym*exp(-fb*(fc)^.x), colour=colinb, size=1.5, xlim=c(0,88), alpha=0.8)
## inbred at F=0.15 

id_mass





