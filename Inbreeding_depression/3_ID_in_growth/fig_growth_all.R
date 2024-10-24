library(tidyverse)
library(brms)
library(viridis)
library(patchwork)
library(grid)
library(ggExtra)
setwd("/Users/ahewett1/Documents")

f=0.5
alpha=0.01
dot_col='black'
colinb="orangered"
colave="palegreen3"
col_none='black'


## assessing models and plotting ##

## Bill ##
## atm neither froh or fgrm show ID in any growth
## fgrm doesnt converge super well -  maybe need to re-run
ib_growth_bill_fgrm=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/3.1.bill_gr_fixedb.RDS")
summary(ib_growth_bill_fgrm)


asym=fixef(ib_growth_bill_fgrm, pars = "asym_Intercept")[,1]
b=fixef(ib_growth_bill_fgrm, pars = "b_Intercept")[,1]
c=fixef(ib_growth_bill_fgrm, pars = "c_Intercept")[,1]

(id_bill_fgrm=ggplot(ib_growth_bill_fgrm$data, aes(x=age_days, y=(BillLength*0.1)))+
  geom_point(alpha=alpha, colour=dot_col, size=2)+
  xlim (0,80)+
  ylim(5,20)+
  theme_classic()+
    theme(axis.title.x=element_blank(),
          #axis.title.y=element_blank()
          )+
  labs(x="Age in days", y="Bill Length (mm)")+
  stat_function(fun=~ (asym*0.1)*exp(-b*(c)^.x), colour=col_none, size=1.5, xlim=c(0,88), alpha=0.8)+ ##male (intercept)
  # (* 0.1 to get bill length in mm)
  stat_function(fun=~(asym*0.1)*exp(-b*(c)^.x), colour=col_none, size=1.5, xlim=c(0,88), alpha=0.8))
  


# converges but shows no affect of inbreeding depression 
ib_growth_bill_froh=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/3.1.bill_gr_FROH_totalfixedb.RDS")
summary(ib_growth_bill_froh)

#b because its fixed eff not random 
post=as_draws_df(ib_growth_bill_froh, variable = c('b_asym_Intercept', 'b_b_Intercept', 'b_c_Intercept'))%>%
  slice_sample(n = 5000)


asym_fr=fixef(ib_growth_bill_froh, pars = "asym_Intercept")[,1]
b_fr=fixef(ib_growth_bill_froh, pars = "b_Intercept")[,1]
c_fr=fixef(ib_growth_bill_froh, pars = "c_Intercept")[,1]

(id_bill_froh=ggplot(ib_growth_bill_froh$data, aes(x=age_days, y=(BillLength*0.1)))+
  geom_point(alpha=alpha, colour=dot_col, size=2)+
  xlim (0,80)+
  ylim(5,20)+
  theme_classic()+
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank())+
  #labs(x="Age in days", y="Bill Length (mm)")+
  stat_function(fun=~ (asym_fr*0.1)*exp(-b_fr*(c_fr)^.x), colour=col_none, size=1.5, xlim=c(0,88), alpha=0.8)+ ##male (intercept)
  stat_function(fun=~(asym_fr*0.1)*exp(-b_fr*(c_fr)^.x), colour=col_none, size=1.5, xlim=c(0,88), alpha=0.8))


# for (i in 1:100) {
#   # Extract parameters from the i-th draw
#   asym <- post$b_asym_Intercept[i]
#   b <- post$b_b_Intercept[i]
#   c <- post$b_c_Intercept[i]
#   
#   # Add a thin line for the i-th draw
#   id_bill_froh <- id_bill_froh +
#     stat_function(fun = ~ (asym * 0.1) * exp(-b * (c)^.x), colour = "yellow")
# }
# 
# id_bill_froh


### Mass ###
# using fgrm shows id in asymptote #using froh shows no id (on the threshold for growth)
ib_growth_mass_fgrm=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/3.2_mass_gr_fixedb.RDS")
summary(ib_growth_mass_fgrm)

asym_m=fixef(ib_growth_mass_fgrm, pars = "asym_Intercept")[,1]
b_m=fixef(ib_growth_mass_fgrm, pars = "b_Intercept")[,1]
c_m=fixef(ib_growth_mass_fgrm, pars = "c_Intercept")[,1]

fasym_m=asym_m+(f*fixef(ib_growth_mass_fgrm, pars = "asym_FuniWE")[,1])

(id_mass_fgrm=ggplot(ib_growth_mass_fgrm$data, aes(x=age_days, y=Mass))+
  geom_point(alpha=alpha, colour=dot_col, size=2)+
  xlim (0,80)+
  ylim(0,500)+
  theme_classic()+
    theme(axis.title.x=element_blank(),
          #axis.title.y=element_blank()
          )+
  labs(x="Age in days", y="Mass (g)")+
  stat_function(fun=~ asym_m*exp(-b_m*(c_m)^.x), colour=colave, size=1.5, xlim=c(0,88), alpha=0.8)+ # line for average male id
  stat_function(fun=~fasym_m*exp(-b_m*(c_m)^.x), colour=colinb, size=1.5, xlim=c(0,88), alpha=0.8)) #line for inbred id
## nothing showing for froh
ib_growth_mass_froh=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/3.2_mass_gr_FROH_totalfixedb.RDS")
summary(ib_growth_mass_froh)

asym_m_fr=fixef(ib_growth_mass_froh, pars = "asym_Intercept")[,1]
b_m_fr=fixef(ib_growth_mass_froh, pars = "b_Intercept")[,1]
c_m_fr=fixef(ib_growth_mass_froh, pars = "c_Intercept")[,1]

(id_mass_froh=ggplot(ib_growth_mass_froh$data, aes(x=age_days, y=Mass))+
  geom_point(alpha=alpha, colour=dot_col, size=2)+
  xlim (0,80)+
  ylim(0,500)+
  theme_classic()+
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank())+
  #labs(x="Age in days", y="Mass (g)")+
  stat_function(fun=~ asym_m_fr*exp(-b_m_fr*(c_m_fr)^.x), colour=col_none, size=1.5, xlim=c(0,88), alpha=0.8)+ ##male (intercept)
  stat_function(fun=~asym_m_fr*exp(-b_m_fr*(c_m_fr)^.x), colour=col_none, size=1.5, xlim=c(0,88), alpha=0.8))





### tarsus ###
# fgrm shows no id - doesnt in lmm too
ib_growth_tarsus_fgrm=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/3.4_tarsus_gr.RDS")
summary(ib_growth_tarsus_fgrm)

asym1_t=fixef(ib_growth_tarsus_fgrm, pars = "asym1_Intercept")[,1]
asym_t=fixef(ib_growth_tarsus_fgrm, pars = "asym_Intercept")[,1]
b_t=fixef(ib_growth_tarsus_fgrm, pars = "b_Intercept")[,1]
c_t=fixef(ib_growth_tarsus_fgrm, pars = "c_Intercept")[,1]


(id_tarsus_fgrm=ggplot(ib_growth_tarsus_fgrm$data, aes(x=age_days, y=LeftTarsus))+
  geom_point(alpha=alpha, colour=dot_col, linewidth=2)+
  xlim (0,70)+
  ylim(100,800)+
  theme_classic()+
   # theme(axis.title.x=element_blank(),
    #      axis.title.y=element_blank())+
  labs(x="Age in days", y="Tarsus length (mm)")+
  stat_function(fun=~ asym1_t + asym_t*exp(-b_t*(c_t)^.x), colour=col_none, size=1.5, xlim=c(0,80), alpha=0.8)+ ##male (intercept)
  stat_function(fun=~ asym1_t +asym_t*exp(-b_t*(c_t)^.x), colour=col_none, size=1.5, xlim=c(0,80), alpha=0.8))


#froh shows id in c
ib_growth_tarsus_froh=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/3.4_tarsus_gr_FROH_fixedb.RDS")
summary(ib_growth_tarsus_froh)


asym1_t_fr=fixef(ib_growth_tarsus_froh, pars = "asym1_Intercept")[,1]
asym_t_fr=fixef(ib_growth_tarsus_froh, pars = "asym_Intercept")[,1]
b_t_fr=fixef(ib_growth_tarsus_froh, pars = "b_Intercept")[,1]
c_t_fr=fixef(ib_growth_tarsus_froh, pars = "c_Intercept")[,1]

fc_t_fr=c_t_fr+(f*fixef(ib_growth_tarsus_froh, pars = "c_FHBD512gen")[,1]) ## c sig


(id_tarsus_froh=ggplot(ib_growth_tarsus_froh$data, aes(x=age_days, y=LeftTarsus))+
    geom_point(alpha=alpha, colour=dot_col, linewidth=2)+
    theme_classic()+
    theme(
          axis.title.y=element_blank())+
    xlim (0,70)+
    ylim(100,800)+
    labs(x="Age in days", y="Tarsus length (mm)")+
    stat_function(fun=~ asym1_t_fr + asym_t_fr*exp(-b_t_fr*(c_t_fr)^.x), colour=colave, size=1.5, xlim=c(0,80), alpha=0.8)+ ##male (intercept)
    stat_function(fun=~ asym1_t_fr +asym_t_fr*exp(-b_t_fr*(fc_t_fr)^.x), colour=colinb, size=1.5, xlim=c(0,80), alpha=0.8))
    



# row_label_1 <- wrap_elements(panel = textGrob('First Row', rot=90))
# row_label_2 <- wrap_elements(panel = textGrob('First Row', rot=90))
# 
# 
# (id_bill_froh+id_mass_froh+id_tarsus_froh)/
#   (id_bill_fgrm+id_mass_fgrm+id_tarsus_fgrm)+plot_annotation(tag_levels = 'a')



col_label_1 <- wrap_elements(panel = textGrob(expression('F'[uniWE])))
col_label_2 <- wrap_elements(panel = textGrob(expression('F'[ROH])))


plot_all=((col_label_1+col_label_2)/
  (id_bill_fgrm+id_bill_froh)/
  (id_mass_fgrm+id_mass_froh)/
  (id_tarsus_fgrm+id_tarsus_froh)+
  plot_layout( 
              heights = c(0.5, 5, 5, 5)))




ggsave(plot_all,
       file = "Inbreeding_depression_owls/Model_outputs/3_growth/growth_plots_all_froh_fgrm.png",
       width = 6,
       height = 6)

