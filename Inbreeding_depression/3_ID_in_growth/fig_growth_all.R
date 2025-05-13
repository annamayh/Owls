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
colinb="darkorange1"
colave="palegreen3"
col_none='black'


## assessing models and plotting ##
###############################################################
                    ## Bill ##
###############################################################
## atm neither froh or fgrm show ID in any growth
## fgrm doesnt converge super well -  maybe need to re-run
ib_growth_bill_fgrm=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth_subset/3.1.bill_gr_Funi_subset.RDS")
summary(ib_growth_bill_fgrm)

asym=fixef(ib_growth_bill_fgrm, pars = "asym_Intercept")[,1]
b=fixef(ib_growth_bill_fgrm, pars = "b_Intercept")[,1]
c=fixef(ib_growth_bill_fgrm, pars = "c_Intercept")[,1]

fc=c+(f*(fixef(ib_growth_bill_fgrm, pars = "c_FuniWE")[,1]))


##  getting upper and lower CIs for average id
asym_l95=fixef(ib_growth_bill_fgrm, pars = "asym_Intercept")[,3]
asym_u95=fixef(ib_growth_bill_fgrm, pars = "asym_Intercept")[,4]
b_l95=fixef(ib_growth_bill_fgrm, pars = "b_Intercept")[,3]
b_u95=fixef(ib_growth_bill_fgrm, pars = "b_Intercept")[,4]
c_l95=fixef(ib_growth_bill_fgrm, pars = "c_Intercept")[,3]
c_u95=fixef(ib_growth_bill_fgrm, pars = "c_Intercept")[,4]
##  getting upper and lower CIs for inbred id
asym_l95_f=asym+(f*(fixef(ib_growth_bill_fgrm, pars = "asym_FuniWE")[,3]))
asym_u95_f=asym+(f*(fixef(ib_growth_bill_fgrm, pars = "asym_FuniWE")[,4]))
c_l95_f=c+(f*(fixef(ib_growth_bill_fgrm, pars = "c_FuniWE")[,3]))
c_u95_f=c+(f*(fixef(ib_growth_bill_fgrm, pars = "c_FuniWE")[,4]))


(id_bill_fgrm=ggplot(ib_growth_bill_fgrm$data, aes(x=age_days, y=(BillLength*0.1)))+
  geom_point(alpha=alpha, colour=dot_col, size=2)+
  xlim (0,80)+
  ylim(5,20)+
  theme_classic()+
    theme(text = element_text(size = 18)#axis.title.x=element_blank(),
          #axis.title.y=element_blank()
          )+
  labs(x="Age in days", y="Bill Length (mm)", tag = 'C')+
  stat_function(fun=~ (asym*0.1)*exp(-b*(c)^.x), colour=colave, size=1.5, xlim=c(0,88), alpha=0.8)+ ##male (intercept)
  # (* 0.1 to get bill length in mm)
  stat_function(fun=~(asym*0.1)*exp(-b*(fc)^.x), colour=colinb, size=1.5, xlim=c(0,88), alpha=0.8)+
    #upper and lower CIs for 
    stat_function(fun=~ (asym_l95*0.1)*exp(-b_l95*(c_u95)^.x), colour=colave, xlim=c(0,88), linetype='dashed')+
    stat_function(fun=~ (asym_u95*0.1)*exp(-b_u95*(c_l95)^.x), colour=colave, xlim=c(0,88), linetype='dashed')+
    # and CIs for inbred ids - use the lower aymptote and upper CI for growth rate (cos lower number =better growing)
    stat_function(fun=~ (asym_l95_f*0.1)*exp(-b*(c_u95_f)^.x), colour=colinb, xlim=c(0,88), linetype='dashed')+
    stat_function(fun=~ (asym_u95_f*0.1)*exp(-b*(c_l95_f)^.x), colour=colinb, xlim=c(0,88), linetype='dashed')

  )
  

############
### FROH ####
############
# converges but shows no affect of inbreeding depression 
ib_growth_bill_froh=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth_subset/3.1.bill_gr_FROH_subset.RDS")
summary(ib_growth_bill_froh)

(((184+4)-184)/184)*100


## when there is no difference only plot the main intercept line as black
asym_fr=fixef(ib_growth_bill_froh, pars = "asym_Intercept")[,1]
b_fr=fixef(ib_growth_bill_froh, pars = "b_Intercept")[,1]
c_fr=fixef(ib_growth_bill_froh, pars = "c_Intercept")[,1]
##  getting upper and lower CIs for average id
asym_fr_l95=fixef(ib_growth_bill_froh, pars = "asym_Intercept")[,3]
asym_fr_u95=fixef(ib_growth_bill_froh, pars = "asym_Intercept")[,4]
b_fr_l95=fixef(ib_growth_bill_froh, pars = "b_Intercept")[,3]
b_fr_u95=fixef(ib_growth_bill_froh, pars = "b_Intercept")[,4]
c_fr_l95=fixef(ib_growth_bill_froh, pars = "c_Intercept")[,3]
c_fr_u95=fixef(ib_growth_bill_froh, pars = "c_Intercept")[,4]
##  getting upper and lower CIs for inbred id
asym_fr_l95_f=asym_fr+(f*(fixef(ib_growth_bill_froh, pars = "asym_FHBD512gen")[,3]))
asym_fr_u95_f=asym_fr+(f*(fixef(ib_growth_bill_froh, pars = "asym_FHBD512gen")[,4]))
c_fr_l95_f=c_fr+(f*(fixef(ib_growth_bill_froh, pars = "c_FHBD512gen")[,3]))
c_fr_u95_f=c_fr+(f*(fixef(ib_growth_bill_froh, pars = "c_FHBD512gen")[,4]))



(id_bill_froh=ggplot(ib_growth_bill_froh$data, aes(x=age_days, y=(BillLength*0.1)))+
  geom_point(alpha=alpha, colour=dot_col, size=2)+
  xlim (0,80)+
  ylim(5,20)+
  theme_classic()+
    theme(text = element_text(size = 18)#axis.title.x=element_blank(),
         # axis.title.y=element_blank()
          )+
  labs(x="Age in days", y="Bill Length (mm)", tag = 'F')+
  stat_function(fun=~ (asym_fr*0.1)*exp(-b_fr*(c_fr)^.x), colour=col_none, size=1.5, xlim=c(0,88), alpha=0.8)+ ##male (intercept)
  #stat_function(fun=~(asym_fr*0.1)*exp(-b_fr*(c_fr)^.x), colour=col_none, size=1, xlim=c(0,88), alpha=0.8)+
  #upper and lower CIs for 
  stat_function(fun=~ (asym_fr_l95*0.1)*exp(-b_fr_l95*(c_fr_l95)^.x), colour=colave, xlim=c(0,88), linetype='dashed')+
  stat_function(fun=~ (asym_fr_u95*0.1)*exp(-b_fr_u95*(c_fr_u95)^.x), colour=colave, xlim=c(0,88), linetype='dashed')+
  # and CIs for inbred ids - use the lower aymptote and upper CI for growth rate (cos lower number =better growing)
  stat_function(fun=~ (asym_fr_l95_f*0.1)*exp(-b_fr*(c_fr_u95_f)^.x), colour=colinb, xlim=c(0,88), linetype='dashed')+
  stat_function(fun=~ (asym_fr_u95_f*0.1)*exp(-b_fr*(c_fr_l95_f)^.x), colour=colinb, xlim=c(0,88), linetype='dashed')
    
)



###############################################################
                ### Mass ###
###############################################################

# using fgrm shows id in asymptote #using froh shows no id (on the threshold for growth)
ib_growth_mass_fgrm=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth_subset/3.2_mass_gr_Funi_subset.RDS")

ib_growth_mass_fgrm2=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/3.2_mass_gr_fixedb.RDS")

summary(ib_growth_mass_fgrm)
summary(ib_growth_mass_fgrm2)


asym_m=fixef(ib_growth_mass_fgrm, pars = "asym_Intercept")[,1]
b_m=fixef(ib_growth_mass_fgrm, pars = "b_Intercept")[,1]
c_m=fixef(ib_growth_mass_fgrm, pars = "c_Intercept")[,1]

fasym_m=asym_m+(f*fixef(ib_growth_mass_fgrm, pars = "asym_FuniWE")[,1])

## upper and lower 95% CIs for intercept 
asym_m_l95=fixef(ib_growth_mass_fgrm, pars = "asym_Intercept")[,3]
asym_m_u95=fixef(ib_growth_mass_fgrm, pars = "asym_Intercept")[,4]
b_m_l95=fixef(ib_growth_mass_fgrm, pars = "b_Intercept")[,3]
b_m_u95=fixef(ib_growth_mass_fgrm, pars = "b_Intercept")[,4]
c_m_l95=fixef(ib_growth_mass_fgrm, pars = "c_Intercept")[,3]
c_m_u95=fixef(ib_growth_mass_fgrm, pars = "c_Intercept")[,4]

# CIs for inbred 
fasym_m_l95=asym_m+(f*fixef(ib_growth_mass_fgrm, pars = "asym_FuniWE")[,3])
fasym_m_u95=asym_m+(f*fixef(ib_growth_mass_fgrm, pars = "asym_FuniWE")[,4])
fc_m_l95=c_m+(f*fixef(ib_growth_mass_fgrm, pars = "c_FuniWE")[,3])
fc_m_u95=c_m+(f*fixef(ib_growth_mass_fgrm, pars = "c_FuniWE")[,4])


(id_mass_fgrm=ggplot(ib_growth_mass_fgrm$data, aes(x=age_days, y=Mass))+
  geom_point(alpha=alpha, colour=dot_col, size=2)+
  xlim (0,80)+
  ylim(0,500)+
  theme_classic()+
    theme(text = element_text(size = 18)#axis.title.x=element_blank(),
          #axis.title.y=element_blank()
          )+
  labs(x="Age in days", y="Mass (g)", tag = 'D')+
  stat_function(fun=~ asym_m*exp(-b_m*(c_m)^.x), colour=colave, size=1, xlim=c(0,88))+ # line for average male id
  stat_function(fun=~fasym_m*exp(-b_m*(c_m)^.x), colour=colinb, size=1, xlim=c(0,88))+
  
  stat_function(fun=~ asym_m_u95*exp(-b_m_u95*(c_m_l95)^.x), colour=colave, linetype='dashed',alpha=0.8, xlim=c(0,88))+ # 
  stat_function(fun=~ asym_m_l95*exp(-b_m_l95*(c_m_u95)^.x), colour=colave, linetype='dashed', alpha=0.8, xlim=c(0,88)) +# 
    
  stat_function(fun=~ fasym_m_u95*exp(-b_m*(fc_m_l95)^.x), colour=colinb, linetype='dashed',alpha=0.8, xlim=c(0,88))+ # 
  stat_function(fun=~ fasym_m_l95*exp(-b_m*(fc_m_u95)^.x), colour=colinb, linetype='dashed', alpha=0.8, xlim=c(0,88)) #
  
  ) 


############
### FROH ####
############

## nothing showing for froh
ib_growth_mass_froh=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth_subset/3.2_mass_gr_FROH_subset.RDS")

ib_growth_mass_froh2=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/3.2_mass_gr_FROH_totalfixedb.RDS")


summary(ib_growth_mass_froh)
summary(ib_growth_mass_froh2)


asym_m_fr=fixef(ib_growth_mass_froh, pars = "asym_Intercept")[,1]
b_m_fr=fixef(ib_growth_mass_froh, pars = "b_Intercept")[,1]
c_m_fr=fixef(ib_growth_mass_froh, pars = "c_Intercept")[,1]

## CIs for intercept 
asym_m_fr_l95=fixef(ib_growth_mass_froh, pars = "asym_Intercept")[,3]
asym_m_fr_u95=fixef(ib_growth_mass_froh, pars = "asym_Intercept")[,4]
b_m_fr_l95=fixef(ib_growth_mass_froh, pars = "b_Intercept")[,3]
b_m_fr_u95=fixef(ib_growth_mass_froh, pars = "b_Intercept")[,4]
c_m_fr_l95=fixef(ib_growth_mass_froh, pars = "c_Intercept")[,3]
c_m_fr_u95=fixef(ib_growth_mass_froh, pars = "c_Intercept")[,4]
# CIs for inbred
fasym_m_fr_l95=asym_m_fr+(f*(fixef(ib_growth_mass_froh, pars = "asym_FHBD512gen")[,3]))
fasym_m_fr_u95=asym_m_fr+(f*(fixef(ib_growth_mass_froh, pars = "asym_FHBD512gen")[,4]))
fc_m_fr_l95=c_m_fr+(f*(fixef(ib_growth_mass_froh, pars = "c_FHBD512gen")[,3]))
fc_m_fr_u95=c_m_fr+(f*(fixef(ib_growth_mass_froh, pars = "c_FHBD512gen")[,4]))





(id_mass_froh=ggplot(ib_growth_mass_froh$data, aes(x=age_days, y=Mass))+
  geom_point(alpha=alpha, colour=dot_col, size=2)+
  xlim (0,80)+
  ylim(0,500)+
  theme_classic()+
    theme(text = element_text(size = 18)#axis.title.x=element_blank(),
          #axis.title.y=element_blank()
          )+
  labs(x="Age in days", y="Mass (g)", tag = 'G')+
  stat_function(fun=~ asym_m_fr*exp(-b_m_fr*(c_m_fr)^.x), colour=col_none, size=1, xlim=c(0,88), alpha=0.8)+ ##male (intercept)
  stat_function(fun=~asym_m_fr_u95*exp(-b_m_fr_u95*(c_m_fr_l95)^.x), colour=colave,linetype='dashed',  xlim=c(0,88))+
  stat_function(fun=~asym_m_fr_l95*exp(-b_m_fr_l95*(c_m_fr_u95)^.x), colour=colave, linetype='dashed', xlim=c(0,88))+
  stat_function(fun=~fasym_m_fr_u95*exp(-b_m_fr*(fc_m_fr_l95)^.x), colour=colinb, linetype='dashed', xlim=c(0,88))+
  stat_function(fun=~asym_m_fr_l95*exp(-b_m_fr*(fc_m_fr_u95)^.x), colour=colinb, linetype='dashed', xlim=c(0,88))
  
  )

(((372+26)-372)/372)*100


###############################################################
                  ### Tarsus ###
###############################################################
# fgrm shows no id - doesnt in lmm too
ib_growth_tarsus_fgrm=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth_subset/3.4_tarsus_gr_Funi_subset.RDS")
summary(ib_growth_tarsus_fgrm)

asym1_t=fixef(ib_growth_tarsus_fgrm, pars = "asym1_Intercept")[,1]
asym_t=fixef(ib_growth_tarsus_fgrm, pars = "asym_Intercept")[,1]
b_t=fixef(ib_growth_tarsus_fgrm, pars = "b_Intercept")[,1]
c_t=fixef(ib_growth_tarsus_fgrm, pars = "c_Intercept")[,1]
## CIs for intercept
asym1_t_l95=fixef(ib_growth_tarsus_fgrm, pars = "asym1_Intercept")[,3]
asym1_t_u95=fixef(ib_growth_tarsus_fgrm, pars = "asym1_Intercept")[,4]
asym_t_l95=fixef(ib_growth_tarsus_fgrm, pars = "asym_Intercept")[,3]
asym_t_u95=fixef(ib_growth_tarsus_fgrm, pars = "asym_Intercept")[,4]
b_t_l95=fixef(ib_growth_tarsus_fgrm, pars = "b_Intercept")[,3]
b_t_u95=fixef(ib_growth_tarsus_fgrm, pars = "b_Intercept")[,4]
c_t_l95=fixef(ib_growth_tarsus_fgrm, pars = "c_Intercept")[,3]
c_t_u95=fixef(ib_growth_tarsus_fgrm, pars = "c_Intercept")[,4]
## CIs for inbred
fasym_t_l95=asym_t+(f*(fixef(ib_growth_tarsus_fgrm, pars = "asym_FuniWE")[,3]))
fasym_t_u95=asym_t+(f*(fixef(ib_growth_tarsus_fgrm, pars = "asym_FuniWE")[,4]))
fc_t_l95=c_t+(f*(fixef(ib_growth_tarsus_fgrm, pars = "c_FuniWE")[,3]))
fc_t_u95=c_t+(f*(fixef(ib_growth_tarsus_fgrm, pars = "c_FuniWE")[,4]))


(id_tarsus_fgrm=ggplot(ib_growth_tarsus_fgrm$data, aes(x=age_days, y=LeftTarsus*0.1))+
  geom_point(alpha=alpha, colour=dot_col, linewidth=2)+
  xlim (0,70)+
  ylim(10,80)+
  theme_classic()+
   theme(text = element_text(size = 18))+#axis.title.x=element_blank(),
    #      axis.title.y=element_blank())+
  labs(x="Age in days", y="Tarsus length (mm)", tag = 'E')+
  stat_function(fun=~ (asym1_t*0.1) + (asym_t*0.1)*exp(-b_t*(c_t)^.x), colour=col_none, size=1, xlim=c(0,80), alpha=0.8)+ ##male (intercept)
 # stat_function(fun=~ asym1_t +asym_t*exp(-b_t*(c_t)^.x), colour=col_none, size=1.5, xlim=c(0,80), alpha=0.8)+
    #ave CIs
    stat_function(fun=~ (asym1_t_u95*0.1) +(asym_t_u95*0.1)*exp(-b_t_u95*(c_t_l95)^.x), colour=colave, linetype='dashed', xlim=c(0,80))+
    stat_function(fun=~ (asym1_t_l95*0.1) + (asym_t_l95*0.1)*exp(-b_t_l95*(c_t_u95)^.x), colour=colave,  linetype='dashed', xlim=c(0,80))+
    #inbred CIs  
    stat_function(fun=~ (asym1_t*0.1)+(fasym_t_l95*0.1)*exp(-b_t*(fc_t_u95)^.x), colour=colinb,  linetype='dashed', xlim=c(0,80))+
    stat_function(fun=~ (asym1_t*0.1) +(fasym_t_u95*0.1)*exp(-b_t*(fc_t_l95)^.x), colour=colinb,  linetype='dashed', xlim=c(0,80))

    )

############
### FROH ###
############
ib_growth_tarsus_froh=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth_subset/3.4_tarsus_gr_FROH_subset.RDS")
summary(ib_growth_tarsus_froh)

(((587-3.16)-587)/587)*100


asym1_t_fr=fixef(ib_growth_tarsus_froh, pars = "asym1_Intercept")[,1]
asym_t_fr=fixef(ib_growth_tarsus_froh, pars = "asym_Intercept")[,1]
b_t_fr=fixef(ib_growth_tarsus_froh, pars = "b_Intercept")[,1]
c_t_fr=fixef(ib_growth_tarsus_froh, pars = "c_Intercept")[,1]

## CIs for intercept
asym1_t_fr_l95=fixef(ib_growth_tarsus_froh, pars = "asym1_Intercept")[,3]
asym1_t_fr_u95=fixef(ib_growth_tarsus_froh, pars = "asym1_Intercept")[,4]
asym_t_fr_l95=fixef(ib_growth_tarsus_froh, pars = "asym_Intercept")[,3]
asym_t_fr_u95=fixef(ib_growth_tarsus_froh, pars = "asym_Intercept")[,4]
b_t_fr_l95=fixef(ib_growth_tarsus_froh, pars = "b_Intercept")[,3]
b_t_fr_u95=fixef(ib_growth_tarsus_froh, pars = "b_Intercept")[,4]
c_t_fr_l95=fixef(ib_growth_tarsus_froh, pars = "c_Intercept")[,3]
c_t_fr_u95=fixef(ib_growth_tarsus_froh, pars = "c_Intercept")[,4]

#inbred with sig growth rate only
fc_t_fr=c_t_fr+(f*fixef(ib_growth_tarsus_froh, pars = "c_FHBD512gen")[,1]) ## c sig
## CIs for inbred
fasym_t_fr_l95=asym_t_fr+(f*(fixef(ib_growth_tarsus_froh, pars = "asym_FHBD512gen")[,3]))
fasym_t_fr_u95=asym_t_fr+(f*(fixef(ib_growth_tarsus_froh, pars = "asym_FHBD512gen")[,4]))
fc_t_fr_l95=c_t_fr+(f*(fixef(ib_growth_tarsus_froh, pars = "c_FHBD512gen")[,3]))
fc_t_fr_u95=c_t_fr+(f*(fixef(ib_growth_tarsus_froh, pars = "c_FHBD512gen")[,4]))



(id_tarsus_froh=ggplot(ib_growth_tarsus_froh$data, aes(x=age_days, y=LeftTarsus*0.1))+
    geom_point(alpha=alpha, colour=dot_col, linewidth=2)+
    theme_classic()+
    theme(text = element_text(size = 18)
          #axis.title.y=element_blank()
          )+
    xlim (0,70)+
    ylim(10,80)+
    labs(x="Age in days", y="Tarsus length (mm)", tag = 'H')+
      stat_function(fun=~ (asym1_t_fr*0.1) + (asym_t_fr*0.1)*exp(-b_t_fr*(c_t_fr)^.x), colour=col_none, size=1, xlim=c(0,80), alpha=0.8)+ ##male (intercept)
      #stat_function(fun=~ (asym1_t_fr*0.1) +(asym_t_fr*0.1)*exp(-b_t_fr*(fc_t_fr)^.x), colour=colinb, size=1, xlim=c(0,80), alpha=0.8)+
        #ave CIs
      stat_function(fun=~ (asym1_t_fr_u95*0.1) +(asym_t_fr_u95*0.1)*exp(-b_t_fr_u95*(c_t_fr_l95)^.x), colour=colave, linetype='dashed', xlim=c(0,80))+
      stat_function(fun=~ (asym1_t_fr_l95*0.1) +(asym_t_fr_l95*0.1)*exp(-b_t_fr_l95*(c_t_fr_u95)^.x), colour=colave,  linetype='dashed', xlim=c(0,80))+
        #inbred CIs
      stat_function(fun=~ (asym1_t_fr*0.1)+(fasym_t_fr_l95*0.1)*exp(-b_t_fr*(fc_t_fr_u95)^.x), colour=colinb,  linetype='dashed', xlim=c(0,80))+
      stat_function(fun=~ (asym1_t_fr*0.1) +(fasym_t_fr_u95*0.1)*exp(-b_t_fr*(fc_t_fr_l95)^.x), colour=colinb,  linetype='dashed', xlim=c(0,80))
  )





# row_label_1 <- wrap_elements(panel = textGrob('First Row', rot=90))
# row_label_2 <- wrap_elements(panel = textGrob('First Row', rot=90))
# 
# 
# (id_bill_froh+id_mass_froh+id_tarsus_froh)/
#   (id_bill_fgrm+id_mass_fgrm+id_tarsus_fgrm)+plot_annotation(tag_levels = 'a')



col_label_1 <- wrap_elements(panel = textGrob(expression('F'[uniWE]),gp = gpar(fontsize = 18)))
col_label_2 <- wrap_elements(panel = textGrob(expression('F'[ROH]),gp = gpar(fontsize = 18)))


plot_all=((col_label_1+col_label_2)/
  (id_bill_fgrm+id_bill_froh)/
  (id_mass_fgrm+id_mass_froh)/
  (id_tarsus_fgrm+id_tarsus_froh)+
  plot_layout( 
              heights = c(0.5, 5, 5, 5)))


plot_all



(bill=(id_bill_fgrm/id_bill_froh))

(mass=(id_mass_fgrm/id_mass_froh))

(tar=(id_tarsus_fgrm/id_tarsus_froh))



plot_all=((col_label_1/col_label_2)|bill|mass|tar)+
  plot_layout(ncol=4, widths = c(1.5, 5, 5, 5))

load(file="Inbreeding_depression_owls/Model_outputs/1_unscaled_ID_wGRM/funiVfroh.RData")

load(file= "Inbreeding_depression_owls/F_dists.RData")

top_plot=both+labs(title='B')

top_left=f_dist+labs(title='A')

layout="AAA#BBB
        AAA#BBB
        AAA#BBB
        #######
        CCCCCCC
        CCCCCCC
        CCCCCCC
        CCCCCCC
        CCCCCCC
        CCCCCCC
        CCCCCCC

"
combined=top_left+top_plot+plot_all+
  plot_layout(design = layout)
combined

ggsave(combined,
       file = "Inbreeding_depression_owls/Model_outputs/beta_and_growth_plots_froh_fgrm2.png",
       width = 12,
       height = 10)




# ggsave(plot_all,
#        file = "Inbreeding_depression_owls/Model_outputs/3_growth/growth_plots_all_froh_fgrm.png",
#        width = 6,
#        height = 7)
# 
