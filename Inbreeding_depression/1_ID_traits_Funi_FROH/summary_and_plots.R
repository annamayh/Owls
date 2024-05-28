
library(tidyverse)
library(brms)
library(viridis)
library(patchwork)
library(ggExtra)

## Basic inbreeding depression models accounting for age ... using 2 types of inbreeding coefficient: FROH and Funi ##
## 3 traits investigated: bill, mass and wing length

setwd("/Users/ahewett1/Documents")
## BILL##
id_bill_funi=readRDS("Inbreeding_depression_owls/Model_outputs/1_unscaled_ID_wGRM/1.1.ID_bill_GRM_Funi_unscaled.RDS")
id_bill_froh=readRDS("Inbreeding_depression_owls/Model_outputs/1_unscaled_ID_wGRM/1.1.ID_bill_GRM_FROH_unscaled.RDS")
## MASS ##
id_mass_funi=readRDS("Inbreeding_depression_owls/Model_outputs/1_unscaled_ID_wGRM/1.2.ID_mass_GRM_Funi_unscaled.RDS")
id_mass_froh=readRDS("Inbreeding_depression_owls/Model_outputs/1_unscaled_ID_wGRM/1.2.ID_mass_GRM_FROH_unscaled.RDS")
## WING ##
id_wing_funi=readRDS("Inbreeding_depression_owls/Model_outputs/1_unscaled_ID_wGRM/1.3.ID_wing_GRM_Funi_unscaled.RDS")
id_wing_froh=readRDS("Inbreeding_depression_owls/Model_outputs/1_unscaled_ID_wGRM/1.3.ID_wing_GRM_FROH_unscaled.RDS")

## check convergence and initial estimates
summary(id_bill_funi) ## 
summary(id_bill_froh)

summary(id_mass_funi) ## running againn 
summary(id_mass_froh)

summary(id_wing_funi)
summary(id_wing_froh)

plot(id_bill_funi) ## 
plot(id_bill_funi)


## function to get model estimates of F and create a df with info of trait and F used 
extract_ID_est=function(model){

  f_ests=fixef(model, pars = paste0(colnames(model$data[2])))%>% ## row 2 of data is the inbreeding coeff used
    as.data.frame()%>%
    add_column(trait=paste0(colnames(model$data[1])))%>% ## row 1 of data is the response variable
    rownames_to_column(var="F")
  f_ests
}


ID_froh_models=list(id_bill_funi,id_bill_froh,id_mass_funi,id_mass_froh,id_wing_funi,id_wing_froh)

f_ests_all_models=NULL
for (i in 1:length(ID_froh_models)){

  ID_model=extract_ID_est(ID_froh_models[[i]])
  f_ests_all_models=rbind(f_ests_all_models,ID_model)
}






funi=f_ests_all_models%>%
  filter(F=="FuniWE")%>%
  ggplot(aes(x=trait, y=Estimate, ymin=Q2.5, ymax=Q97.5,colour=trait)) +
  theme_minimal() +
  geom_pointrange(fatten = 7) +
  coord_flip()+
  geom_hline(yintercept = 0, lty=2)+
  scale_color_viridis(discrete = TRUE, end = 0.8) +
  theme(text = element_text(size = 18), 
        axis.title.y = element_blank(), 
        legend.position = "none")+
  labs(title = "Funi", y="Inbreeding depression (\u03B2)")

funi



froh=f_ests_all_models%>%
  filter(F=="FHBD512gen")%>%
  ggplot(aes(x=trait, y=Estimate, ymin=Q2.5, ymax=Q97.5,colour=trait)) +
  theme_minimal() +
  geom_pointrange(fatten = 7) +
  coord_flip()+
  geom_hline(yintercept = 0, lty=2)+
  scale_color_viridis(discrete = TRUE, end = 0.8) +
  theme(text = element_text(size = 18), 
        axis.title.y = element_blank(), 
        legend.position = "none")+
  labs(title = "FROH", y="Inbreeding depression (\u03B2)")

froh


funi+froh+plot_annotation(tag_levels = 'A')



both=f_ests_all_models%>%
  ggplot(aes(x=trait, y=Estimate, ymin=Q2.5, ymax=Q97.5,colour=trait, shape=F)) +
  theme_minimal() +
  geom_pointrange(fatten = 7, position = position_dodge2(width = 0.5) )+
  coord_flip()+
  geom_hline(yintercept = 0, lty=2)+
  scale_color_viridis(discrete = TRUE, end = 0.8) +
  theme(text = element_text(size = 18), 
        axis.title.y = element_blank())+
  labs( y="Inbreeding depression (\u03B2)")+
  guides(colour='none', shape = guide_legend(reverse = TRUE))+
  scale_shape_discrete(name = "Inbreeding Coeff", labels = c("FROH", "Funi"))

both


ggsave(both, 
       file= "Inbreeding_depression_owls/Model_outputs/unscaled_ID_wGRM/funiVfroh.png",
       width = 7,
       height = 5,
       bg = 'white'
       )
