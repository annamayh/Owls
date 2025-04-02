
library(tidyverse)
library(brms)
library(viridis)
library(patchwork)
library(ggExtra)

## Basic inbreeding depression models accounting for age ... using 2 types of inbreeding coefficient: FROH and Funi ##
## 3 traits investigated: bill, mass and tarsus length

setwd("/Users/ahewett1/Documents")
## BILL##
id_bill_funi=readRDS("Inbreeding_depression_owls/Model_outputs/1_unscaled_ID_wGRM/1.1.ID_bill_GRM_Funi_unscaled.RDS")
id_bill_froh=readRDS("Inbreeding_depression_owls/Model_outputs/1_unscaled_ID_wGRM/1.1.ID_bill_GRM_FROH_unscaled.RDS")
## MASS ##
id_mass_funi=readRDS("Inbreeding_depression_owls/Model_outputs/1_unscaled_ID_wGRM/1.2.ID_mass_GRM_Funi_unscaled.RDS")
id_mass_froh=readRDS("Inbreeding_depression_owls/Model_outputs/1_unscaled_ID_wGRM/1.2.ID_mass_GRM_FROH_unscaled.RDS")
## tarsus ##
id_tarsus_funi=readRDS("Inbreeding_depression_owls/Model_outputs/1_unscaled_ID_wGRM/1.4.ID_tarsus_GRM_Funi.RDS")
id_tarsus_froh=readRDS("Inbreeding_depression_owls/Model_outputs/1_unscaled_ID_wGRM/1.4.ID_tarsus_GRM_FROH.RDS")

#bill_data=read.table("Inbreeding_depression_owls/pheno_df/bill_all_pheno_df.txt",sep=",", header=T)

## check convergence and initial estimates
summary(id_bill_funi) ## 
summary(id_bill_froh)

summary(id_mass_funi) ##
summary(id_mass_froh)

summary(id_tarsus_funi)
summary(id_tarsus_froh)

0.25*-43.41


#plot(id_tarsus_froh)


#plot(id_bill_funi) ## 
#plot(id_bill_funi)


## function to get model estimates of F and create a df with info of trait and F used 
extract_ID_est=function(model){

  f_ests=fixef(model, pars = paste0(colnames(model$data[2])))%>% ## row 2 of data is the inbreeding coeff used
    as.data.frame()%>%
    add_column(trait=paste0(colnames(model$data[1])))%>% ## row 1 of data is the response variable
    rownames_to_column(var="F")
  f_ests
}


ID_froh_models=list(id_bill_funi,id_bill_froh,id_mass_funi,id_mass_froh,id_tarsus_funi,id_tarsus_froh)

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
  mutate(trait = factor(trait, levels=c('LeftTarsus','Mass', 'BillLength'))) %>%
  ggplot(aes(x=trait, y=Estimate, ymin=Q2.5, ymax=Q97.5,colour=trait, shape=F)) +
  theme_minimal() +
  geom_pointrange(fatten = 7, position = position_dodge2(width = 0.5) )+
  coord_flip()+
  geom_hline(yintercept = 0, lty=2)+
  scale_color_viridis(discrete = TRUE, end = 0.8) +
  theme(text = element_text(size = 18), 
        axis.title.y = element_blank(), 
        legend.title = element_blank())+
  labs( y="Inbreeding depression (\u03B2)")+
  guides(colour='none', shape = guide_legend(reverse = TRUE))+
  scale_shape_discrete(labels = c(expression('F'[ROH]), expression('F'[uniWE])))+
  scale_x_discrete(labels = c('BillLength' = 'Bill \nLength', 'Mass' = 'Mass', 'LeftTarsus' = 'Tarsus \nLength'))
both


ggsave(both, 
       file= "Inbreeding_depression_owls/Model_outputs/1_unscaled_ID_wGRM/funiVfroh.png",
       width = 5,
       height = 4,
       bg = 'white'
       )



save(both, 
       file= "Inbreeding_depression_owls/Model_outputs/1_unscaled_ID_wGRM/funiVfroh.RData")


mean(id_bill_froh$data$BillLength)

mean(id_mass_froh$data$Mass)

