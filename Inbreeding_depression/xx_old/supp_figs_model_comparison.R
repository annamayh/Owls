
setwd("/Users/ahewett1/Documents/Inbreeding_depression_owls")##

library(tidyverse)
library(matrixStats)
library(viridis)
library(MCMCglmm)
library(patchwork)


load(file="Funi_all_model_outputs.RData")##load model outputs for fledgelings
load(file = "FROH_all_model_outputs_Juveniles.RData")

## function to select the ibc from the model outputs in MCMCglmm and find Fhat(ID measure) and CIs
F_hat=function(model_input, ibc){
  F_est=as.data.frame(model_input$Sol)%>%select(starts_with(ibc))
  #ibc_eff=F_est*0.25 
  mean<-apply(F_est,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
  CI_upper<-apply(F_est,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
  CI_lower<-apply(F_est,2,quantile,probs = c(0.025)) #gets lower CI for all solutions
  pred=as.data.frame(matrix(c(mean,CI_upper,CI_lower), nrow = 1, ncol = 3))
  pred
}


all_f=list(mod_tarsus.1.2,mod_tarsus.2.ped, 
              mod_mass.1.2, mod_mass.2.ped, 
              mod_bill.1.2, mod_bill.2.ped,
           mod_tarsus.1.2.FROH,mod_tarsus.2.ped.FROH, 
           mod_mass.1.2.FROH, mod_mass.2.ped.FROH, 
           mod_bill.1.2.FROH, mod_bill.2.ped.FROH)

all_f_hat_juve=NULL
for (i in all_f){
  
  hat=F_hat(i, "F")
  all_f_hat_juve=rbind(all_f_hat_juve,hat)
}



all_f_hat_juve_df=all_f_hat_juve%>%
  rename("est"="V1","upr"="V2","lwr"="V3")%>%
  add_column(mod_name=c("Tarsus.simple.Funi","Tarsus.ped.Funi", 
                        "Mass.simple.Funi", "Mass.ped.Funi", 
                        "Bill.simple.Funi", "Bill.ped.Funi",
                        "Tarsus.simple.Froh","Tarsus.ped.Froh", 
                        "Mass.simple.Froh", "Mass.ped.Froh", 
                        "Bill.simple.Froh", "Bill.ped.Froh"))%>%
  separate_wider_delim(mod_name, delim = ".", names = c("trait","model", "F"), cols_remove = F)%>%
  mutate(mod_name=as.character(mod_name))
  

## just getting the models in the right order
all_f_hat_juve_df$mod_name<-ordered(all_f_hat_juve_df$mod_name,
                          levels = c(
                                     
                                     "Bill.simple.Funi", "Bill.ped.Funi",
                                     "Bill.simple.Froh", "Bill.ped.Froh",
                                     "Mass.simple.Funi", "Mass.ped.Funi",
                                     "Mass.simple.Froh", "Mass.ped.Froh",
                                     "Tarsus.simple.Funi","Tarsus.ped.Funi",
                                     "Tarsus.simple.Froh","Tarsus.ped.Froh",

                                     ))



juve_model_com=ggplot(data=all_f_hat_juve_df, aes(x=mod_name, y=est, ymin=upr, ymax=lwr, col=trait, shape=F)) +
  geom_pointrange(position = position_dodge2(width = 0.5)) +
  scale_colour_brewer(palette = "Dark2")+ ## actually looks alright just black?? to get cols just change fill to colour
  theme_bw() +
  geom_hline(yintercept=0, lty=2) +
  coord_flip()+
  theme(text = element_text(size = 12),
              axis.title.y = element_blank())+
  labs(y="Inbreeding depression (\u03B2)", shape="Inbreeding coefficient", col="Trait", 
       x="Model used", title = "Estimates in juveniles")

juve_model_com



ggsave(juve_model_com,file="beta_juveniles_models.png",
       height = 5,
       width = 8)

summary(mod_tarsus.1.2.FROH)
summary(mod_tarsus.1.2)
## only thing that changes slightly is the effect of the sex??

plot(mod_bill.1.2.FROH)
plot(mod_bill.1.2)

### Now Elo's adult models 



load(file="FuniWE_SimpleIDModels_ADULTS.RData")##load model outputs for fledgelings
load(file = "FHBD_SimpleIDModels_ADULTS (1).RData")

summary(EggsHatchingIDMprobitLINK.FuniWE)

all_f_adults=list(TarsuslengthIDM.FuniWE, TarsuslengthIDM.FHBD,
                  MassIDM.FuniWE, MassIDM.FHBD,
                  BilllengthIDM.FuniWE, BilllengthIDM.FHBD)



all_f_hat_adults=NULL
for (i in all_f_adults){
  
  hat=F_hat(i, "F")
  all_f_hat_adults=rbind(all_f_hat_adults,hat)
}


all_f_hat_adults_df=all_f_hat_adults%>%
  rename("est"="V1","upr"="V2","lwr"="V3")%>%
  add_column(mod_name=c("Tarsus.simple.Funi","Tarsus.simple.Froh", 
                        "Mass.simple.Funi", "Mass.simple.Froh", 
                        "Bill.simple.Funi", "Bill.simple.Froh"
                        ))%>%
  separate_wider_delim(mod_name, delim = ".", names = c("trait","model", "F"), cols_remove = F)%>%
  mutate(mod_name=as.character(mod_name))



adult_model_com=ggplot(data=all_f_hat_adults_df, aes(x=mod_name, y=est, ymin=upr, ymax=lwr, col=trait, shape=F)) +
  geom_pointrange(position = position_dodge2(width = 0.5)) +
  scale_colour_brewer(palette = "Dark2")+ ## actually looks alright just black?? to get cols just change fill to colour
  theme_bw() +
  geom_hline(yintercept=0, lty=2) +
  coord_flip()+
  theme(text = element_text(size = 12),
              axis.title.y = element_blank())+
  labs(y="Inbreeding depression (\u03B2)", shape="Inbreeding coefficient", col="Trait", 
       x="Model used", title = "Estimates in adults")

adult_model_com


# 
# juve_model_com/adult_model_com+plot_layout(guides = "collect")+
#   plot_annotation(tag_levels = 'A')
# 
# ggsave(juve_model_com,file="beta_juveniles_models.png",
#        height = 5,
#        width = 8)

## also eggs 

summary(NBEggslaidIDM.FuniWE)
summary(EggsHatchingIDMprobitLINK.FHBD)

## have to modify function slightly as interaction means it picks up interaction too
F_hat_2=function(model_input, ibc){
  F_est=as.data.frame(model_input$Sol)%>%select(all_of(ibc))
  #ibc_eff=F_est*0.25 
  mean<-apply(F_est,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
  CI_upper<-apply(F_est,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
  CI_lower<-apply(F_est,2,quantile,probs = c(0.025)) #gets lower CI for all solutions
  pred=as.data.frame(matrix(c(mean,CI_upper,CI_lower), nrow = 1, ncol = 3))
  pred
}


hat_eggslaid_Funi=F_hat_2(NBEggslaidIDM.FuniWE, "FuniWE")
hat_eggslaid_Froh=F_hat_2(NBEggslaidIDM.FHBD, "FHBD")
hat_eggs_hatching_mum_Funi=F_hat_2(EggsHatchingIDMprobitLINK.FuniWE, "MotherFuniWE")
hat_eggs_hatching_dad_Funi=F_hat_2(EggsHatchingIDMprobitLINK.FuniWE, "FatherFuniWE")
hat_eggs_hatching_mum_Froh=F_hat_2(EggsHatchingIDMprobitLINK.FHBD, "MotherFHBD")
hat_eggs_hatching_dad_Froh=F_hat_2(EggsHatchingIDMprobitLINK.FHBD, "FatherFHBD")






eggies=rbind(hat_eggslaid_Funi,hat_eggslaid_Froh)%>%
  rbind(hat_eggs_hatching_mum_Funi)%>%
  rbind(hat_eggs_hatching_dad_Funi)%>%
  rbind(hat_eggs_hatching_mum_Froh)%>%
  rbind(hat_eggs_hatching_dad_Froh)


fhat_eggies=eggies%>%
  rename("est"="V1","upr"="V2","lwr"="V3")%>%
  add_column(mod_name=c(  "Eggslaid.simple.Funi", "Eggslaid.simple.Froh",
                          "Eggshatched\nMumibc.simple.Funi", "Eggshatched\nDadibc.simple.Funi",
                          "Eggshatched\nMumibc.simple.Froh", "Eggshatched\nDadibc.simple.Froh"
                          
  ))%>%
  separate_wider_delim(mod_name, delim = ".", names = c("trait","model", "F"), cols_remove = F)%>%
  mutate(mod_name=as.character(mod_name))

  
fhat_eggies

adult_model_com_eggies=ggplot(data=fhat_eggies, aes(x=mod_name, y=est, ymin=upr, ymax=lwr, col=trait, shape=F)) +
  geom_pointrange(position = position_dodge2(width = 0.5)) +
  scale_colour_brewer(palette = "Accent")+ ## actually looks alright just black?? to get cols just change fill to colour
  theme_bw() +
  geom_hline(yintercept=0, lty=2) +
  coord_flip()+
  theme(text = element_text(size = 12),
        axis.title.y = element_blank())+
  labs(y="Inbreeding depression (\u03B2)", shape="Inbreeding coefficient", col="Trait", 
       x="Model used", title = "Estimates in adults (egg traits)")

adult_model_com_eggies

all_comparison=juve_model_com+adult_model_com/adult_model_com_eggies+
  plot_layout(guides = "collect")+
  plot_annotation(tag_levels = 'A')

all_comparison

ggsave(all_comparison,file="beta_all_models.png",
       height = 6,
       width = 11)


Funi=read.table("All3085_FuniWE_FHBD512g.txt", header = T)%>%
  rename(RingId=INDVs)



plot(fledge_tarsus_df$FuniWE, fledge_tarsus_df$FHBD512gen)
