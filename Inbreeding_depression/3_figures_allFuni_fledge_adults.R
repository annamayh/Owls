

setwd("/Users/ahewett1/Documents/Inbreeding_depression_owls")##

library(tidyverse)
library(matrixStats)
library(viridis)
library(MCMCglmm)
library(patchwork)


load(file="Funi_all_model_outputs.RData")##load model outputs for fledgelings
load(file="FuniWE_SimpleIDModels_ADULTS.RData") ## Elos models using Funi (simple)
load(file="noF_AnimalModels_ADULTS.RData") ## Elos models including ped


#check random effects and their order 
posterior.mode(mod_tarsus.2.ped$VCV)
posterior.mode(TarsuslengthAM$VCV)

summary(TarsuslengthAM)


posterior.mode(mod_mass.2.ped$VCV)
as.data.frame(posterior.mode(MassAM$VCV))

summary(MassAM)


posterior.mode(mod_bill.1.2$VCV)
posterior.mode(BilllengthAM$VCV)



## this function only works because all of the random effects are in the same order and are all the same 
## will need to edit if other things are included / they are out of order 
get_var_partitions=function(model_input){

      tar_randoms=as.data.frame(posterior.mode(model_input$VCV))
      
      total_V_tar=as.numeric(colSums(tar_randoms))
      va=as.numeric(tar_randoms[1,]/total_V_tar)
      
      pe=as.numeric(tar_randoms[2,]/total_V_tar)
      obs=as.numeric(tar_randoms[3,]/total_V_tar)
      clutch=as.numeric(tar_randoms[4,]/total_V_tar)
      year=as.numeric(tar_randoms[5,]/total_V_tar)
      res=as.numeric(tar_randoms[6,]/total_V_tar)
      
      
      Var_ped=rbind(va,pe)%>%
        rbind(obs)%>%
        rbind(clutch)%>%
        rbind(year)%>%
        rbind(res)%>%
        as.data.frame()%>%
        tibble::rownames_to_column("Variance_element")
      
      Var_ped

}



##tarsus

flede_tarsus=get_var_partitions(mod_tarsus.2.ped)%>%
  mutate(age="Juvenile")%>%
  mutate(trait="Tarsus Length")
adult_tarsus=get_var_partitions(TarsuslengthAM)%>%
  mutate(age="Adult")%>%
  mutate(trait="Tarsus Length")

#bill length

flede_bill=get_var_partitions(mod_bill.2.ped)%>%
  mutate(age="Juvenile")%>%
  mutate(trait="Bill Length")

adult_bill=get_var_partitions(BilllengthAM)%>%
  mutate(age="Adult")%>%
  mutate(trait="Bill Length")

#Mass

flede_mass=get_var_partitions(mod_mass.2.ped)%>%
  mutate(age="Juvenile")%>%
  mutate(trait="Mass")

adult_mass=get_var_partitions(MassAM)%>%
  mutate(age="Adult")%>%
  mutate(trait="Mass")

all_var_comp=rbind(flede_tarsus, adult_tarsus)%>%
  rbind(flede_bill)%>%rbind(adult_bill)%>%
  rbind(flede_mass)%>%rbind(adult_mass)%>%
  mutate()

all_var_comp$Variance_element<-as.factor(all_var_comp$Variance_element)
all_var_comp$Variance_element<-ordered(all_var_comp$Variance_element, 
                            levels = c("year", "pe", "clutch", 
                                       "obs","res","va"))



var_plot_lengths=ggplot(all_var_comp, aes(fill=Variance_element, x=factor(age, level = c('Juvenile', 'Adult')), y=V1)) +
  geom_bar(position="stack",stat="identity")+
  theme_bw() +
  facet_wrap(~trait)+
  scale_fill_viridis_d(direction = -1, 
                       labels=c(bquote (V["Year born"]), 
                                bquote (V["pe"]),
                                bquote (V["Clutch"]),
                                bquote (V["Observer"]),
                                bquote (V["Residual"]),
                                bquote (V["a"])
                                ))+
  theme(text = element_text(size = 15), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 15)
        )+
  labs(fill="Variance \ncomponents")


var_plot_lengths


# ggsave(var_plot_lengths,file="variance_plots.png",
#        height = 5,
#        width = 9)
# 

######## TO GO IN SUPP MAT #####

## now adult only traits
# 
# posterior.mode(NBEggslaidAM$VCV)
# summary(NBEggslaidAM)
# 
# eggs_laid_randoms=as.data.frame(posterior.mode(NBEggslaidAM$VCV))
# 
# total_V_tar=as.numeric(colSums(eggs_laid_randoms))
# va=as.numeric(eggs_laid_randoms[1,]/total_V_tar)
# 
# pe=as.numeric(eggs_laid_randoms[2,]/total_V_tar)
# year=as.numeric(eggs_laid_randoms[3,]/total_V_tar)
# site=as.numeric(eggs_laid_randoms[4,]/total_V_tar)
# res=as.numeric(eggs_laid_randoms[5,]/total_V_tar)
# 
# 
# Var_ped_eggs_laid=rbind(va,pe)%>%
#   rbind(site)%>%
#   rbind(year)%>%
#   rbind(res)%>%
#   as.data.frame()%>%
#   tibble::rownames_to_column("Variance_element")
# 
# Var_ped_eggs_laid
# 
# Var_ped_eggs_laid$Variance_element<-as.factor(Var_ped_eggs_laid$Variance_element)
# Var_ped_eggs_laid$Variance_element<-ordered(Var_ped_eggs_laid$Variance_element, 
#                                        levels = c("year", "pe", "site", 
#                                                   "res","va"))
# 
# 
# var_plot_lengths_eggs_laid=ggplot(Var_ped_eggs_laid, aes(fill=Variance_element, y=V1, x=1)) +
#   geom_bar(position="stack",stat="identity")+
#   theme_bw() +
#   scale_fill_viridis_d(direction = -1, 
#                        labels=c(bquote (V["Year of clutch"]), 
#                                 bquote (V["pe"]),
#                                 bquote (V["Site"]),
#                                 bquote (V["Residual"]),
#                                 bquote (V["a"])
#                        ))+
#   theme(text = element_text(size = 15), 
#         axis.title.y = element_blank(),
#         axis.title.x = element_blank(),
#         legend.text = element_text(size = 15)
#   )+
#   labs(fill="Variance \ncomponents")
# 
# var_plot_lengths_eggs_laid
# 
# 


####################################################################################################
####################################################################################################
###   BETA ESTIMATE PLOTS #######
####################################################################################################

## ONLY USING 'SIMPLE' MODEL for 

summary(mod_tarsus.1.2)
summary(TarsuslengthIDM.FuniWE)

summary(mod_mass.1.2)
summary(MassIDM.FuniWE)

summary(mod_bill.1.2)
summary(BilllengthIDM.FuniWE)


estimate_F_effect=function(model_input){
    F_est=as.data.frame(model_input$Sol)%>%dplyr::select(matches("FuniWE"))
    
    #ibc_eff=F_est*-0.25 
    mean<-apply(F_est,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
    CI_upper<-apply(F_est,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
    CI_lower<-apply(F_est,2,quantile,probs = c(0.025)) #gets lower CI for all solutions
    
    pred=as.data.frame(matrix(c(mean,CI_upper,CI_lower), nrow = 1, ncol = 3))
    pred
}


flede_tarsus_est=estimate_F_effect(mod_tarsus.1.2)%>%
  mutate(age="Juvenile")%>%
  mutate(trait="Tarsus \nLength")
adult_tarsus_est=estimate_F_effect(TarsuslengthIDM.FuniWE)%>%
  mutate(age="Adult")%>%
  mutate(trait="Tarsus \nLength")

flede_mass_est=estimate_F_effect(mod_mass.1.2)%>%
  mutate(age="Juvenile")%>%
  mutate(trait="Mass")
adult_mass_est=estimate_F_effect(MassIDM.FuniWE)%>%
  mutate(age="Adult")%>%
  mutate(trait="Mass")

flede_bill_est=estimate_F_effect(mod_bill.1.2)%>%
  mutate(age="Juvenile")%>%
  mutate(trait="Bill \nLength")
adult_bill_est=estimate_F_effect(BilllengthIDM.FuniWE)%>%
  mutate(age="Adult")%>%
  mutate(trait="Bill \nLength")

all_ests_F=rbind(flede_tarsus_est, adult_tarsus_est)%>%
  rbind(flede_mass_est)%>%rbind(adult_mass_est)%>%
  rbind(flede_bill_est)%>%rbind(adult_bill_est)

all_ests_F$trait<-ordered(all_ests_F$trait, 
                       levels = c("Tarsus \nLength", "Mass", "Bill \nLength"))

beta_ests_Funi=ggplot(data=all_ests_F, aes(x=trait, y=V1, ymin=V2, ymax=V3,shape=age)) +
  geom_pointrange(position = position_dodge2(width = 0.5)) +
  scale_colour_brewer(palette = "Dark2")+ ## actually looks alright just black?? to get cols just change fill to colour
  theme_bw() +
  geom_hline(yintercept=0, lty=2) +
  coord_flip()+
  theme(title = element_text(size = 12),
        text = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_blank(),
        legend.text =  element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "left")+
  labs(y="Inbreeding depression (\u03B2)", shape="Life stage", title = "Adult and juvenile traits")+
  guides(colour = "none")+
  geom_text(aes(x=3,y=-5 ,label="ns"),colour="black",size=3)+
  geom_text(aes(x=3.35,y=-17 ,label="*** \n(p<0.01)"),colour="black", size=3)+
  geom_text(aes(x=2,y=18 ,label="ns"),colour="black",size=3)+
  geom_text(aes(x=2.35,y=-28 ,label=". \n(p=0.08)"),colour="black", size=3)+
  geom_text(aes(x=1,y=31 ,label="ns"),colour="black",size=3)+
  geom_text(aes(x=1.35,y=-28 ,label=". \n(p=0.07)"),colour="black" ,size=3)

beta_ests_Funi


# ggsave(beta_ests_Funi,file="Owls_temp/beta_Funi_traits.png",
#        height = 5,
#        width = 8)


  

##and for eggs 

summary(NBEggslaidIDM.FuniWE)
number_eggs=estimate_F_effect(NBEggslaidIDM.FuniWE)%>%
  mutate(age="# eggs laid")%>%
  mutate(trait="# eggs laid")


summary(EggsHatchingIDMprobitLINK.FuniWE)
F_est_mother=as.data.frame(EggsHatchingIDMprobitLINK.FuniWE$Sol)%>%dplyr::select(matches("MotherFuniWE"))
#ibc_eff=F_est*-0.25 
mean<-apply(F_est_mother,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
CI_upper<-apply(F_est_mother,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
CI_lower<-apply(F_est_mother,2,quantile,probs = c(0.025)) #gets lower CI for all solutions
pred_mother=as.data.frame(matrix(c(mean,CI_upper,CI_lower), nrow = 1, ncol = 3))%>%
  add_column(age="Probability of \nEggs hatching \n(Female)")%>%
  add_column(trait="Probability of \nEggs hatching")


F_est_father=as.data.frame(EggsHatchingIDMprobitLINK.FuniWE$Sol)%>%dplyr::select(matches("FatherFuniWE"))
#ibc_eff=F_est*-0.25 
mean<-apply(F_est_father,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
CI_upper<-apply(F_est_father,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
CI_lower<-apply(F_est_father,2,quantile,probs = c(0.025)) #gets lower CI for all solutions
pred_father=as.data.frame(matrix(c(mean,CI_upper,CI_lower), nrow = 1, ncol = 3))%>%
  add_column(age="Probability of \nEggs hatching \n(Male)")%>%
  add_column(trait="Probability of \nEggs hatching")

eggs=pred_father%>%
  rbind(pred_mother)%>%
  rbind(number_eggs)
eggs

my_colors <- RColorBrewer::brewer.pal(2, "Dark2")[4:5]

beta_ests_Funi_eggs=ggplot(data=eggs, aes(x=age, y=V1, ymin=V2, ymax=V3)) +
  geom_pointrange(position = position_dodge2(width = 0.5)) +
 # scale_colour_brewer(palette = "Greens", direction = -1)+ ## actually looks alright just black?? to get cols just change fill to colour
  theme_bw() +
  geom_hline(yintercept=0, lty=2) +
  coord_flip()+
  theme(title = element_text(size = 12),
        text = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_blank())+
  labs(y="Inbreeding depression (\u03B2)", title = "Adult only traits")+
  guides(color = "none")+
  geom_text(aes(x=1.15,y=0.9 ,label="ns"),colour="black",size=3)+
  geom_text(aes(x=2.15,y=1.8 ,label="ns"),colour="black",size=3)+
  geom_text(aes(x=3.15,y=4 ,label="ns"),colour="black",size=3)

beta_ests_Funi_eggs




all_beta_traits=beta_ests_Funi+beta_ests_Funi_eggs+ plot_layout(guides = "collect")+
  plot_annotation(tag_levels = 'A')

all_beta_traits


ggsave(all_beta_traits,file="Owls_temp/beta_Funi_traits_all.png",
       height = 5,
       width = 10)


