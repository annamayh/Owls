
library(tidyverse)
library(brms)
library(viridis)
library(patchwork)
library(ggExtra)


setwd("/Users/ahewett1/Documents")

## models run in brms without including pedigree (froh and Funi)
load(file = "Inbreeding_depression_owls/Model_outputs/brms_funi_froh_simple_all.RData")


funi_models=list(mod_tarsus_all.1.Funi, mod_mass_all.1.Funi, mod_bill_all.1.Funi,mod_wing_all.1.Funi)




f_ests_tar=fixef(mod_tarsus_all.1.Funi, pars = "FuniWE")%>%
  as.data.frame()%>%
  add_column(trait="Tarsus")%>%
  rownames_to_column(var="F")

f_ests_mass=fixef(mod_mass_all.1.Funi, pars = "FuniWE")%>%
  as.data.frame()%>%
  add_column(trait="Mass")%>%
  rownames_to_column(var="F")

f_ests_bill=fixef(mod_bill_all.1.Funi, pars = "FuniWE") %>%
  as.data.frame()%>%
  add_column(trait="Bill")%>%
  rownames_to_column(var="F")

f_ests_wing=fixef(mod_wing_all.1.Funi, pars = "FuniWE") %>%
  as.data.frame()%>%
  add_column(trait="Wing")%>%
  rownames_to_column(var="F")


f_ests_tar_roh=fixef(mod_tarsus_all.1.FROH, pars = "FHBD512gen")%>%
as.data.frame()%>%
  add_column(trait="Tarsus")%>%
  rownames_to_column(var="F")

f_ests_mass_roh=fixef(mod_mass_all.1.FROH, pars = "FHBD512gen")%>%
as.data.frame()%>%
  add_column(trait="Mass")%>%
  rownames_to_column(var="F")

f_ests_bill_roh=fixef(mod_bill_all.1.FROH, pars = "FHBD512gen")%>%
as.data.frame()%>%
  add_column(trait="Bill")%>%
  rownames_to_column(var="F")

f_ests_wing_roh=fixef(mod_wing_all.1.FROH, pars = "FHBD512gen")%>%
as.data.frame()%>%
  add_column(trait="Wing")%>%
  rownames_to_column(var="F")






all_f_ests=f_ests_tar%>%
  rbind(f_ests_mass, f_ests_bill, f_ests_wing)


funi=ggplot(all_f_ests, aes(x=trait, y=Estimate, ymin=Q2.5, ymax=Q97.5,colour=trait)) +
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

ggsave(funi, file="Inbreeding_depression_owls/plots/id_simple_funi.png", 
       width = 4, 
       height=5)




all_froh_ests=f_ests_tar_roh%>%
  rbind(f_ests_mass_roh, f_ests_bill_roh, f_ests_wing_roh)

froh=ggplot(all_froh_ests, aes(x=trait, y=Estimate, ymin=Q2.5, ymax=Q97.5,colour=trait)) +
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

ggsave(froh, file="Inbreeding_depression_owls/plots/id_simple_froh.png", 
       width = 4, 
       height=5)



funi+froh



both=all_froh_ests%>%
  rbind(all_f_ests)


ggplot(both, aes(x=trait, y=Estimate, ymin=Q2.5, ymax=Q97.5,colour=trait, shape=F)) +
  theme_minimal() +
  geom_pointrange(fatten = 7, position = position_dodge2(width = 0.5)) +
  coord_flip()+
  geom_hline(yintercept = 0, lty=2)+
  scale_color_viridis(discrete = TRUE, end = 0.8) +
  theme(text = element_text(size = 18), 
        axis.title.y = element_blank(), 
        )+
  guides(colour = "none")+
  labs( y="Inbreeding depression (\u03B2)")






###############################################################################
### why some sig in froh and some sig in funi ?? ##############################
#############################################################################

## box plots of F


funib=ggplot(tarsus_df, aes(y=FuniWE))+
  geom_boxplot(color="darkblue", fill="lightblue")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


frohb=ggplot(tarsus_df, aes(y=FHBD512gen))+
  geom_boxplot(color="darkblue", fill="lightblue")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

library(patchwork)

funi+froh

funib+frohb


sd(tarsus_df$FuniWE) ## sd higher for Funi
IQR(tarsus_df$FuniWE) ## But IQR higher for Froh

sd(tarsus_df$FHBD512gen)
IQR(tarsus_df$FHBD512gen)

median(tarsus_df$FuniWE)

## IQR measures spread of middle 50% of values
## but SD measures typical devation of individual values from mean so can be more affected by outlires


fvf=tarsus_df%>%
  distinct(RingId, .keep_all = TRUE)%>%
  ggplot(aes(y=FuniWE, x=FHBD512gen))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(text = element_text(size = 18))+
  labs(x="FROH", y="FuniWE")+
  geom_hline(yintercept = mean(tarsus_df$FuniWE), alpha=0.5, size=1)+
  geom_vline(xintercept = mean(tarsus_df$FHBD512gen), alpha=0.5, size=1)


ggMarginal(fvf, fill = "lightblue", size = 3)


