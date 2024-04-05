### id in mass growth using all available data ###

library(tidyverse)
library(brms)
library(corpcor)
library(viridis)

setwd("/Users/ahewett1/Documents/")

##################################################################################
########################### ~~ tarsus ~~ #############################################
#####################################################################################

tarsus_df=read.table("Inbreeding_depression_owls/pheno_df/tarsus_all_pheno_df.txt",sep=",", header=T)%>%
  select(-min_mes)%>%
  mutate(tarsus_scale=LeftTarsus/100)

tarsus_records=ggplot(tarsus_df, aes(x=age_days, y=LeftTarsus))+
  geom_point(alpha=0.2, colour="black", size=2)+
  #xlim (0,1000)+
  ylim(0,800)+
  theme_bw()+
  labs(x="Age in days", y="Tarsus length")+
  theme(text = element_text(size = 18))
tarsus_records

# ggsave(mass_records, file="Inbreeding_depression_owls/plots/mass_records.png", 
#        width = 8, 
#        height=7)
# 


mass_records_cut=ggplot(mass_df, aes(x=age_days, y=Mass))+
  geom_point(alpha=0.2, colour="black", size=2)+
  xlim (0,90)+
  ylim(0,600)+
  theme_bw()+
  labs(x="Age in days", y="Mass (g)")+
  theme(text = element_text(size = 18))+
  stat_function(fun=~ 369.4840*exp(-4.491982*(0.889770)^.x), colour="red", size=1.5, xlim=c(0,88), alpha=0.8)
  

mass_records_cut  

# ggsave(mass_records_cut, file="Inbreeding_depression_owls/plots/mass_records90_line.png", 
#        width = 8, 
#        height=7)
# 


check= tarsus_df%>%
  filter(age_days<30)
#1661 different ids measured before 25 days old 

#3153 records in 2083 ids before 30 days mark 

n_distinct(check$RingId) 

tarsus_df$clutch_merge=as.factor(tarsus_df$clutch_merge)
tarsus_df$GeneticSex=as.factor(tarsus_df$GeneticSex)
tarsus_df$RingId=as.factor(tarsus_df$RingId)
tarsus_df$year=as.factor(tarsus_df$year)
tarsus_df$Observer=as.factor(tarsus_df$Observer)
tarsus_df$rank=as.numeric(tarsus_df$rank)

n_distinct(tarsus_df$RingId)


prior_tarsus<- c(
  prior(normal(3.5, 2), nlpar = "asym",  class="b"),##
  prior(normal(4.5, 1), nlpar = "b",  class="b"), ## 
  prior(normal(1, 0.5), nlpar = "c",  class="b"), ## 
  
  prior(student_t(3, 0, 2.5), class = "sigma", lb=0),
  prior(student_t(3, 0, 2.5), class="sd",nlpar = "asym", lb=0), # 
  prior(student_t(3, 0, 2.5),  class="sd", nlpar = "b", lb=0),
  prior(student_t(3, 0, 2.5),  class="sd", nlpar = "c", lb=0), 
  #prior(student_t(3, 0, 2.5),  class="sd", nlpar = "error", lb=0), 
  
  
  prior(cauchy(0, 1), class="sd", group="RingId", nlpar = "asym", lb=0), #
  prior(cauchy(0, 0.5),  class="sd", group="RingId", nlpar = "b", lb=0),
  prior(cauchy(0, 0.2),  class="sd", group="RingId", nlpar = "c", lb=0)

)

growth_tarsus.mod=brm(
  ## model 
  bf(tarsus_scale ~ asym * exp(-b*(c)^age_days),
     asym ~ 1 + FuniWE + rank + GeneticSex + (1|clutch_merge) + (1|Observer) + (1|year) + (1|RingId),
      b ~ 1 + FuniWE + rank + GeneticSex + (1|clutch_merge) + (1|RingId),
      c ~ 1 + FuniWE + rank + GeneticSex + (1|RingId), 
     nl=TRUE),
  data=tarsus_df, 
  family = gaussian(),
  chains = 4,
  prior = prior_tarsus,
  control = list(adapt_delta = 0.85),
  init = 0, 
  cores = 4,
  iter = 10000, 
  warmup = 4000, 
  thin = 5
  
)

summary(growth_tarsus.mod)
plot(growth_tarsus.mod)


saveRDS(growth_tarsus.mod,file="Inbreeding_depression_owls/Model_outputs/Growth_tarsus_funi_alldata.RDS")

growth_tarsus.mod=readRDS(file="Inbreeding_depression_owls/Model_outputs/Growth_tarsus_funi_alldata.RDS")

## plot of CIs for F

F_eff=fixef(growth_tarsus.mod, pars = c("asym_FuniWE"))%>%
  rbind(fixef(growth_tarsus.mod, pars = c("b_FuniWE")))%>%
  rbind(fixef(growth_tarsus.mod, pars = c("c_FuniWE")))%>%
  as.data.frame()%>%
  rownames_to_column(var="parameter")%>%
  separate_wider_delim(parameter, delim = "_", names = c("param","explanatory"))
  

F_eff

feff=ggplot(F_eff, aes(x=explanatory, y=Estimate, ymin=Q2.5, ymax=Q97.5, colour=param))+
  geom_pointrange(size=1.5, linewidth = 1)+
  geom_hline(yintercept = 0, lty=2)+
  facet_wrap(~param, scales = "free_y")+
  theme_bw() +
  theme(text = element_text(size = 18), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        title = element_text(size = 12), 
        legend.position = "none")+
  scale_colour_brewer(palette = "Dark2")+
  labs(title="Inbreeding depression in tarsus paramters")

feff
 
ggsave(feff, file="Inbreeding_depression_owls/plots/id_tarsus_growth_estimateCIs.png", 
       width = 7, 
       height=5)

  


## plotting predicted slopes   
library(showtext)
female = intToUtf8(9792)
male = intToUtf8(9794)
colinb="goldenrod1"
colave="chocolate"

id_mass=ggplot(mass_df, aes(x=age_days, y=mass_scale))+
  geom_point(alpha=0.1, colour="black", size=2)+
  xlim (0,93)+
  ylim(0,6)+
  theme_classic()+
  labs(x="Age in days", y="Mass (0.01 grams)")+
  stat_function(fun=~ 3.80*exp(-3.68*(0.89)^.x), colour=colave, size=1.5, xlim=c(0,88), alpha=0.8)+ ##male (intercept)
  stat_function(fun=~ 4.09*exp(-3.19*(0.9)^.x), colour=colave, size=1.5, xlim=c(0,88), alpha=0.8)+ ## female
  
  stat_function(fun=~3.61*exp(-4.565*(0.88)^.x), colour=colinb, size=1.5, xlim=c(0,88), alpha=0.8)+ #inbred male 
  stat_function(fun=~3.9*exp(-4.075*(0.89)^.x), colour=colinb, size=1.5, xlim=c(0,88), alpha=0.8)+ #inbred female 
  
  theme(text = element_text(size = 18))+
  annotate("text", x=80, y=4.4, label="Average", colour=colave, size=5)+
  annotate("text", x=80, y=3.3, label="Inbred", colour=colinb, size=5)+
  annotate("text", x=90, y=4.09, label=female, colour=colave, size=6)+
  annotate("text", x=90, y=3.80, label=male, colour=colave, size=6)+
  annotate("text", x=93, y=3.9, label=female, colour=colinb, size=6)+
  annotate("text", x=93, y=3.61, label=male, colour=colinb, size=6)

id_mass
  
ggsave(id_mass, file="Inbreeding_depression_owls/plots/id_mass_growth.png", 
       width = 8, 
       height=7)
