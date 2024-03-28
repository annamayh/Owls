### id in mass growth using all available data ###

library(tidyverse)
library(brms)
library(corpcor)

##################################################################################
########################### ~~ MASS ~~ #############################################
#####################################################################################

mass_df=read.table("Inbreeding_depression_owls/pheno_df/mass_all_pheno_df.txt",sep=",", header=T)%>%
  select(-min_mes)%>%
  mutate(mass_scale=Mass/100)%>%
  filter(mass_scale<10) ## remove one id that obviously has wrong mass record 4201g (most likely 420g i hope)

mass_records=ggplot(mass_df, aes(x=age_days, y=Mass))+
  geom_point(alpha=0.2, colour="black", size=2)+
  #xlim (0,1000)+
  ylim(0,600)+
  theme_bw()+
  labs(x="Age in days", y="Mass (g)")+
  theme(text = element_text(size = 18))

ggsave(mass_records, file="Inbreeding_depression_owls/plots/mass_records.png", 
       width = 8, 
       height=7)



mass_records_cut=ggplot(mass_df, aes(x=age_days, y=Mass))+
  geom_point(alpha=0.2, colour="black", size=2)+
  xlim (0,90)+
  ylim(0,600)+
  theme_bw()+
  labs(x="Age in days", y="Mass (g)")+
  theme(text = element_text(size = 18))+
  stat_function(fun=~ 369.4840*exp(-4.491982*(0.889770)^.x), colour="red", size=1.5, xlim=c(0,88), alpha=0.8)
  

mass_records_cut  

ggsave(mass_records_cut, file="Inbreeding_depression_owls/plots/mass_records90_line.png", 
       width = 8, 
       height=7)



check= mass_df%>%
  filter(age_days<30)
#1661 different ids measured before 25 days old 

#3153 records in 2083 ids before 30 days mark 

n_distinct(check$RingId) 

mass_df$clutch_merge=as.factor(mass_df$clutch_merge)
mass_df$GeneticSex=as.factor(mass_df$GeneticSex)
mass_df$RingId=as.factor(mass_df$RingId)
mass_df$year=as.factor(mass_df$year)
mass_df$Observer=as.factor(mass_df$Observer)
mass_df$rank=as.numeric(mass_df$rank)

n_distinct(mass_df$RingId)


prior_mass<- c(
  prior(normal(3.5, 2), nlpar = "asym",  class="b"),##
  prior(normal(4.5, 1), nlpar = "b",  class="b"), ## 
  prior(normal(1, 0.5), nlpar = "c",  class="b"), ## 
  
  prior(student_t(3, 0, 2.5), class = "sigma", lb=0),
  prior(student_t(3, 0, 2.5), class="sd",nlpar = "asym", lb=0), # 
  prior(student_t(3, 0, 2.5),  class="sd", nlpar = "b", lb=0)
  #prior(student_t(3, 0, 2.5),  class="sd", nlpar = "c", lb=0)
  
)

growth_mass.mod=brm(
  ## model 
  bf(mass_scale ~ asym * exp(-b*(c)^age_days),
     asym ~ 1 + FuniWE + rank + GeneticSex +(1|clutch_merge) + (1|Observer)+(1|year),
     b ~ 1 + FuniWE + rank + GeneticSex + (1|clutch_merge),
     c ~ 1 + FuniWE + rank + GeneticSex, 
     nl=TRUE),
  data=mass_df, 
  family = gaussian(),
  chains = 4,
  prior = prior_mass,
  control = list(adapt_delta = 0.85),
  init = 0, 
  cores = 4,
   iter = 10000, 
   warmup = 3000
  
)

summary(growth_mass.mod)
plot(growth_mass.mod)


saveRDS(growth_mass.mod,file="Inbreeding_depression_owls/Model_outputs/Growth_mass_funi_alldata.RDS")

growth_mass.mod=readRDS(file="Inbreeding_depression_owls/Model_outputs/Growth_mass_funi_alldata.RDS")

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


###############################################################################
## same model using FROH #################################################
##########################################################################

growth_mass.mod.froh=brm(
  ## model 
  bf(mass_scale ~ asym * exp(-b*(c)^age_days),
     asym ~ 1 + FHBD512gen + rank + GeneticSex +(1|clutch_merge) + (1|Observer)+(1|year),
     b ~ 1 + FHBD512gen + rank + GeneticSex + (1|clutch_merge),
     c ~ 1 + FHBD512gen + GeneticSex, 
     nl=TRUE),
  data=mass_df, 
  family = gaussian(),
  chains = 4,
  prior = prior_mass,
  control = list(adapt_delta = 0.85),
  init = 0, 
  cores = 4,
  iter = 10000, 
  warmup = 3000
  
)

summary(growth_mass.mod.froh)
plot(growth_mass.mod.froh)


saveRDS(growth_mass.mod.froh,file="Inbreeding_depression_owls/Model_outputs/Growth_mass_FROH_alldata.RDS")













####

prior_mass_rep<- c(
  prior(normal(3.5, 2), nlpar = "asym",  class="b"),##
  prior(normal(4.5, 1), nlpar = "b",  class="b"), ## 
  prior(normal(1, 0.5), nlpar = "c",  class="b"), ## 
  
  prior(student_t(3, 0, 2.5), class = "sigma", lb=0),
  prior(student_t(3, 0, 2.5), class="sd",nlpar = "asym", lb=0), # 
  prior(student_t(3, 0, 2.5),  class="sd", nlpar = "b", lb=0),
  prior(student_t(3, 0, 2.5),  class="sd", nlpar = "c", lb=0),

  prior(cauchy(0, 0.2), class="sd", group="RingId", nlpar = "asym", lb=0), #
  prior(cauchy(0, 0.2),  class="sd", group="RingId", nlpar = "b", lb=0),
  prior(cauchy(0, 0.2),  class="sd", group="RingId", nlpar = "c", lb=0)

)
growth_mass.mod.repm=brm(
  ## model 
  bf(mass_scale ~ asym * exp(-b*(c)^age_days),
     asym ~ 1 + FuniWE + rank + GeneticSex +(1|RingId) + (1|clutch_merge) + (1|Observer)+(1|year),
     b ~ 1 + FuniWE + rank + GeneticSex + (1|RingId) + (1|clutch_merge),
     c ~ 1 + FuniWE + GeneticSex + (1|RingId), 
     nl=TRUE),
  data=mass_df, 
  family = gaussian(),
  chains = 4,
  prior = prior_mass,
  control = list(adapt_delta = 0.85),
  init = 0, 
  cores = 4,
   iter = 10000, 
   warmup = 3000
  
)

summary(growth_mass.mod.repm)
plot(growth_mass.mod.repm)

saveRDS(growth_mass.mod.repm,file="Inbreeding_depression_owls/Model_outputs/Growth_mass_funi_alldata_ringid.RDS")


### for ids one year old 

mass_df_1yr=mass_df%>%
  filter(age_days<600)

plot(mass_df_1yr$age_days, mass_df_1yr$mass_scale)

plot(mass_df$age_days, mass_df$mass_scale)


prior_mass<- c(
  prior(normal(3.5, 2), nlpar = "asym",  class="b"),##
  prior(normal(4.5, 1), nlpar = "b",  class="b"), ## 
  prior(normal(1, 0.5), nlpar = "c",  class="b"), ## 
  
  prior(student_t(3, 0, 2.5), class = "sigma", lb=0),
  prior(student_t(3, 0, 2.5), class="sd",nlpar = "asym", lb=0), # 
  prior(student_t(3, 0, 2.5),  class="sd", nlpar = "b", lb=0)
  #prior(student_t(3, 0, 2.5),  class="sd", nlpar = "c", lb=0
  # ) #,
  
)

growth_mass.mod=brm(
  ## model 
  bf(mass_scale ~ asym * exp(-b*(c)^age_days),
     asym ~ 1 + FuniWE + rank + GeneticSex +(1|clutch_merge) + (1|Observer)+(1|year),
     b ~ 1 + FuniWE + rank + GeneticSex + (1|clutch_merge),
     c ~ 1 + FuniWE + GeneticSex, 
     nl=TRUE),
  data=mass_df, 
  family = gaussian(),
  chains = 4,
  prior = prior_mass,
  control = list(adapt_delta = 0.85),
  init = 0, 
  cores = 4,
  iter = 15000, 
  warmup = 5000
  
)