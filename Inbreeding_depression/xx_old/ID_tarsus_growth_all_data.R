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
  mutate(tarsus_scale=LeftTarsus/100)%>%
  mutate(agedays2=age_days)

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



# ggsave(mass_records_cut, file="Inbreeding_depression_owls/plots/mass_records90_line.png", 
#        width = 8, 
#        height=7)
# 


check= tarsus_df%>%
  filter(age_days<25)
#1661 different ids measured before 25 days old 

#3153 records in 2083 ids before 30 days mark 

n_distinct(check$RingId) 
table(check$year) 
table(check$age_days) 

check%>%
  group_by(age_days)%>%
  summarise(n=n())%>%
ggplot(aes(x=age_days, y=n))+
  geom_point()

check%>%
  group_by(year)%>%
  summarise(n=n())%>%
  ggplot(aes(x=year, y=n))+
  geom_point()



tarsus_df$clutch_merge=as.factor(tarsus_df$clutch_merge)
tarsus_df$GeneticSex=as.factor(tarsus_df$GeneticSex)
tarsus_df$RingId=as.factor(tarsus_df$RingId)
tarsus_df$year=as.factor(tarsus_df$year)
tarsus_df$Observer=as.factor(tarsus_df$Observer)
tarsus_df$nestboxID=as.factor(tarsus_df$nestboxID)
tarsus_df$rank=as.numeric(tarsus_df$rank)
tarsus_df$agedays2=as.numeric(tarsus_df$agedays2)
tarsus_df$age_days=as.numeric(tarsus_df$age_days)
tarsus_df$LeftTarsus=as.numeric(tarsus_df$LeftTarsus)


n_distinct(tarsus_df$RingId)

## standard model with only intercept to get prior

growth.tarsus.inter=brm(
  ## model
  bf(tarsus_scale ~  asym2 * exp(-b*(c)^agedays),
      asym2 + b + c ~ 1 ,
     agedays ~ 0 + me(julian_hatchdate, 2, gr=RingId),
     nl=TRUE),
  data=tarsus_df,
  family = gaussian(),
  control = list(adapt_delta = 0.85),
  init = 0,
  chains = 4


)

 summary(growth.tarsus.inter)

###set priors##
prior_tarsus<- c(
 # prior(normal(1.9, 1), nlpar = "asym1",  class="b"),##
  prior(normal(7, 1), nlpar = "asym2",  class="b"),##
  prior(normal(2.3, 1), nlpar = "b",  class="b"), ## 
  prior(normal(0.9, 0.5), nlpar = "c",  class="b"), ## 
  
  prior(student_t(3, 0, 2.5), class = "sigma", lb=0),
   prior(student_t(3, 0, 2.5), class="sd",nlpar = "asym2", lb=0), # 
 # prior(student_t(3, 0, 2.5), class="sd",nlpar = "asym1", lb=0), # 
  prior(student_t(3, 0, 2.5),  class="sd", nlpar = "b", lb=0),
  prior(student_t(3, 0, 2.5),  class="sd", nlpar = "c", lb=0),

# prior(cauchy(0, 0.5), class="sd", group="RingId", nlpar = "asym1", lb=0), #
 prior(normal(0, 1), class="sd", group="RingId", nlpar = "asym2", lb=0), #
  prior(normal(0, 1),  class="sd", group="RingId", nlpar = "b", lb=0),
  prior(normal(0, 0.5),  class="sd", group="RingId", nlpar = "c", lb=0)

)

growth_tarsus.mod=brm(
  ## model 
  bf(tarsus_scale ~  asym2 * exp(-b*(c)^age_days),
     asym2 ~ 1 + FuniWE + rank + GeneticSex + (1|clutch_merge) + (1|Observer) + (1|year) + (1|RingId),
      b ~ 1 + FuniWE + rank + GeneticSex + (1|clutch_merge) + (1|RingId),
      c ~ 1 + FuniWE + rank + GeneticSex + (1|clutch_merge) + (1|RingId),
     nl=TRUE),
  data=tarsus_df, 
  family = gaussian(),
  chains = 4,
  prior = prior_tarsus,
  control = list(adapt_delta = 0.95),
  init = 0, 
  cores = 4
  # iter = 15000, 
  # warmup = 5000, 
  # thin = 10
  
)

summary(growth_tarsus.mod)
plot(growth_tarsus.mod)


#saveRDS(growth_tarsus.mod,file="Inbreeding_depression_owls/Model_outputs/Growth_tarsus_funi_alldata.RDS")

growth_tarsus.mod=readRDS(file="Inbreeding_depression_owls/Model_outputs/Growth_tarsus_funi_alldata.RDS")

summary(growth_tarsus.mod)
plot(growth_tarsus.mod)

## plot of CIs for F

F_eff=fixef(growth_tarsus.mod, pars = c("asym_FuniWE", "asym_rank", "asym_GeneticSex2"))%>%
  rbind(fixef(growth_tarsus.mod, pars = c("b_FuniWE", "b_rank", "b_GeneticSex2")))%>%
  rbind(fixef(growth_tarsus.mod, pars = c("c_FuniWE", "c_rank", "c_GeneticSex2")))%>%
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

id_tarsus=ggplot(tarsus_df, aes(x=age_days, y=tarsus_scale))+
  geom_point(alpha=0.1, colour="black", size=2)+
  xlim (0,93)+
  ylim(0,8)+
  theme_classic()+
  labs(x="Age in days", y="Tarsus")+
  stat_function(fun=~ 1.86 + 5.35*exp(-7.20*(0.84)^.x), colour=colave, size=1.5, xlim=c(0,88), alpha=0.8)+ ##male (intercept)
  stat_function(fun=~ 7.24*exp(-2.45*(0.87)^.x), colour=colave, size=1.5, xlim=c(0,88), alpha=0.8)+ ## female
  
  stat_function(fun=~1.86 + 5.35*exp(-8.295*(0.84)^.x), colour=colinb, size=1.5, xlim=c(0,88), alpha=0.8)+ #inbred male 
  #stat_function(fun=~3.9*exp(-4.075*(0.89)^.x), colour=colinb, size=1.5, xlim=c(0,88), alpha=0.8)+ #inbred female 
  
  theme(text = element_text(size = 18))+
  annotate("text", x=80, y=7.5, label="Average", colour=colave, size=5)+
  annotate("text", x=80, y=7, label="Inbred", colour=colinb, size=5)
  # annotate("text", x=90, y=4.09, label=female, colour=colave, size=6)+
  # annotate("text", x=90, y=3.80, label=male, colour=colave, size=6)+
  # annotate("text", x=93, y=3.9, label=female, colour=colinb, size=6)+
  # annotate("text", x=93, y=3.61, label=male, colour=colinb, size=6
           

id_tarsus
  
ggsave(id_tarsus, file="Inbreeding_depression_owls/plots/id_tarsus_growth.png", 
       width = 8, 
       height=7)
