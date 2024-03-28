## df for modelling inbreeding depression in growth using estimated hatch date##

## https://encyclopedia.pub/entry/626 

library(tidyverse)
library(brms)

setwd("/Users/ahewett1/Documents/")

##Fgrm and FROH from elo
inb=read.table("Inbreeding_depression_owls/All3085_FuniWE_FHBD512g.txt", header = T)%>%
  rename(RingId=INDVs)
## reading in phenotype info
owl_mes=read.csv("Barn_owl_general/BarnOwls_Legacy_20231010153920/BarnOwls_Legacy_20231010153920_BirdMeasurement.csv",header = T)%>%
select(RingId,LeftTarsus, Mass, LeftWing, BillLength, Observer,ObservationDate,EstimatedGrowthStage)%>% ##only selecting the paramters we need
  separate_wider_delim(ObservationDate, delim = " ", names = c("obs_date",NA))%>%
  mutate(obs_date=as.Date(obs_date, format = "%d/%m/%Y"))%>%
  filter(EstimatedGrowthStage=="Fledgling"|EstimatedGrowthStage=="Nestling or Fledgling")%>% ##only taking nestlings/fledgelings i.e. first year of life 
  unique()


## genetic sex
owl_sex=read.table("Inbreeding_depression_owls/GeneticSex_3K_RP548.txt", header = T) # 
## info on clutch 
owl_bird=read.csv("Barn_owl_general/BarnOwls_Legacy_20231010153920/BarnOwls_Legacy_20231010153920_Bird.csv", header = T)%>%
  select(RingId, BornClutchId, RaisedClutchId, SiteId, HatchDate)%>%
  na.omit()%>%
  separate_wider_delim(HatchDate, delim = " ", names = c("hatch_date",NA))%>%
  mutate(hatch_date=as.Date(hatch_date, format = "%d/%m/%Y"))%>%
  unite("clutch_merge", BornClutchId:RaisedClutchId, remove = F)# merging raised and born clutch to one variable 

fledge_merge=inb%>%
  inner_join(owl_mes)%>%
  inner_join(owl_sex)%>%
  inner_join(owl_bird,relationship = "many-to-many")%>%
  separate_wider_delim(hatch_date, delim = "-", cols_remove=F,names = c("year",NA,NA))

fledge_merge$date_diff <- as.Date(as.character(fledge_merge$obs_date))-
  as.Date(as.character(fledge_merge$hatch_date)) ## getting rough age of id when observation was taken based on hatch date

## also getting info on rank i.e. what order they were born in in the nest
ranked=fledge_merge%>%
  select(RingId,hatch_date,RaisedClutchId)%>%
  distinct(RingId, .keep_all = TRUE)%>%
  group_by(RaisedClutchId)%>%
  mutate(rank = rank(hatch_date, ties.method= "min")) %>%
  arrange(hatch_date)%>%
  ungroup()

fledge_mes_rep=fledge_merge%>%
  inner_join(ranked,relationship = "many-to-many")

###  
fledge_pheno_df=fledge_mes_rep%>%
  mutate(age_days=as.numeric(date_diff))%>%
select(RingId, LeftTarsus, Mass, LeftWing, BillLength,age_days, FuniWE ,FHBD512gen, GeneticSex,rank, clutch_merge, Observer, year)%>% ##only seleting info we need 
 group_by(RingId)%>%
 mutate(min_mes=min(age_days))%>%
 filter(!min_mes<0) %>%## removing ids where measurememnt is somehow before estimated hatch date??
 ungroup()%>%
  na.omit(FuniWE)

n_distinct(fledge_pheno_df$RingId) ## 2252

fledge_pheno_df=fledge_pheno_df%>%
  mutate(tarsus_scale=LeftTarsus/100)%>%
  mutate(mass_scale=Mass/100)%>%
  mutate(bill_scale=BillLength/100)%>%
  mutate(wing_scale=LeftWing/100)

fledge_pheno_df$clutch_merge=as.factor(fledge_pheno_df$clutch_merge)
fledge_pheno_df$GeneticSex=as.factor(fledge_pheno_df$GeneticSex)
fledge_pheno_df$RingId=as.factor(fledge_pheno_df$RingId)
fledge_pheno_df$year=as.factor(fledge_pheno_df$year)
fledge_pheno_df$Observer=as.factor(fledge_pheno_df$Observer)
fledge_pheno_df$rank=as.numeric(fledge_pheno_df$rank)
fledge_pheno_df$LeftTarsus=as.numeric(fledge_pheno_df$LeftTarsus)



##################################################################################
##### Running models #######################################################
##########################################################################

## overall growth curve for all left tarsus #####

## using a more simplistic version of the gompertz growth curve to get model priors
fm1 <- nls(tarsus_scale ~ SSgompertz(age_days, asym, b2, b3),
           data = fledge_pheno_df)
summary(fm1)
# Formula: tarsus_scale ~ SSgompertz(age_days, asym, b2, b3)
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# asym 7.1627307  0.0068270 1049.17   <2e-16 ***
#   b2   2.8121146  0.0478616   58.76   <2e-16 ***
#   b3   0.8790505  0.0009756  901.02   <2e-16 ***


# priors_base <- c(
#   prior(normal(2, 2), nlpar = "asym2",  class="b"),##
#   prior(normal(3, 2), nlpar = "b",  class="b"), ## 
#   prior(normal(1, 2), nlpar = "c",  class="b"), ##
#   
#   prior(student_t(3, 0, 2.5), class = "sigma", lb=0))
# 
# growth_model.base=brm(
#   ## model 
#   bf(tarsus_scale ~   asym2 * exp(-b*(c)^age_days), 
#       asym2 + b + c ~ 1, 
#      nl=TRUE),
#   data=fledge_pheno_df, 
#   family = gaussian(),
#   chains = 2,
#   prior = priors_base,
#   control = list(adapt_delta = 0.85),
#   init = 0, 
#   cores = 2
#   
# )

#summary(growth_model.base)
## base model recovers the parameter from SSgompertz which is encouraging
## and plot against real data looks like it fits

priors_1 <- c(
  prior(normal(7, 2), nlpar = "asym2",  class="b"),##
  prior(normal(3, 1), nlpar = "b",  class="b"), ## 
  prior(normal(1, 0.5), nlpar = "c",  class="b"), ## 

  prior(student_t(3, 0, 2.5), class = "sigma", lb=0),
  prior(student_t(3, 0, 2.5), class="sd",nlpar = "asym2", lb=0), # 
  prior(student_t(3, 0, 2.5),  class="sd", nlpar = "b", lb=0),
  prior(student_t(3, 0, 2.5),  class="sd", nlpar = "c", lb=0)
  
   # 
   # prior(cauchy(0, 0.2), class="sd", group="RingId", nlpar = "asym", lb=0), # 
   # prior(cauchy(0, 0.2),  class="sd", group="RingId", nlpar = "b", lb=0),
   # prior(cauchy(0, 0.2),  class="sd", group="RingId", nlpar = "c", lb=0)
)

## should i include RingID as random variable?? its kind of a repeated measure model but not

## adding various other parameters into the model to determine th e best one 
growth_model.1=brm(
  ## model 
  bf(tarsus_scale ~  asym2 * exp(-b*(c)^age_days),
      asym2  ~ 1 + FuniWE + rank + GeneticSex +(1|clutch_merge) + (1|Observer) + (1|year),
       b ~ 1 + FuniWE + rank + GeneticSex + (1|clutch_merge) ,
       c ~ 1 + FuniWE + (1|clutch_merge),
     nl=TRUE),
  data=fledge_pheno_df, 
  family = gaussian(),
  chains = 4,
  prior = priors_1,
  control = list(adapt_delta = 0.85),
  init = 0, 
  cores = 4, 
  iter = 7500, 
  warmup = 2500
  
)


summary(growth_model.1)
plot(growth_model.1)

## converges well and Rhats all = 1

######################## ########################
################### MASS ########################
######################## ########################
fm2 <- nls(mass_scale ~ SSgompertz(age_days, asym, b2, b3),
           data = fledge_pheno_df)
summary(fm2)


plot(fledge_pheno_df$age_days, fledge_pheno_df$mass_scale, xlim = c(0,90))
curve(3.694840*exp(-4.491982*(0.889770)^x), from = 0, to=90, add = TRUE, col="red", lwd=2)


# 
# priors_base_mass <- c(
#   prior(normal(3.5, 2), nlpar = "asym",  class="b"),##
#   prior(normal(4.5, 1), nlpar = "b",  class="b"), ## 
#   prior(normal(1, 0.5), nlpar = "c",  class="b"), ##
#   prior(student_t(3, 0, 2.5), class = "sigma", lb=0))
# 
# growth_model.mass.base=brm(
#   ## model 
#   bf(mass_scale ~ asym * exp(-b*(c)^age_days),
#      asym + b + c ~ 1, 
#      nl=TRUE),
#   data=fledge_pheno_df, 
#   family = gaussian(),
#   prior = priors_base_mass, 
#   chains = 2,
#   control = list(adapt_delta = 0.85),
#   init = 0, 
#   cores = 2
#   
# )
# 
# summary(growth_model.mass.base)




priors_2 <- c(
  prior(normal(3.5, 2), nlpar = "asym",  class="b"),##
  prior(normal(4.5, 1), nlpar = "b",  class="b"), ## 
  prior(normal(1, 0.5), nlpar = "c",  class="b"), ## 
  
  prior(student_t(3, 0, 2.5), class = "sigma", lb=0),
  prior(student_t(3, 0, 2.5), class="sd",nlpar = "asym", lb=0), # 
  prior(student_t(3, 0, 2.5),  class="sd", nlpar = "b", lb=0),
  prior(student_t(3, 0, 2.5),  class="sd", nlpar = "c", lb=0) #,

)

## should i include RingID as random variable?? its kind of a repeated measure model but not

## adding various other parameters into the model to determine th e best one 
growth_model.2=brm(
  ## model 
  bf(mass_scale ~ asym * exp(-b*(c)^age_days),
     asym ~ 1 + FuniWE + rank + GeneticSex +(1|clutch_merge) + (1|Observer)+(1|year),
      b ~ 1 + FuniWE + rank + GeneticSex + (1|clutch_merge),
      c ~ 1 + FuniWE + (1|clutch_merge),
     nl=TRUE),
  data=fledge_pheno_df, 
  family = gaussian(),
  chains = 4,
  prior = priors_2,
  control = list(adapt_delta = 0.85),
  init = 0, 
   cores = 4,
   iter = 7500, 
   warmup = 2500
  
)


summary(growth_model.2)
plot(growth_model.2)


ggplot(fledge_pheno_df, aes(x=age_days, y=mass_scale))+
  geom_point(alpha=0.2)+
  xlim (0,90)+
  theme_classic()+
  labs(x="Age in days", y="Mass (0.01 grams)")+
  stat_function(fun=~ 3.79*exp(-3.20*(0.90)^.x), colour="coral", size=1)+
  stat_function(fun=~3.79*exp(-6.44*(0.90)^.x), colour="blue", size=1)+
  theme(text = element_text(size = 18))+
  annotate("text", x=80, y=4, label="Average individual", colour="coral", size=4)+
  annotate("text", x=80, y=3.3, label="Inbred individual", colour="blue", size=4)
  



################################
 ### Bill ###



## now mass

fm3 <- nls(bill_scale ~ SSgompertz(age_days, asym, b2, b3),
           data = fledge_pheno_df)
summary(fm3)

# 
# priors_base_mass <- c(
#   prior(normal(2, 2), nlpar = "asym",  class="b"),##
#   prior(normal(0, 2), nlpar = "b",  class="b"), ## 
#   prior(normal(0, 2), nlpar = "c",  class="b"), ##
#   prior(student_t(3, 0, 2.5), class = "sigma", lb=0))
# 
# growth_model.bill.base=brm(
#   ## model 
#   bf(bill_scale ~ asym * exp(-b*(c)^age_days),
#      asym + b + c ~ 1, 
#      nl=TRUE),
#   data=fledge_pheno_df, 
#   family = gaussian(),
#   prior = priors_base_mass, 
#   chains = 2,
#   control = list(adapt_delta = 0.85),
#   init = 0, 
#   cores = 2
#   
# )
# 
# summary(growth_model.bill.base)
# 


priors_3 <- c(
  prior(normal(2, 2), nlpar = "asym",  class="b"),##
  prior(normal(2, 1), nlpar = "b",  class="b"), ## 
  prior(normal(1, 0.5), nlpar = "c",  class="b"), ## 
  
  prior(student_t(3, 0, 2.5), class = "sigma", lb=0),
  prior(student_t(3, 0, 2.5), class="sd",nlpar = "asym", lb=0), # 
  prior(student_t(3, 0, 2.5),  class="sd", nlpar = "b", lb=0)

)

## should i include RingID as random variable?? its kind of a repeated measure model but not

## adding various other parameters into the model to determine th e best one 
growth_model.3=brm(
  ## model 
  bf(bill_scale ~ asym * exp(-exp(b)*exp(c)^age_days),
     asym ~ 1 + FuniWE + rank + GeneticSex +(1|clutch_merge) + (1|Observer)+(1|year),
     b ~ 1 + FuniWE + rank + GeneticSex + (1|clutch_merge),
     c ~ 1 + FuniWE, # clutch here has no effect and converges badly so remove 
     nl=TRUE),
  data=fledge_pheno_df, 
  family = gaussian(),
  chains = 4,
  prior = priors_3,
  control = list(adapt_delta = 0.85),
  init = 0, 
  cores = 4, 
  iter = 7500, 
  warmup = 2500
  
)


summary(growth_model.3)
plot(growth_model.3)





##############################################################


fm4 <- nls(wing_scale ~ SSgompertz(age_days, asym, b2, b3),
           data = fledge_pheno_df)
summary(fm4)

priors_4 <- c(
  prior(normal(3, 2), nlpar = "asym2",  class="b"),##
  prior(normal(4, 1), nlpar = "b",  class="b"), ## 
  prior(normal(1, 0.5), nlpar = "c",  class="b"), ## 
  
  prior(student_t(3, 0, 2.5), class = "sigma", lb=0),
  prior(student_t(3, 0, 2.5), class="sd",nlpar = "asym2", lb=0), # 
  prior(student_t(3, 0, 2.5),  class="sd", nlpar = "b", lb=0),
  prior(student_t(3, 0, 2.5),  class="sd", nlpar = "c", lb=0)
  

)


growth_model.4=brm(
  ## model 
  bf(wing_scale ~  asym2 * exp(-b*(c)^age_days),
     asym2  ~ 1 + FuniWE + rank + GeneticSex +(1|clutch_merge) + (1|Observer) + (1|year),
     b ~ 1 + FuniWE + rank + GeneticSex + (1|clutch_merge) ,
     c ~ 1 + FuniWE,
     nl=TRUE),
  data=fledge_pheno_df, 
  family = gaussian(),
  chains = 4,
  prior = priors_4,
  control = list(adapt_delta = 0.85),
  init = 0, 
  cores = 4,
   iter = 7500, 
  warmup = 2500
  
)


summary(growth_model.4)
plot(growth_model.4)




# saveRDS(growth_model.1,file="Inbreeding_depression_owls/Model_outputs/Growth_tarsus_funi.RDS")
# 
# saveRDS(growth_model.2,file="Inbreeding_depression_owls/Model_outputs/Growth_mass_funi.RDS")
# 
# saveRDS(growth_model.3,file="Inbreeding_depression_owls/Model_outputs/Growth_bill_funi.RDS")
# 
# saveRDS(growth_model.4,file="Inbreeding_depression_owls/Model_outputs/Growth_wing_funi.RDS")
