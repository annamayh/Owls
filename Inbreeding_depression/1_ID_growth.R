## df for modelling inbreeding depression in growth using estimated hatch date##

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
select(RingId, LeftTarsus, age_days, FuniWE ,FHBD512gen, GeneticSex,rank, clutch_merge, Observer, year)%>% ##only seleting info we need 
 group_by(RingId)%>%
 mutate(min_mes=min(age_days))%>%
 filter(!min_mes<0) %>%## removing ids where measurememnt is somehow before estimated hatch date??
 #filter(!min_mes>30)%>% ## remove ids we dont have info for before 30 days i.e. the growth period
 ungroup()%>%
  na.omit(LeftTarsus,FuniWE)

n_distinct(fledge_pheno_df$RingId) ## will just stick to juveniles because sa,e number of ids 

fledge_pheno_df=fledge_pheno_df%>%
  mutate(tarsus_scale=LeftTarsus/100)%>%
  mutate(funi_scale=FuniWE+1)

#fledge_pheno_df[99,]




## overall growth curve for all left tarsus
fm1 <- nls(tarsus_scale ~ SSgompertz(age_days, asym, b2, b3),
           data = fledge_pheno_df)
summary(fm1)
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# asym 7.1948339  0.0070856 1015.41   <2e-16 ***
#   b2   2.2608293  0.0263329   85.86   <2e-16 ***
#   b3   0.8899650  0.0007425 1198.61   <2e-16 ***

plot(fledge_pheno_df$age_days, fledge_pheno_df$tarsus_scale, xlim = c(0,90))
curve(7.19*exp(-2.26*0.89^x), from = 0, to=90, add = TRUE, col="red", lwd=2)

fledge_pheno_df$clutch_merge=as.factor(fledge_pheno_df$clutch_merge)
fledge_pheno_df$GeneticSex=as.factor(fledge_pheno_df$GeneticSex)
fledge_pheno_df$RingId=as.factor(fledge_pheno_df$RingId)
fledge_pheno_df$year=as.factor(fledge_pheno_df$year)
fledge_pheno_df$Observer=as.factor(fledge_pheno_df$Observer)
fledge_pheno_df$rank=as.numeric(fledge_pheno_df$rank)
fledge_pheno_df$LeftTarsus=as.numeric(fledge_pheno_df$LeftTarsus)


## start with non-linear model of growth
# Priors


# prior1=c( 
#   prior(cauchy(0, 20), class="sd", group="RingId", nlpar = "asym"),
#   prior(cauchy(0, 0.5),  class="sd", group="RingId", nlpar = "b"),
#   prior(cauchy(0, 0.5),  class="sd", group="RingId", nlpar = "c"),
#   
#   prior(normal(700, 20), nlpar = "asym", coef="Intercept"),
#   prior(normal(0, 5), nlpar = "b",coef="Intercept"),
#   prior(normal(0, 1), nlpar = "c",coef="Intercept"))
# 
# 
# 
# ## with upper and lower bounds on the growth curve parameters
# priors_2 <- c(
#   prior(normal(700, 20), nlpar = "asym", lb=500, ub=950, coef="Intercept"),## resonable weights 
#   prior(normal(2.5, 0.5), nlpar = "b", lb=0, coef="Intercept"), ## lower bound of 0 because is -ve in the expression
#   prior(normal(0.8, 1), nlpar = "c", lb=0, ub=0.99,coef="Intercept"), 
# 
#   prior(cauchy(0, 20), class="sd", group="RingId", nlpar = "asym"),
#   prior(cauchy(0, 0.5),  class="sd", group="RingId", nlpar = "b"),
#   prior(cauchy(0, 0.5),  class="sd", group="RingId", nlpar = "c")
# )


priors_3 <- c(
  prior(normal(7, 2), nlpar = "asym"),## resonable weight 
  prior(normal(2, 1), nlpar = "b", lb=0), ## lower bound of 0 because is -ve in the expression
  prior(normal(0.8, 0.5), nlpar = "c", lb=0, ub=0.99), 
  
  prior(cauchy(0, 5), class="sd", group="RingId", nlpar = "asym"), # assume minimal variation within ids
  prior(cauchy(0, 0.5),  class="sd", group="RingId", nlpar = "b"),
  prior(cauchy(0, 0.5),  class="sd", group="RingId", nlpar = "c")
)


growth_model=brm(
  ## model 
  bf(tarsus_scale ~ asym * exp(-b*c^age_days),
              asym ~ 1 + funi_scale + GeneticSex+ (1|RingId),
              b ~  1 + funi_scale + (1|RingId), 
              c ~  1 + funi_scale + (1|RingId),
              nl=TRUE),
              
  data=fledge_pheno_df, 
  family = gaussian(),
  chains = 2,
  prior = priors_3,
  control = list(adapt_delta = 0.85),
  init = 0

)


summary(growth_model)
plot(growth_model)

## base model to build up from 
# bf(tarsus_scale ~ asym * exp(-b*c^age_days),
#    asym ~ 1 + FHBD512gen + GeneticSex+ (1|RingId),
#    b ~  1 + FHBD512gen + (1|RingId), 
#    c ~  1 + FHBD512gen + (1|RingId),
#    nl=TRUE)


## add sex in
growth_model.1clutch=brm(
  ## model 
  bf(tarsus_scale ~ asym * exp(-b*c^age_days),
     asym ~ 1 + FHBD512gen + (1|clutch_merge) + (1|RingId),
     b ~  1 + FHBD512gen + (1|clutch_merge) + (1|RingId), 
     c ~  1 + FHBD512gen + (1|RingId), #clutch and sex have no effect here
     nl=TRUE),
  
  data=fledge_pheno_df, 
  family = gaussian(),
  chains = 2,
  prior = priors_3,
  control = list(adapt_delta = 0.85),
  init = 0
  
)


summary(growth_model.1sex)
plot(growth_model.1clutch)

