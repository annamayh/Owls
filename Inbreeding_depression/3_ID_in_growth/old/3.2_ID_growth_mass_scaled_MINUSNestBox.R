### ID in mass growth ###

library(tidyverse)
library(brms)
library(corpcor)

##################################################################################
########################### ~~ mass ~~ #############################################
#####################################################################################

mass_df=read.table("./input_dfs/mass_all_pheno_df.txt",sep=",", header=T)%>% ##
  mutate(mass_scale=Mass/100) ## scale mass to make model converage better

mass_df$clutch_merge=as.factor(mass_df$clutch_merge)
mass_df$sex=as.factor(mass_df$sex)
mass_df$RingId=as.factor(mass_df$RingId)
mass_df$year=as.factor(mass_df$year)
mass_df$Observer=as.factor(mass_df$Observer)
mass_df$nestboxID=as.factor(mass_df$nestboxID)

mass_df$rank=as.numeric(mass_df$rank)


prior_mass_gr<- c(
  prior(normal(3.6, 2), nlpar = "asym",  class="b"),##
  prior(normal(4, 0.5), nlpar = "b",  class="b"), ## 
  prior(normal(0.89, 0.2), nlpar = "c",  class="b"), ## 
  
  prior(student_t(3, 0, 6), class = "sigma", lb=0),
  prior(student_t(3, 0, 6), class="sd",nlpar = "asym", lb=0), # 
  prior(student_t(3, 0, 6),  class="sd", nlpar = "b", lb=0),
  prior(student_t(3, 0, 6),  class="sd", nlpar = "c", lb=0), 
  
  
  prior(cauchy(0, 1), class="sd", group="RingId", nlpar = "asym", lb=0), #
  prior(cauchy(0, 0.2),  class="sd", group="RingId", nlpar = "b", lb=0),
  prior(cauchy(0, 0.2),  class="sd", group="RingId", nlpar = "c", lb=0)
  
)



growth_mass.mod=brm(
  ## model 
  bf(mass_scale ~ asym * exp(-b*(c)^age_days),
     asym ~ 1 + FuniWE + rank + sex + (1|clutch_merge) + (1|Observer) + (1|year) + (1|month) + (1|RingId),
     b ~ 1 + FuniWE + rank + sex + (1|clutch_merge) + (1|RingId),
     c ~ 1 + FuniWE + rank + sex + (1|RingId), 
     nl=TRUE),
  data=mass_df, 
  family = gaussian(),
  chains = 4,
  prior = prior_mass_gr,
  control = list(adapt_delta = 0.95),
  init = 0, 
  cores = 4,
  iter = 15000, 
  warmup = 5000, 
  thin=5
  
)


summary(growth_mass.mod)


### save into outputs folder
saveRDS(growth_mass.mod,file="./outputs/3_growth/3.2_mass_scaled_gr_MINUS_nestbox.RDS") ##

