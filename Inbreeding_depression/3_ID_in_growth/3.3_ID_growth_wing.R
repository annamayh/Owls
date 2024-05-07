### ID in growth rate of wing ###

library(tidyverse)
library(brms)
library(corpcor)

##################################################################################
########################### ~~ wing ~~ #############################################
#####################################################################################

wing_df=read.table("./input_dfs/wing_all_pheno_df.txt",sep=",", header=T)%>% ##
  mutate(wing_scale=LeftWing/100) ## scale wing to make model converage better

wing_df$clutch_merge=as.factor(wing_df$clutch_merge)
wing_df$sex=as.factor(wing_df$sex)
wing_df$RingId=as.factor(wing_df$RingId)
wing_df$year=as.factor(wing_df$year)
wing_df$Observer=as.factor(wing_df$Observer)
wing_df$nestboxID=as.factor(wing_df$nestboxID)

wing_df$rank=as.numeric(wing_df$rank)



prior_wing_gr<- c(
  prior(normal(298, 10), nlpar = "asym",  class="b"),##
  prior(normal(4, 0.5), nlpar = "b",  class="b"), ## 
  prior(normal(0.95, 0.2), nlpar = "c",  class="b"), ## 
  
  prior(student_t(3, 0, 90), class = "sigma", lb=0),
  prior(student_t(3, 0, 90), class="sd",nlpar = "asym", lb=0), # 
  prior(student_t(3, 0, 90),  class="sd", nlpar = "b", lb=0),
  prior(student_t(3, 0, 90),  class="sd", nlpar = "c", lb=0), 
  
  
  prior(cauchy(0, 10), class="sd", group="RingId", nlpar = "asym", lb=0), #
  prior(cauchy(0, 0.2),  class="sd", group="RingId", nlpar = "b", lb=0),
  prior(cauchy(0, 0.2),  class="sd", group="RingId", nlpar = "c", lb=0)
  
)



growth_wing.mod=brm(
  ## model 
  bf(LeftWing ~ asym * exp(-b*(c)^age_days),
     asym ~ 1 + FuniWE + rank + sex + (1|clutch_merge) + (1|Observer) + (1|year) + (1|month) + (1|RingId) + (1|nestboxID),
     b ~ 1 + FuniWE + rank + sex + (1|clutch_merge) + (1|RingId),
     c ~ 1 + FuniWE + rank + sex + (1|RingId), 
     nl=TRUE),
  data=wing_df, 
  family = gaussian(),
  chains = 4,
  prior = prior_wing_gr,
  control = list(adapt_delta = 0.95),
  init = 0, 
  cores = 4,
  iter = 15000, 
  warmup = 5000, 
  thin=5
  
)


summary(growth_wing.mod)

plot(growth_wing.mod)


### save into outputs folder
saveRDS(growth_wing.mod,file="./outputs/3_growth/3.3_wing_unscaled_gr_15k.RDS") ##


