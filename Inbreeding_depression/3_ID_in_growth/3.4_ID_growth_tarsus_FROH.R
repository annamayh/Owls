### ID in growth rate of tarsus ###

.libPaths(c("/work/FAC/FBM/DEE/jgoudet/barn_owl/ahewett/R", .libPaths()))

library(brms)
library(dplyr)

##################################################################################
########################### ~~ tarsus ~~ #############################################
#####################################################################################

args = commandArgs(trailingOnly = TRUE)
scratch = as.character(args[1]) 

tarsus_df=read.table("./input_dfs/tarsus_all_pheno_df.txt",sep=",", header=T)

# sorting out variables
tarsus_df$clutch_merge=as.factor(tarsus_df$clutch_merge)
tarsus_df$sex=as.factor(tarsus_df$sex)
tarsus_df$RingId=as.factor(tarsus_df$RingId)
tarsus_df$year=as.factor(tarsus_df$year)
tarsus_df$Observer=as.factor(tarsus_df$Observer)
tarsus_df$nestboxID=as.factor(tarsus_df$nestboxID)
tarsus_df$rank=as.numeric(tarsus_df$rank)


prior_tarsus_gr<- c(
  prior(normal(50, 5), nlpar = "asym1",  class="b", coef="Intercept"),## priors for intercept expectations
  prior(normal(660, 20), nlpar = "asym",  class="b", coef="Intercept"),## priors for intercept expectations
  prior(normal(2, 1), nlpar = "b",  class="b", coef="Intercept"), ## 
  prior(normal(0.9, 0.2), nlpar = "c",  class="b", coef="Intercept"), ## 
  
  prior(normal(0, 20), nlpar = "asym",  class="b"),## priors for fixed effects (deviaition from intercept)
  #prior(normal(0, 5), nlpar = "b",  class="b"), ## 
  prior(normal(0, 1), nlpar = "c",  class="b"), ## 
  
  prior(student_t(3, 0, 20), class = "sigma", lb=0), ## priors for random effects 
  prior(student_t(3, 0, 20), class="sd",nlpar = "asym", lb=0), # 
  #prior(student_t(3, 0, 5),  class="sd", nlpar = "b", lb=0),
  prior(student_t(3, 0, 1),  class="sd", nlpar = "c", lb=0), 
  
  #prior(cauchy(0, 5), class="sd", group="RingId", nlpar = "asym1", lb=0),
  prior(cauchy(0, 5), class="sd", group="RingId", nlpar = "asym", lb=0), # priors for within id effect 
  #prior(cauchy(0, 0.5),  class="sd", group="RingId", nlpar = "b", lb=0),
  prior(cauchy(0, 0.2),  class="sd", group="RingId", nlpar = "c", lb=0)
  
)


growth_tarsus.mod=brm(
  ## model 
  bf(LeftTarsus ~ asym1 + (asym * exp(-b*(c)^age_days)),
     asym1 ~ 1, ## assuming there is no effect at age 0 
     asym ~ 1 + FROH + rank + sex + (1|clutch_merge) + (1|Observer) + (1|year) + (1|month) + (1|RingId) + (1|nestboxID),
     b ~ 1 ,
     c ~ 1 + FROH + rank + sex + (1|RingId), 
     nl=TRUE),
  data=tarsus_df, 
  family = gaussian(),
  chains = 4,
  prior = prior_tarsus_gr,
  control = list(adapt_delta = 0.99),
  init = 0, 
  cores = 4,
  iter = 55000, 
  warmup = 15000, 
  thin=10
  
)





### save into outputs folder
saveRDS(growth_tarsus.mod,file=paste0(scratch,"3.4_tarsus_gr_FROH_subset.RDS")) ##

summary(growth_tarsus.mod)

