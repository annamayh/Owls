### ID in growth rate of wing ###

library(brms,  lib = "/users/ahewett1/R")
library(corpcor,  lib = "/users/ahewett1/R")
library(readr,  lib = "/users/ahewett1/R")

##################################################################################
########################### ~~ wing ~~ #############################################
#####################################################################################

args = commandArgs(trailingOnly = TRUE)
scratch = as.character(args[1]) 

wing_df=read.table("./input_dfs/wing_all_pheno_df.txt",sep=",", header=T)


wing_df$clutch_merge=as.factor(wing_df$clutch_merge)
wing_df$sex=as.factor(wing_df$sex)
wing_df$RingId=as.factor(wing_df$RingId)
wing_df$year=as.factor(wing_df$year)
wing_df$Observer=as.factor(wing_df$Observer)
wing_df$nestboxID=as.factor(wing_df$nestboxID)

wing_df$rank=as.numeric(wing_df$rank)



prior_wing_gr<- c(
  prior(normal(3, 0.5), nlpar = "asym",  class="b", coef="Intercept"),##
  prior(normal(4, 1), nlpar = "b",  class="b", coef="Intercept"), ## 
  prior(normal(0.95, 0.5), nlpar = "c",  class="b", coef="Intercept"), ## 
  
  prior(normal(0, 0.5), nlpar = "asym",  class="b"),## more stringent priors for indivduals effects
  prior(normal(0, 1), nlpar = "b",  class="b"), ## 
  prior(normal(0, 1), nlpar = "c",  class="b"), ## 
  
  prior(student_t(3, 0, 2), class = "sigma", lb=0),
  prior(student_t(3, 0, 0.5), class="sd",nlpar = "asym", lb=0), # 
  prior(student_t(3, 0, 1),  class="sd", nlpar = "b", lb=0),
  prior(student_t(3, 0, 1),  class="sd", nlpar = "c", lb=0), 
  
  
  prior(cauchy(0, 0.5), class="sd", group="RingId", nlpar = "asym", lb=0), #
  prior(cauchy(0, 0.5),  class="sd", group="RingId", nlpar = "b", lb=0),
  prior(cauchy(0, 0.2),  class="sd", group="RingId", nlpar = "c", lb=0)
  
)



growth_wing.mod=brm(
  ## model 
  bf(wing_scale ~ asym * exp(-b*(c)^age_days),
     asym ~ 1 + FuniWE + rank + sex + (1|clutch_merge) + (1|Observer) + (1|year) + (1|month) + (1|RingId) + (1|nestboxID),
     b ~ 1 + FuniWE + rank + sex + (1|clutch_merge) + (1|RingId),
     c ~ 1 + FuniWE + rank + sex + (1|RingId), 
     nl=TRUE),
  data=wing_df, 
  family = gaussian(),
  chains = 4,
  prior = prior_wing_gr,
  control = list(adapt_delta = 0.99),
  init = 0, 
  cores = 4,
  iter = 30000, 
  warmup = 10000, 
  thin=10
  
)





### save into outputs folder
saveRDS(growth_wing.mod,file=paste0(scratch,"3.3_wing_scaled_gr_30k_stringentpr.RDS")) ##

#summary(growth_wing.mod)

