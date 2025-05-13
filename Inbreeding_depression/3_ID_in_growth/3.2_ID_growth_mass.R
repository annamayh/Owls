### ID in mass growth ###

.libPaths(c("/work/FAC/FBM/DEE/jgoudet/barn_owl/ahewett/R", .libPaths())) #specify library cluster path

library(brms)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)
scratch = as.character(args[1]) 

##################################################################################
########################### ~~ mass ~~ #############################################
#####################################################################################

mass_df=read.table("./input_dfs/mass_all_pheno_df.txt",sep=",", header=T)

mass_df$clutch_merge=as.factor(mass_df$clutch_merge)
mass_df$sex=as.factor(mass_df$sex)
mass_df$RingId=as.factor(mass_df$RingId)
mass_df$year=as.factor(mass_df$year)
mass_df$Observer=as.factor(mass_df$Observer)
mass_df$nestboxID=as.factor(mass_df$nestboxID)

mass_df$rank=as.numeric(mass_df$rank)


prior_mass_gr<- c(
  prior(normal(360, 20), nlpar = "asym",  class="b", coef="Intercept"),## 
  prior(normal(4, 1), nlpar = "b",  class="b", coef="Intercept"), ## 
  prior(normal(0.89, 0.5), nlpar = "c",  class="b", coef="Intercept"), ## 
  
  prior(normal(0, 20), nlpar = "asym",  class="b"),## more stringent priors for indivduals effects
  #prior(normal(0, 5), nlpar = "b",  class="b"), ## 
  prior(normal(0, 1), nlpar = "c",  class="b"), ## 
  
  
  prior(student_t(3, 0, 20), class = "sigma", lb=0),
  prior(normal(0, 20), class="sd",nlpar = "asym", lb=0), # 
  #prior(normal(0, 5),  class="sd", nlpar = "b", lb=0),
  prior(normal(0, 1),  class="sd", nlpar = "c", lb=0), 
  
  
  prior(normal(0, 10), class="sd", group="RingId", nlpar = "asym", lb=0), #
  # prior(normal(0, 1),  class="sd", group="RingId", nlpar = "b", lb=0),
  prior(normal(0, 0.2),  class="sd", group="RingId", nlpar = "c", lb=0)
  
)



growth_mass.mod=brm(
  ## model 
  bf(Mass ~ asym * exp(-b*(c)^age_days),
     asym ~ 1 + FuniW + rank + sex + (1|clutch_merge) + (1|Observer) + (1|year) + (1|month) + (1|RingId) + (1|nestboxID),
     b ~ 1 ,
     c ~ 1 + FuniW + rank + sex + (1|RingId), 
     nl=TRUE),
  data=mass_df, 
  family = gaussian(),
  chains = 4,
  prior = prior_mass_gr,
  control = list(adapt_delta = 0.99),
  init = 0, 
  cores = 4,
  iter = 70000, 
  warmup = 30000, 
  thin=5
  
)



### save into outputs folder
saveRDS(growth_mass.mod,file=paste0(scratch,"3.2_mass_gr_Funi_subset.RDS")) ##

summary(growth_mass.mod)
