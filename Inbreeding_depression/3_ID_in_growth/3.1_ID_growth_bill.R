### ID in growth rate of bill ###
.libPaths(c("/work/FAC/FBM/DEE/jgoudet/barn_owl/ahewett/R", .libPaths())) #specify library cluster path

library(brms)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)
scratch = as.character(args[1]) 

##################################################################################
########################### ~~ bill ~~ #############################################
#####################################################################################

bill_df=read.table("./input_dfs/bill_all_pheno_df.txt",sep=",", header=T)

bill_df$clutch_merge=as.factor(bill_df$clutch_merge)
bill_df$sex=as.factor(bill_df$sex)
bill_df$RingId=as.factor(bill_df$RingId)
bill_df$year=as.factor(bill_df$year)
bill_df$Observer=as.factor(bill_df$Observer)
bill_df$nestboxID=as.factor(bill_df$nestboxID)

bill_df$rank=as.numeric(bill_df$rank)



prior_bill_gr<- c(
  prior(normal(184, 10), nlpar = "asym",  class="b", coef="Intercept"),##
  prior(normal(1, 2), nlpar = "b",  class="b", coef="Intercept"), ## 
  prior(normal(0.93, 0.5), nlpar = "c",  class="b", coef="Intercept"), ## 
  
  prior(normal(0, 10), nlpar = "asym",  class="b"),## 
  #prior(normal(0, 2), nlpar = "b",  class="b"), ## 
  prior(normal(0, 1), nlpar = "c",  class="b"), ## 
  
  prior(student_t(3, 0, 10), class = "sigma", lb=0),
  prior(student_t(3, 0, 10), class="sd",nlpar = "asym", lb=0), # 
 # prior(student_t(3, 0, 2),  class="sd", nlpar = "b", lb=0),
  prior(student_t(3, 0, 1),  class="sd", nlpar = "c", lb=0), 

  
  prior(normal(0, 2), class="sd", group="RingId", nlpar = "asym", lb=0), #
  #prior(cauchy(0, 0.5),  class="sd", group="RingId", nlpar = "b", lb=0),
  prior(normal(0, 0.1),  class="sd", group="RingId", nlpar = "c", lb=0)
  
)



growth_bill.mod=brm(
  ## model 
  bf(BillLength ~ asym * exp(-b*(c)^age_days),
     asym ~ 1 + FuniW + rank + sex + (1|clutch_merge) + (1|Observer) + (1|year) + (1|month) + (1|RingId) + (1|nestboxID),
     b ~ 1 ,
     c ~ 1 + FuniW + rank + sex + (1|RingId), 
     nl=TRUE),
  data=bill_df, 
  family = gaussian(),
  chains = 4,
  prior = prior_bill_gr,
  control = list(adapt_delta = 0.98),
  init = 0, 
  cores = 4,
  iter = 70000, 
  warmup = 30000, 
  thin=5
  
)


saveRDS(growth_bill.mod,file=paste0(scratch,"3.1.bill_gr_Funi_subset.RDS")) ##

summary(growth_bill.mod)

