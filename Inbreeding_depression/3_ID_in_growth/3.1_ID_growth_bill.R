### ID in growth rate of bill ###

library(tidyverse)
library(brms)
library(corpcor)

##################################################################################
########################### ~~ bill ~~ #############################################
#####################################################################################

bill_df=read.table("Inbreeding_depression_owls/pheno_df/bill_all_pheno_df.txt",sep=",", header=T)%>% ##
  mutate(bill_scale=BillLength/100) ## scale bill to make model converage better

bill_df$clutch_merge=as.factor(bill_df$clutch_merge)
bill_df$sex=as.factor(bill_df$sex)
bill_df$RingId=as.factor(bill_df$RingId)
bill_df$year=as.factor(bill_df$year)
bill_df$Observer=as.factor(bill_df$Observer)
bill_df$nestboxID=as.factor(bill_df$nestboxID)

bill_df$rank=as.numeric(bill_df$rank)

#check dist 

plot(bill_df$age_days, bill_df$BillLength, xlim = c(0,90))
curve(184*exp(-0.99*0.932^x), from = 0, to=90, add = TRUE, col="red", lwd=2)





prior_bill_gr<- c(
  prior(normal(184, 10), nlpar = "asym",  class="b"),##
  prior(normal(1, 0.5), nlpar = "b",  class="b"), ## 
  prior(normal(0.93, 0.2), nlpar = "c",  class="b"), ## 
  
  prior(student_t(3, 0, 17), class = "sigma", lb=0),
  prior(student_t(3, 0, 17), class="sd",nlpar = "asym", lb=0), # 
  prior(student_t(3, 0, 17),  class="sd", nlpar = "b", lb=0),
  prior(student_t(3, 0, 17),  class="sd", nlpar = "c", lb=0), 

  
  prior(cauchy(0, 10), class="sd", group="RingId", nlpar = "asym", lb=0), #
  prior(cauchy(0, 0.2),  class="sd", group="RingId", nlpar = "b", lb=0),
  prior(cauchy(0, 0.2),  class="sd", group="RingId", nlpar = "c", lb=0)
  
)



growth_bill.mod=brm(
  ## model 
  bf(BillLength ~ asym * exp(-b*(c)^age_days),
     asym ~ 1 + FuniWE + rank + sex + (1|clutch_merge) + (1|Observer) + (1|year) + (1|month) + (1|RingId) + (1|nestboxID),
     b ~ 1 + FuniWE + rank + sex + (1|clutch_merge) + (1|RingId),
     c ~ 1 + FuniWE + rank + sex + (1|RingId), 
     nl=TRUE),
  data=bill_df, 
  family = gaussian(),
  chains = 4,
  prior = prior_bill_gr,
  control = list(adapt_delta = 0.95),
  init = 0, 
  cores = 4,
  iter = 15000, 
  warmup = 5000, 
  thin=5
  
)


summary(growth_bill.mod)

plot(growth_bill.mod)



saveRDS(growth_bill.mod,file="Inbreeding_depression_owls/Model_outputs/growth/bill_unscaled_gr_15k.RDS") ##





## now filter for ages <100 days 
## ####.   results pretty simular ######
# 
# bill_df_filt=bill_df%>%
#   filter(age_days<100) ## removes ~1000 records
# 
# growth_bill.mod.filter=brm(
#   ## model 
#   bf(BillLength ~ asym * exp(-b*(c)^age_days),
#      asym ~ 1 + FuniWE + rank + sex + (1|clutch_merge) + (1|Observer) + (1|year) + (1|month) + (1|RingId) + (1|nestboxID),
#      b ~ 1 + FuniWE + rank + sex + (1|clutch_merge) + (1|RingId),
#      c ~ 1 + FuniWE + rank + sex + (1|RingId), 
#      nl=TRUE),
#   data=bill_df_filt, 
#   family = gaussian(),
#   chains = 4,
#   prior = prior_bill_gr,
#   control = list(adapt_delta = 0.95),
#   init = 0, 
#   cores = 4,
#   iter = 15000, 
#   warmup = 5000, 
#   thin=5
#   
# )
# 
# 
# summary(growth_bill.mod.filter)
# 
# plot(growth_bill.mod.filter)
