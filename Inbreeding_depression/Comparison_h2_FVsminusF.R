library(brms)

setwd("/Users/ahewett1/Documents")

## combined jueniles and adult records in to one big model 

## for tarsus and bill 100k iter was a bit of overkill .. can probable run for much less 

grm_funi_bill=readRDS("Inbreeding_depression_owls/Model_outputs/1.1.ID_bill_GRM_Funi_unscaled.RDS")

plot(grm_funi_bill)

summary(grm_funi_bill)

vt_animal <- (VarCorr(grm_funi_bill, summary = FALSE)$RingId$sd)^2

vt_year <- (VarCorr(grm_funi_bill, summary = FALSE)$year$sd)^2
vt_clutch <- (VarCorr(grm_funi_bill, summary = FALSE)$clutch_merge$sd)^2
# vt_obs <- (VarCorr(grm_funi_bill, summary = FALSE)$Observer$sd)^2
vt_rank <- (VarCorr(grm_funi_bill, summary = FALSE)$rank$sd)^2
vt_nestb <- (VarCorr(grm_funi_bill, summary = FALSE)$nestboxID$sd)^2

vt_pe <- (VarCorr(grm_funi_bill, summary = FALSE)$RingId_pe$sd)^2

vt_r <- (VarCorr(grm_funi_bill, summary = FALSE)$residual$sd)^2

h.bill <- as.mcmc(vt_animal / (vt_animal + vt_year + vt_clutch  + vt_rank+vt_pe + vt_nestb+vt_r))
summary(h.bill)

# Iterations = 1:32000
# Thinning interval = 1 
# Number of chains = 1 
# Sample size per chain = 32000 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#   Mean             SD       Naive SE Time-series SE 
# 0.1874903      0.0214911      0.0001201      0.0001858 
# 
# 2. Quantiles for each variable:
#   
#   2.5%    25%    50%    75%  97.5% 
# 0.1442 0.1734 0.1880 0.2021 0.2283 

### now h2 without Funi ###


grm_bill=readRDS("Inbreeding_depression_owls/Model_outputs/1.1.ID_bill_GRM_MINUSFuni_unscaled.RDS")

plot(grm_bill)

summary(grm_bill)

vt_animal <- (VarCorr(grm_bill, summary = FALSE)$RingId$sd)^2

vt_year <- (VarCorr(grm_bill, summary = FALSE)$year$sd)^2
vt_clutch <- (VarCorr(grm_bill, summary = FALSE)$clutch_merge$sd)^2
#vt_obs <- (VarCorr(grm_bill, summary = FALSE)$Observer$sd)^2
vt_rank <- (VarCorr(grm_bill, summary = FALSE)$rank$sd)^2
vt_nestb <- (VarCorr(grm_bill, summary = FALSE)$nestboxID$sd)^2
vt_pe <- (VarCorr(grm_bill, summary = FALSE)$RingId_pe$sd)^2

vt_r <- (VarCorr(grm_bill, summary = FALSE)$residual$sd)^2

h.bill_noF <- as.mcmc(vt_animal / (vt_animal + vt_year + vt_clutch + vt_pe + vt_r+ vt_rank+ vt_nestb))
summary(h.bill_noF)

# Iterations = 1:32000
# Thinning interval = 1 
# Number of chains = 1 
# Sample size per chain = 32000 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#   Mean             SD       Naive SE Time-series SE 
# 0.1873828      0.0208905      0.0001168      0.0001870 
# 
# 2. Quantiles for each variable:
#   
#   2.5%    25%    50%    75%  97.5% 
# 0.1450 0.1737 0.1877 0.2017 0.2269 