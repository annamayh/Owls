
setwd("/Users/ahewett1/Documents")

## combined jueniles and adult records in to one big model 

## for tarsus and bill 100k iter was a bit of overkill .. can probable run for much less 

grm_funi_tarsus=readRDS("Inbreeding_depression_owls/Model_outputs/tarsus_all_GRM_Funi.RDS")

plot(grm_funi_tarsus)

summary(grm_funi_tarsus)

vt_animal <- (VarCorr(grm_funi_tarsus, summary = FALSE)$RingId$sd)^2

vt_year <- (VarCorr(grm_funi_tarsus, summary = FALSE)$year$sd)^2
vt_clutch <- (VarCorr(grm_funi_tarsus, summary = FALSE)$clutch_merge$sd)^2
vt_obs <- (VarCorr(grm_funi_tarsus, summary = FALSE)$Observer$sd)^2
vt_pe <- (VarCorr(grm_funi_tarsus, summary = FALSE)$RingId_pe$sd)^2

vt_r <- (VarCorr(grm_funi_tarsus, summary = FALSE)$residual$sd)^2

h.tarsus <- as.mcmc(vt_animal / (vt_animal + vt_year + vt_clutch + vt_obs + vt_pe + vt_r))
summary(h.tarsus)

# Mean             SD       Naive SE Time-series SE 
# 0.1066648      0.0131436      0.0001073      0.0001088 
# 


## heritability of bill length ##

grm_funi_bill=readRDS("Inbreeding_depression_owls/Model_outputs/bill_all_GRM_Funi.RDS")

plot(grm_funi_bill)

summary(grm_funi_bill)

vb_animal<- (VarCorr(grm_funi_bill, summary = FALSE)$RingId$sd)^2

vb_year <- (VarCorr(grm_funi_bill, summary = FALSE)$year$sd)^2
vb_clutch <- (VarCorr(grm_funi_bill, summary = FALSE)$clutch_merge$sd)^2
vb_obs <- (VarCorr(grm_funi_bill, summary = FALSE)$Observer$sd)^2
vb_pe <- (VarCorr(grm_funi_bill, summary = FALSE)$RingId_pe$sd)^2

vb_r <- (VarCorr(grm_funi_bill, summary = FALSE)$residual$sd)^2

h.bill <- as.mcmc(vb_animal / (vb_animal + vb_year + vb_clutch + vb_obs + vb_pe + vb_r))
summary(h.bill)

# Mean             SD       Naive SE Time-series SE 
# 0.1745728      0.0191503      0.0001564      0.0001661 
