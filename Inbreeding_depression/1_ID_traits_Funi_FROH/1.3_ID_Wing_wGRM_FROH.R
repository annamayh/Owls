library(tidyverse)
library(brms)
library(corpcor)

##################################################################################
########################### ~~ WING LENGTH ~~ #############################################
#####################################################################################

wing_df=read.table("./input_dfs/wing_all_pheno_df.txt",sep=",", header=T)%>%
  mutate(age_acc=298*exp(-4.19*0.945^age_days)) %>% ## account for age unsing gompertz growth
  mutate(RingId_pe=RingId) # add permanent environment for repeated measures 


wing_df$clutch_merge=as.factor(wing_df$clutch_merge)
wing_df$sex=as.factor(wing_df$sex)
wing_df$RingId=as.factor(wing_df$RingId)
wing_df$year=as.factor(wing_df$year)
wing_df$Observer=as.factor(wing_df$Observer)
wing_df$nestboxID=as.factor(wing_df$nestboxID)

wing_df$rank=as.numeric(wing_df$rank)
wing_df$CH1903X=as.numeric(wing_df$CH1903X) # locations of the nestboxes
wing_df$CH1903Y=as.numeric(wing_df$CH1903Y) ## fitting this means the model takes ages to run - probably because its nested within nestboxID


## read in GRM
grm=read_rds("./input_dfs/All3085_AUTOSAUMES_RP502SNPs.RDS")

grm_filt=grm[rownames(grm)%in%wing_df[["RingId"]],colnames(grm)%in%wing_df[["RingId"]]]##filtering grm for ids only in df

grm_filt_pd <- make.positive.definite(grm_filt)
GRM <- as(grm_filt_pd, "dgCMatrix")

## same as priors for simple model + prior for pe ~= 0 because we expect little var 
prior_wing=c(prior(student_t(3, 200, 90), class = "Intercept"), ## 
             prior(student_t(3, 0, 90), class = "sd"),
             prior(student_t(3, 0, 90), class = "sigma"),
             prior(cauchy(0, 5), class = "sd", group="RingId_pe"))


## slight trouble converging when using default number of itts so increased and using priors
mod_wing_GRM.Froh <- brm(LeftWing ~  1 + FHBD512gen+sex+rank+age_acc+
                           (1|gr(RingId, cov=Amat)) + (1|RingId_pe) + (1|Observer) + (1|clutch_merge) +
                           (1|year) + (1|month) + (1|nestboxID) + (1|rank),
                         data = wing_df,
                         control=list(adapt_delta=0.85),
                         data2 = list(Amat = GRM),
                         chains = 4,
                         cores=4,
                         prior=prior_wing, 
                         iter = 50000,
                         warmup = 10000,
                         thin=5
)

summary(mod_wing_GRM.Froh)

saveRDS(mod_wing_GRM.Froh,file="./outputs/1.3.ID_wing_GRM_FROH_unscaled.RDS")


