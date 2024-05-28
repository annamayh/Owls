library(brms,  lib = "/users/ahewett1/R")
library(corpcor,  lib = "/users/ahewett1/R")
library(readr,  lib = "/users/ahewett1/R")

##################################################################################
########################### ~~ Mass ~~ #############################################
#####################################################################################

args = commandArgs(trailingOnly = TRUE)
scratch = as.character(args[1]) 

mass_df=read.table("./input_dfs/mass_all_pheno_df.txt",sep=",", header=T)

mass_df$clutch_merge=as.factor(mass_df$clutch_merge)
mass_df$sex=as.factor(mass_df$sex)
mass_df$RingId=as.factor(mass_df$RingId)
mass_df$year=as.factor(mass_df$year)
mass_df$Observer=as.factor(mass_df$Observer)
mass_df$nestboxID=as.factor(mass_df$nestboxID)

mass_df$rank=as.numeric(mass_df$rank)


## read in GRM
grm=read_rds("./input_dfs/All3085_AUTOSAUMES_RP502SNPs.RDS")
grm_filt=grm[rownames(grm)%in%mass_df[["RingId"]],colnames(grm)%in%mass_df[["RingId"]]]##filtering grm for ids only in df

grm_filt_pd <- make.positive.definite(grm_filt)
GRM <- as(grm_filt_pd, "dgCMatrix")

## same as priors for simple model + prior for pe ~= 0 because we expect little var 
prior_mass=c(prior(student_t(3, 330,60), class = "Intercept"), ## 
             prior(student_t(3,0,60), class = "sd"),
             prior(cauchy(0, 10), class = "sd", group="RingId_pe"))


mod_mass_GRM.split_stage_h2 <- brm(
                      
                      bf(Mass ~  1 +sex+age_acc+ # F removed as a fixed effect 
                           (0 + gr_stage||gr(RingId, cov=Amat)) + (0 + gr_stage||RingId_pe) + (0 + gr_stage||Observer) + (0 + gr_stage||clutch_merge) +
                           (0 + gr_stage||year) + (0 + gr_stage||month) + (0 + gr_stage||nestboxID) + (0 + gr_stage||rank),
                         sigma ~ gr_stage-1), 
                      
                         data = mass_df,
                         control=list(adapt_delta=0.97),
                         data2 = list(Amat = GRM),
                         chains = 4,
                         cores=4,
                         prior=prior_mass, #
                         iter = 40000,
                         warmup = 10000,
                         thin=5
)

summary(mod_mass_GRM.split_stage_h2) ###

saveRDS(mod_mass_GRM.split_stage_h2,file=paste0(scratch,"2.2.Mass_split_stage_h2_rank_rand.RDS")) ##