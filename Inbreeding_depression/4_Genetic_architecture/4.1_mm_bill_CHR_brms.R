# library(MCMCglmm,  lib = "/users/ahewett1/R")
# library(readr,  lib = "/users/ahewett1/R")
# 
# args = commandArgs(trailingOnly = TRUE)
# scratch = as.character(args[1])

library(MCMCglmm)
 
#load("./input_dfs/linked_bill_df_OnlyGoodSS.RData")

bill_df=readRDS(file="Inbreeding_depression_owls/ROH_regions/df_with_HBD_mat/CHR_bill_df_OnlyGoodSS.RDS")

head(bill_df)


chr_mm_names <- colnames(bill_df)[2:39]
print(chr_mm_names)  # Check the extracted window names

brms_formula <- as.formula(paste0("BillLength ~ 1 + sex + age_acc + rank + hbd_sum_div +
                          (1|RingId) + (1|Observer) + (1|clutch_merge) + (1|year) + 
                          (1|month) + (1|nestboxID) +
                          (1|mm(",
                                 paste(chr_mm_names, collapse = " , "), "))"))
print(brms_formula)  # Validate the constructed formula string

prior_bill_gen_arch<- c(
  prior(student_t(3, 170, 18), class="Intercept"),##
  prior(student_t(3, 0 , 18),  class="sd"), ## 
  prior(student_t(3, 0 , 18),  class="sigma")
)
#### RUNNING MULTI-MEMBERSHIP MODEL ####
bill_arch_model<-brm(brms_formula, ## FROH per window matrix 
                    data=bill_df,
                    prior = prior_bill_gen_arch,
                    cores = 3, 
                    chains = 3
                          )


saveRDS(bill_arch_model,file="Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4.1.bill_IDarch_model_CHR.RDS") ##

summary(bill_arch_model)

