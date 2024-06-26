library(MCMCglmm,  lib = "/users/ahewett1/R")
library(readr,  lib = "/users/ahewett1/R")


args = commandArgs(trailingOnly = TRUE)
scratch = as.character(args[1])

 
load("./input_dfs/linked_bill_df_OnlyGoodSS.RData")


nrow(bill_df$froh_mat_linked)
ncol(bill_df$froh_mat_linked)


## need to play around with priors 
k<-10000
prior<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G2=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G3=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G4=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G5=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G6=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G7=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#### RUNNING MULTI-MEMBERSHIP MODEL ####
bill_arch_model<-MCMCglmm(BillLength ~ 1 +sex+age_acc+rank +hbd_sum_div,
                          random= ~ RingId + Observer + clutch_merge + year + month + nestboxID +
                            idv(froh_mat_linked), ## FROH per window matrix 
                          data=bill_df,
                          prior = prior,
                          pr=TRUE,#
                          nitt=65000,burnin=10000, thin = 10)


saveRDS(bill_arch_model,file=paste0(scratch,"4.1.bill_IDarch_model.RDS")) ##

summary(bill_arch_model)

