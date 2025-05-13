.libPaths(c("/work/FAC/FBM/DEE/jgoudet/barn_owl/ahewett/R", .libPaths())) #specify library cluster path

library(brms)
library(janitor)
library(dplyr)
library(tibble)


args <- commandArgs(trailingOnly = TRUE)

input_file=args[1]
input_dir=args[2]
scratch=args[3]


hbd_segs_list_chr <- readRDS(paste0(input_dir,"/",input_file))%>%
  as.data.frame()%>%
  #rownames_to_column(var='RingId')%>%
  select_if(~ !any(is.na(.))) # remove the columns with Nas at the end of chromosomes

## This is to skip the chromosomes that have bad quality
if(ncol(hbd_segs_list_chr)<3){
  print(paste0("Exiting due to NA values"))
  stop()  # Stop program
} else{
  
  start_time=Sys.time() #starting this chromosome
  print(paste0("Starting super scaffold ",input_file," at: ", start_time))
  
  
names(hbd_segs_list_chr)=make_clean_names(names(hbd_segs_list_chr)) ## removing weird syntax for pasting into model later

# reading in the 
mass_df=read.table("./input_dfs/mass_all_pheno_df.txt",sep=",", header=T)


mass_with_IBDinfo=hbd_segs_list_chr%>%
  rename(RingId=ring_id)%>%
  right_join(mass_df, by = 'RingId')

mass_with_IBDinfo$clutch_merge=as.factor(mass_with_IBDinfo$clutch_merge)
mass_with_IBDinfo$sex=as.factor(mass_with_IBDinfo$sex)
mass_with_IBDinfo$RingId=as.factor(mass_with_IBDinfo$RingId)
mass_with_IBDinfo$year=as.factor(mass_with_IBDinfo$year)
mass_with_IBDinfo$Observer=as.factor(mass_with_IBDinfo$Observer)
mass_with_IBDinfo$nestboxID=as.factor(mass_with_IBDinfo$nestboxID)
mass_with_IBDinfo$rank=as.numeric(mass_with_IBDinfo$rank)


windows=colnames(hbd_segs_list_chr)[2:ncol(hbd_segs_list_chr)]
# mean centering age to give a meaningful intercept
mass_with_IBDinfo$mc_age_acc <- mass_with_IBDinfo$age_acc - mean(mass_with_IBDinfo$age_acc)

prior_mass=c(prior(student_t(3, 330, 60), class = "Intercept"), ## 
             prior(student_t(3, 0, 60), class = "sd"),
             prior(student_t(3, 0, 60), class = "sigma"),
             prior(cauchy(0, 5), class = "sd", group="RingId_pe"))


gwas_out=NULL
counter=0
start_time=Sys.time()

for (wind in windows){
    
     counter=counter+1
  
    form_wind=as.formula(paste0("Mass ~  1 +", wind," + 
        sex + mc_age_acc + rank + FROH +
       (1|RingId_pe) + (1|Observer) + (1|clutch_merge) +
        (1|year) + (1|month) + (1|nestboxID)"))

    gwas_mod=brm(
      formula = form_wind,
      data = mass_with_IBDinfo,
      chains = 3,
      cores=3,
      prior=prior_mass ## default itts      
      )
    
    
    # output estimate of window for each chain
    f_ests_all_chains=fixef(gwas_mod, pars = paste0(wind), summary = F)
    
    blw_zero=length(which(sign(f_ests_all_chains)%in%'-1'))## number of times estimate is below 0
    abv_zero=length(which(sign(f_ests_all_chains)%in%'1'))## number of times estimate is above 0
    
    total=nrow(f_ests_all_chains) ## total number of chains 
    eps <- 1e-10 # very small number so that we arent dividing by 0
    
    # estimate 
    f_est_sign=fixef(gwas_mod, pars = paste0(wind))[1]
    
    # Calulate P-value = the proportion of the MCMC chain that is the opposite sign to the point estimate. 
    if (sign(f_est_sign) < 0 ){
      p_val=max(abv_zero, eps)/total
    }
    if (sign(f_est_sign) > 0 ){
      p_val=max(blw_zero,eps)/total
    }
    
    ## combine into df
    f_ests <- data.frame(
      Window = wind,
      estimate = fixef(gwas_mod, pars = paste0(wind))[[1]],
      lwr_CI = fixef(gwas_mod, pars = paste0(wind))[[3]],
      upr_CI = fixef(gwas_mod, pars = paste0(wind))[[4]],
      p_val = p_val
    )
    # merge with previous window info
    gwas_out=rbind(gwas_out, f_ests)
    
    if (as.numeric(counter) %% 10 == 0) {
      print(paste0(">>> FINISHED WINDOW ", counter, " OF ",length(windows)))
    }

}


saveRDS(gwas_out,file=paste0(scratch,"gwas_out_",input_file,".RDS")) ##
}


end_time=Sys.time()
total_elapsed = difftime(end_time, start_time, units = 'mins')
print(paste0("Finished at: ", Sys.time()))
print(paste0("Total elapsed time: ", round(total_elapsed, 2), " minutes"))
