
library(brms,lib = "/users/ahewett1/R")
library(janitor,lib = "/users/ahewett1/R")
library(dplyr,lib = "/users/ahewett1/R")
library(tibble,lib = "/users/ahewett1/R")


args <- commandArgs(trailingOnly = TRUE)

input_file=args[1]
input_dir=args[2]
scratch=args[3]


hbd_segs_list_chr <- readRDS(paste0(input_dir,"/",input_file))%>%
  as.data.frame()%>%
  rownames_to_column(var='RingId')%>%
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
bill_df=read.table("./input_dfs/bill_all_pheno_df.txt",sep=",", header=T)

bill_with_IBDinfo=hbd_segs_list_chr%>%
  rename(RingId=ring_id)%>%
  right_join(bill_df, by = 'RingId')

bill_with_IBDinfo$clutch_merge=as.factor(bill_with_IBDinfo$clutch_merge)
bill_with_IBDinfo$sex=as.factor(bill_with_IBDinfo$sex)
bill_with_IBDinfo$RingId=as.factor(bill_with_IBDinfo$RingId)
bill_with_IBDinfo$year=as.factor(bill_with_IBDinfo$year)
bill_with_IBDinfo$Observer=as.factor(bill_with_IBDinfo$Observer)
bill_with_IBDinfo$nestboxID=as.factor(bill_with_IBDinfo$nestboxID)
bill_with_IBDinfo$rank=as.numeric(bill_with_IBDinfo$rank)


windows=colnames(hbd_segs_list_chr)[2:ncol(hbd_segs_list_chr)]
# mean centering age to give a meaningful intercept
bill_with_IBDinfo$mc_age_acc <- bill_with_IBDinfo$age_acc - mean(bill_with_IBDinfo$age_acc)

prior_bill=c(prior(student_t(3, 180,15), class = "Intercept"), ## 
             prior(student_t(3,0,20), class = "sd"),
             prior(student_t(3,0,20), class = "sigma"),
             prior(cauchy(0, 5), class = "sd", group="RingId_pe"))


gwas_out=NULL
counter=0
start_time=Sys.time()

for (wind in windows){
    
     counter=counter+1
  
    form_wind=as.formula(paste0("BillLength ~  1 +", wind," + 
        sex + mc_age_acc + rank + FHBD512gen +
       (1|RingId_pe) + (1|Observer) + (1|clutch_merge) +
        (1|year) + (1|month) + (1|nestboxID)"))

    gwas_mod=brm(
      formula = form_wind,
      data = bill_with_IBDinfo,
      chains = 3,
      cores=3,
      prior=prior_bill ## default itts      
      )
    
    
    f_ests=fixef(gwas_mod, pars = paste0(wind))
    
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

