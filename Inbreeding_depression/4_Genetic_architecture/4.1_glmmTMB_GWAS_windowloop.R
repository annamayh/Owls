
library(glmmTMB,lib = "/users/ahewett1/R")
library(janitor,lib = "/users/ahewett1/R")
library(dplyr,lib = "/users/ahewett1/R")
library(tibble,lib = "/users/ahewett1/R")

args <- commandArgs(trailingOnly = TRUE)

input_file=args[1]
input_dir=args[2]
scratch=args[3]

# input_dir="Inbreeding_depression_owls/ROH_regions/HBD_segs_out/"
#  input_file="HBD_perID_200-snp-wind_Super-Scaffold_15.RDS"

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
  

# cleaning the names for formatting in model
names(hbd_segs_list_chr)=make_clean_names(names(hbd_segs_list_chr))

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

bill_with_IBDinfo$mc_age_acc <- bill_with_IBDinfo$age_acc - mean(bill_with_IBDinfo$age_acc)


gwas_out=NULL
counter=0
start_time=Sys.time()

for (wind in windows){
    counter=counter+1
  
    form_wind=as.formula(paste0("BillLength ~  1 +", wind," + 
        sex + mc_age_acc + rank + FHBD512gen +
       (1|RingId_pe) + (1|Observer) + (1|clutch_merge) +
        (1|year) + (1|month) + (1|nestboxID)"))

    gwas_mod=glmmTMB(
      formula = form_wind,
      data = bill_with_IBDinfo
      
      )
    
    
    f_ests=summary(gwas_mod)$coefficients$cond[2,]
    
    f_ests=c(window = wind,f_ests)
    
    gwas_out=rbind(gwas_out, f_ests)
    
    if (as.numeric(counter) %% 10 == 0) {
      print(paste0("Finished window ", counter, " of ",length(windows)))}
  
}


saveRDS(gwas_out,file=paste0(scratch,"/glmmTMB_gwas_out_",input_file,".RDS")) ##

}

loop_end_time=Sys.time()
total_elapsed = difftime(loop_end_time, start_time, units = 'mins')
print(paste0("Finished super scaffold"))
print(paste0("Current elapsed time: ", round(total_elapsed, 2), " minutes"))


