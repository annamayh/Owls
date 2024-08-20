
library(tidyverse)
library(data.table)
library(brms)
library(janitor)
library(glmmTMB)

setwd("/Users/ahewett1/Documents")

hbd_segs_list <- lapply(Sys.glob("Inbreeding_depression_owls/ROH_regions/HBD_segs_out/HBD_perID_200-snp-wind_Super-Scaffold_*.RDS"), readRDS)


start_time=Sys.time()

all_gwas_out=list()

for (i in 1:length(hbd_segs_list)){
  

        hbd_segs_list_chr <- hbd_segs_list[i]%>%
        as.data.frame()%>%
        rownames_to_column(var='RingId')%>%
        select_if(~ !any(is.na(.))) # remove the columns with Nas at the end of chromosomes
      
        if(ncol(hbd_segs_list_chr)<2){
          print(paste("Skipping due to NA values"))
          next  # Skip to the next iteration
        } else{
      
      names(hbd_segs_list_chr)=make_clean_names(names(hbd_segs_list_chr))
      
      
      bill_df=read.table("Inbreeding_depression_owls/pheno_df/bill_all_pheno_df.txt",sep=",", header=T)
      
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
      
      
      all_gwas_out[[i]]=gwas_out
      
      loop_end_time=Sys.time()
      total_elapsed = difftime(loop_end_time, start_time, units = 'mins')
      print(paste0("Finished super scaffold ",i, " of ", length(hbd_segs_list) ," at: ", Sys.time()))
      print(paste0("Current elapsed time: ", round(total_elapsed, 2), " minutes"))
      }
}



unlisted_all_windows_gwas=do.call(rbind,all_gwas_out)%>%
  as.data.frame()


saveRDS(unlisted_all_windows_gwas,file=paste0("Inbreeding_depression_owls/Model_outputs/4_Gen_arch/glmTMB_GWAS_HBD_winds7500.RDS")) ##




