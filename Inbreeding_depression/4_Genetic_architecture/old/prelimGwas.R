
library(tidyverse)
library(data.table)
library(brms)
library(janitor)
library(glmmTMB)

setwd("/Users/ahewett1/Documents")

# hbd_segs_list <- lapply(Sys.glob("Inbreeding_depression_owls/ROH_regions/HBD_segs_out/HBD_perID_200-snp-wind_Super-Scaffold_*.RDS"), readRDS)
# 
# hbd_segs_matrix <- map(hbd_segs_list, ~ as.data.frame(.x))%>% # convert to tibble 
#   map(~ rownames_to_column(.x, var="RingId"))%>% # make column of ring id
#   reduce(full_join, by="RingId")%>% # reducing list into df with cols of ids
#   select_if(~ !any(is.na(.))) # remove the columns with Nas at the end of chromosomes

hbd_segs_list_chr <- readRDS("Inbreeding_depression_owls/ROH_regions/HBD_segs_out/HBD_perID_7500-snp-wind_Super-Scaffold_14.RDS")%>%
  as.data.frame()%>%
  rownames_to_column(var='RingId')%>%
  select_if(~ !any(is.na(.))) # remove the columns with Nas at the end of chromosomes
  

names(hbd_segs_list_chr)=make_clean_names(names(hbd_segs_list_chr))


# ## filtering out scaffolds with bad qualti
good_sc=read.table("Barn_owl_general/ss_goodQ.txt")
scaffold_names=good_sc$V1


scaffolds_to_keep <- paste0(scaffold_names, collapse = "|")

## linking HBD segs to phenotype df and creating a matrix
hbd_segs_matrix_rm=hbd_segs_matrix%>%
  select(RingId, matches(scaffolds_to_keep)) # removes about windows in bad scaffolds

  


bill_df_read=read.table("Inbreeding_depression_owls/pheno_df/bill_all_pheno_df.txt",sep=",", header=T)

bill_with_IBDinfo=hbd_segs_list_chr%>%
  rename(RingId=ring_id)%>%
  right_join(bill_df_read, by = 'RingId')

bill_with_IBDinfo$clutch_merge=as.factor(bill_with_IBDinfo$clutch_merge)
bill_with_IBDinfo$sex=as.factor(bill_with_IBDinfo$sex)
bill_with_IBDinfo$RingId=as.factor(bill_with_IBDinfo$RingId)
bill_with_IBDinfo$year=as.factor(bill_with_IBDinfo$year)
bill_with_IBDinfo$Observer=as.factor(bill_with_IBDinfo$Observer)
bill_with_IBDinfo$nestboxID=as.factor(bill_with_IBDinfo$nestboxID)
bill_with_IBDinfo$rank=as.numeric(bill_with_IBDinfo$rank)


windows=colnames(hbd_segs_list_chr)[2:ncol(hbd_segs_list_chr)]

bill_with_IBDinfo$mc_age_acc <- bill_with_IBDinfo$age_acc - mean(bill_with_IBDinfo$age_acc)


prior_bill=c(prior(student_t(3, 180,20), class = "Intercept"), ##
             prior(student_t(3,0,20), class = "sd"),
             prior(student_t(3,0,20), class = "sigma"),
             prior(cauchy(0, 5), class = "sd", group="RingId_pe"))


wind=windows[1]

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
      prior=prior_bill, ##
      iter = 5000,
      warmup = 2000,
      thin=5
    )
    
    f_ests=fixef(gwas_mod, pars = paste0(wind))
    
    gwas_out=rbind(gwas_out, f_ests)
    
    if (as.numeric(counter) %% 10 == 0) {
      elapsed_time = difftime(Sys.time(), start_time, units = 'mins')
      print(paste0("Finished window ", counter, " of ",length(windows)))
      print(paste0("Elapsed time: ", round(elapsed_time, 2)))
    }

}
