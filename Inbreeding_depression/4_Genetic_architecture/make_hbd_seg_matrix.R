
library(tidyverse)
library(data.table)


setwd("/Users/ahewett1/Documents")


hbd_segs_list <- lapply(Sys.glob("Inbreeding_depression_owls/ROH_regions/HBD_segs_out/HBD_perID_200-snp-wind_Super-Scaffold*.RDS"), readRDS)

hbd_segs_matrix <- map(hbd_segs_list, ~ as.data.frame(.x))%>% # convert to tibble 
  map(~ rownames_to_column(.x, var="RingId"))%>% # make column of ring id
  reduce(full_join, by="RingId") 
  


bill_df=read.table("Inbreeding_depression_owls/pheno_df/bill_all_pheno_df.txt",sep=",", header=T)%>%
  mutate(age_acc=184*exp(-0.99*0.932^age_days)) ## account for age using gompertz growth
  
## linking HBD segs to phenotype df and creating a matrix 
froh_mat=hbd_segs_matrix%>%
  right_join(bill_df, by = "RingId") %>% # to account for the repeated records
  select(-RingId) ## removing ring id to convert to df (BUT BE CAREFUL IT IS EXACTLY THE SAME AS DF)

#HBD segs as matrix for multi-memebership
bill_df$froh_mat_linked=as.matrix(froh_mat)

sum_hbd_segs=hbd_segs_matrix%>%
  mutate(hbd_sum = rowSums(.[2:ncol(hbd_segs_matrix)], na.rm = TRUE))%>% ## remove nas from last col of each chr
  mutate(hbd_sum_div=hbd_sum/(ncol(hbd_segs_matrix)-1))%>%
  select(RingId, hbd_sum, hbd_sum_div)


bill_df=bill_df%>%
  left_join(sum_hbd_segs)


plot(bill_df$hbd_sum_div, bill_df$FHBD512gen)


saveRDS(bill_df, file="Inbreeding_depression_owls/ROH_regions/df_with_HBD_mat/bill_df.RDS")
