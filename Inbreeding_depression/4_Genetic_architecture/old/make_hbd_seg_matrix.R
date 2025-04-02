library(tidyverse)
library(data.table)

setwd("/Users/ahewett1/Documents")

hbd_segs_list <- lapply(Sys.glob("Inbreeding_depression_owls/ROH_regions/HBD_segs_out/HBD_perID_7500-snp-wind_Super-Scaffold_*.RDS"), readRDS)

hbd_segs_matrix <- map(hbd_segs_list, ~ as.data.frame(.x))%>% # convert to tibble 
  map(~ rownames_to_column(.x, var="RingId"))%>% # make column of ring id
  reduce(full_join, by="RingId")%>% # reducing list into df with cols of ids
  select_if(~ !any(is.na(.))) # remove the columns with Nas at the end of chromosomes


hbd_segs_matrix[1:5,1:4]


bill_df_read=read.table("Inbreeding_depression_owls/pheno_df/bill_all_pheno_df.txt",sep=",", header=T)
  
## linking HBD segs to phenotype df and creating a matrix 
# froh_mat=hbd_segs_matrix%>%
#   right_join(bill_df_read, by = "RingId") %>% # to account for the repeated records
#   select(-RingId) ## removing ring id to convert to df (BUT BE CAREFUL IT IS EXACTLY THE SAME AS DF)
# 
# 
# sum_hbd_segs=hbd_segs_matrix%>%
#   mutate(hbd_sum = rowSums(.[2:ncol(hbd_segs_matrix)], na.rm = TRUE))%>% ## remove nas from last col of each chr
#   mutate(hbd_sum_div=hbd_sum/(ncol(hbd_segs_matrix)-1))%>%
#   select(RingId, hbd_sum, hbd_sum_div)
# 
# 
# bill_df=bill_df_read%>%
#   left_join(sum_hbd_segs)
# 
# 
# #HBD segs as matrix for multi-memebership
# bill_df$froh_mat_linked=as.matrix(froh_mat)
# 
# #plot(bill_df$hbd_sum_div, bill_df$FHBD512gen)
# 
# 
# saveRDS(bill_df, file="Inbreeding_depression_owls/ROH_regions/df_with_HBD_mat/linked_bill_df.RDS")


#######################################################################################################
######        now only for chr that have good coverage #######
##########################################################################

good_sc=read.table("Barn_owl_general/ss_goodQ.txt")
scaffold_names=good_sc$V1

scaffolds_to_keep <- paste0(scaffold_names, collapse = "|")

## linking HBD segs to phenotype df and creating a matrix 
hbd_segs_matrix_rm=hbd_segs_matrix%>%
  select(RingId, matches(scaffolds_to_keep)) ## removes about windows in bad scaffolds 


froh_mat=hbd_segs_matrix_rm%>%
    right_join(bill_df_read, by = "RingId") %>% # to account for the repeated records
    select(-RingId) %>%## removing ring id to convert to df (BUT BE CAREFUL IT IS EXACTLY THE SAME AS DF)
    data.matrix()

froh_mat[1:3, 1:3]

sum_hbd_segs=hbd_segs_matrix_rm%>%
    mutate(hbd_sum = rowSums(.[2:ncol(hbd_segs_matrix_rm)], na.rm = TRUE))%>% ## remove nas from last col of each chr
    mutate(hbd_sum_div=hbd_sum/(ncol(hbd_segs_matrix_rm)-1))%>%
    select(RingId, hbd_sum, hbd_sum_div)
  
bill_df=bill_df_read%>%
    left_join(sum_hbd_segs)%>%
    as.data.frame()

nrow(bill_df)
ncol(bill_df)
  
#HBD segs as matrix for multi-memebership
nrow(froh_mat)
bill_df$froh_mat_linked=froh_mat


any(is.na(froh_mat))

save(bill_df, file="Inbreeding_depression_owls/ROH_regions/df_with_HBD_mat/linked_bill_df_OnlyGoodSS.RDS")
